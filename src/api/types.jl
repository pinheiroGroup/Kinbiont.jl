# =============================================================================
# Unified Matrix API — Type Definitions
# =============================================================================
# Three public structs that form the entire API surface:
#   GrowthData   — the data (immutable, matrix-centric)
#   FitOptions   — all configuration with sensible defaults (@kwdef)
#   ModelSpec    — which models to fit and their initial parameters
#
# Plus the model trait types (NLModel, ODEModel, LogLinModel) used for
# dispatch, and the result types returned by fit() and preprocess().
# =============================================================================

using OptimizationBBO: BBO_adaptive_de_rand_1_bin_radiuslimited
using OrdinaryDiffEq: Tsit5

# ---------------------------------------------------------------------------
# 1. Data container
# ---------------------------------------------------------------------------

"""
    GrowthData(curves, times, labels[, clusters[, centroids]])

Immutable container for a set of growth curves measured at common time points.

# Fields
- `curves::Matrix{Float64}`: `n_curves × n_timepoints` matrix; each row is one curve.
- `times::Vector{Float64}`: time points shared by all curves (length `n_timepoints`).
- `labels::Vector{String}`: identifier for each curve (length `n_curves`).
- `clusters::Union{Nothing, Vector{Int}}`: cluster assignment per curve (1-based, within
  `1..n_clusters`), populated by [`preprocess`](@ref). `nothing` until clustering is run.
- `centroids::Union{Nothing, Matrix{Float64}}`: `n_clusters × n_timepoints` matrix of
  per-cluster shape centroids in **z-normalised space**, populated alongside `clusters`.
  Row `k` is the mean of the z-scored curves assigned to cluster `k` — a scale-independent
  shape prototype. To obtain original-space prototypes, compute
  `mean(data.curves[data.clusters .== k, :], dims=1)` for each `k`.
- `wcss::Union{Nothing, Float64}`: within-cluster sum of squares (total SSE) returned
  by k-means, populated alongside `clusters`. Use this to construct an elbow plot and
  choose `n_clusters`. `nothing` when clustering was not performed.
"""
struct GrowthData
    curves::Matrix{Float64}
    times::Vector{Float64}
    labels::Vector{String}
    clusters::Union{Nothing, Vector{Int}}
    centroids::Union{Nothing, Matrix{Float64}}
    wcss::Union{Nothing, Float64}

    function GrowthData(curves, times, labels, clusters=nothing, centroids=nothing, wcss=nothing)
        n_curves, n_tp = size(curves)
        length(times)  == n_tp     || error("times length $(length(times)) ≠ n_timepoints $n_tp")
        length(labels) == n_curves || error("labels length $(length(labels)) ≠ n_curves $n_curves")
        clusters === nothing || length(clusters) == n_curves ||
            error("clusters length $(length(clusters)) ≠ n_curves $n_curves")
        new(curves, times, labels, clusters, centroids, wcss)
    end
end

"""
    GrowthData(path::String)

Load growth curves from a CSV file and return a `GrowthData`.

The file must follow the Kinbiont column convention:
- **First column**: time points (numeric).
- **Remaining columns**: one curve per column; the column header becomes the curve label.

# Example
```julia
data = GrowthData("/path/to/experiment.csv")
```
"""
function GrowthData(path::String)
    tbl    = CSV.File(path)
    cols   = propertynames(tbl)
    times  = Float64.(tbl[cols[1]])
    labels = String.(cols[2:end])
    curves = Matrix{Float64}(reduce(hcat, [Float64.(tbl[c]) for c in cols[2:end]])')
    return GrowthData(curves, times, labels)
end

# ---------------------------------------------------------------------------
# 2. All configuration in one place, every field has a default
# ---------------------------------------------------------------------------

"""
    FitOptions(; kwargs...)

All configuration for preprocessing and fitting in a single struct.
Every field has a sensible default so users only override what they need.

# Preprocessing fields
- `smooth::Bool = false`: apply smoothing before fitting.
- `smooth_method::Symbol = :lowess`: `:lowess`, `:rolling_avg`, `:gaussian`, or `:none`.
- `smooth_pt_avg::Int = 7`: window size for `:rolling_avg`.
- `lowess_frac::Float64 = 0.05`: bandwidth fraction for `:lowess`.
- `gaussian_h_mult::Float64 = 2.0`: bandwidth multiplier for Gaussian smoothing
  (bandwidth = `gaussian_h_mult × median(Δt)`).
- `gaussian_time_grid::Union{Nothing,Vector{Float64}} = nothing`: optional target
  time grid for Gaussian smoothing; when set, smoothed curves are evaluated at
  these times (interpolation). `nothing` keeps the original time grid.
- `blank_subtraction::Bool = false`: subtract a blank value from all curves.
- `blank_value::Float64 = 0.0`: constant blank to subtract when `blank_subtraction=true`.
- `correct_negatives::Bool = false`: handle negative values after blank subtraction.
- `negative_method::Symbol = :remove`: `:remove`, `:thr_correction`, or `:blank_correction`.
- `negative_threshold::Float64 = 0.01`: floor value used by `:thr_correction`.
- `scattering_correction::Bool = false`: apply multiple-scattering OD correction.
- `calibration_file::String = ""`: path to calibration CSV (required when `scattering_correction=true`).
- `scattering_method::Symbol = :interpolation`: `:interpolation` or `:exp_fit`.

# Stationary phase fields
- `cut_stationary_phase::Bool = false`: truncate each curve at the onset of stationary
  phase before fitting. Uses a specific-growth-rate threshold to locate the cutoff,
  then snaps it to the OD peak within the next `stationary_win_size` time points.
- `stationary_percentile_thr::Float64 = 0.05`: SGR drops below this fraction of its
  maximum to trigger the cutoff.
- `stationary_pt_smooth_derivative::Int = 10`: smoothing window for SGR evaluation.
- `stationary_win_size::Int = 5`: look-ahead window for OD-peak snapping.
- `stationary_thr_od::Float64 = 0.02`: minimum OD value; time points below this are
  excluded before detecting stationary phase.

# Clustering fields
- `cluster::Bool = false`: cluster curves after preprocessing.
- `n_clusters::Int = 3`: total number of cluster labels (1..`n_clusters`). Labels are
  always within this range regardless of other options.
- `cluster_trend_test::Bool = true`: reserve label `n_clusters` for flat/non-growing
  curves (identified by an OLS slope t-test, p ≥ 0.05). K-means then runs with
  `n_clusters - 1` dynamic groups, so all labels remain in `1..n_clusters`.
  Requires `n_clusters ≥ 2`. Ignored when `cluster_prescreen_constant=true`.
- `cluster_prescreen_constant::Bool = false`: before running k-means, identify
  non-growing wells using a quantile-ratio criterion (high tail / low tail ≤
  `cluster_tol_const`) and pin them to label `n_clusters`. K-means then runs with
  `n_clusters - 1` groups on the remaining dynamic curves only. More biologically
  meaningful than the post-hoc `cluster_trend_test` because k-means never sees
  flat wells. Requires `n_clusters ≥ 2`.
- `cluster_tol_const::Float64 = 1.5`: threshold for constant pre-screening.
  A curve is flagged as non-growing when its `cluster_q_high` quantile is ≤
  `cluster_tol_const × cluster_q_low` quantile. Increase to be more permissive
  (fewer constant calls); decrease to be stricter.
- `cluster_q_low::Float64 = 0.05`: lower quantile used to estimate the baseline
  tail in constant pre-screening.
- `cluster_q_high::Float64 = 0.95`: upper quantile used to estimate the signal
  tail in constant pre-screening.
- `cluster_exp_prototype::Bool = false`: add a dedicated "exponential" cluster
  (label `n_clusters - 1` when used together with constant pre-screening, or
  `n_clusters` when used alone). For each curve, the distance to the nearest
  exponential shape prototype (z-scored, bases 2⁶..2¹⁶) is compared to the
  distance to the nearest k-means centroid; the exponential label wins when
  it is closer. Requires `n_clusters ≥ 3` when combined with constant pre-screening
  (`n_clusters - 1` for exponential, `n_clusters` for constant), or ≥ 2 otherwise.

After clustering, `processed.wcss` holds the within-cluster sum of squares. Run
`preprocess` for `n_clusters = 2, 3, 4, ...` and plot `wcss` vs `n_clusters` to
find the elbow and choose the optimal number of clusters.

# Fitting fields
- `loss::String = "RE"`: loss function — `"RE"`, `"L2"`, `"L2_derivative"`, etc.
- `optimizer`: BBO optimizer instance (default: `BBO_adaptive_de_rand_1_bin_radiuslimited()`).
- `integrator`: ODE integrator (default: `Tsit5()`).
- `multistart::Bool = false`: enable multistart optimization.
- `n_restart::Int = 50`: number of restarts when `multistart=true`.
- `aic_correction::Bool = true`: use AICc (corrected AIC) for model selection.
- `pt_smooth_derivative::Int = 7`: window for specific-growth-rate evaluation.
- `auto_diff_method`: autodiff backend passed to Optimization.jl (`nothing` = no autodiff).
- `cons`: constraint function for Optimization.jl (`nothing` = unconstrained).
- `opt_params::NamedTuple`: extra keyword arguments forwarded verbatim to
  `OptimizationProblem` / `solve` — e.g. `(maxiters=1_000_000, abstol=1e-9)`.
  Accepts anything the underlying Optimization.jl solver understands.
"""
@kwdef struct FitOptions
    # --- preprocessing ---
    smooth::Bool                = false
    smooth_method::Symbol       = :lowess
    smooth_pt_avg::Int          = 7
    lowess_frac::Float64        = 0.05
    gaussian_h_mult::Float64    = 2.0
    gaussian_time_grid::Union{Nothing, Vector{Float64}} = nothing
    blank_subtraction::Bool     = false
    blank_value::Float64        = 0.0
    correct_negatives::Bool     = false
    negative_method::Symbol     = :remove
    negative_threshold::Float64 = 0.01
    scattering_correction::Bool = false
    calibration_file::String    = ""
    scattering_method::Symbol   = :interpolation

    # --- stationary phase detection ---
    cut_stationary_phase::Bool             = false
    stationary_percentile_thr::Float64     = 0.05
    stationary_pt_smooth_derivative::Int   = 10
    stationary_win_size::Int               = 5
    stationary_thr_od::Float64             = 0.02

    # --- clustering ---
    cluster::Bool                    = false
    n_clusters::Int                  = 3
    cluster_trend_test::Bool         = true
    cluster_prescreen_constant::Bool = false
    cluster_tol_const::Float64       = 1.5
    cluster_q_low::Float64           = 0.05
    cluster_q_high::Float64          = 0.95
    cluster_exp_prototype::Bool      = false

    # --- fitting ---
    loss::String                = "RE"
    optimizer                   = BBO_adaptive_de_rand_1_bin_radiuslimited()
    integrator                  = Tsit5()
    multistart::Bool            = false
    n_restart::Int              = 50
    aic_correction::Bool        = true
    pt_smooth_derivative::Int   = 7
    auto_diff_method            = nothing
    cons                        = nothing
    opt_params::NamedTuple      = (;)
end

# ---------------------------------------------------------------------------
# 3. Model trait types — dispatch replaces string matching
# ---------------------------------------------------------------------------

"""
Abstract base for all growth models. Concrete subtypes: [`NLModel`](@ref),
[`ODEModel`](@ref), [`LogLinModel`](@ref).
"""
abstract type AbstractGrowthModel end

"""
    NLModel(name, func, param_names[, guess])

A non-linear (closed-form) growth model.

# Fields
- `name::String`: unique identifier, used in results.
- `func::Function`: model function with signature `(p, t) -> y` or `(p, times) -> ys`.
- `param_names::Vector{String}`: human-readable name for each parameter.
- `guess::Union{Function,Nothing}`: `(data::Matrix{Float64}) -> Vector{Float64}` returning
  a starting guess; `nothing` means the guess must be supplied in [`ModelSpec`](@ref).
"""
struct NLModel <: AbstractGrowthModel
    name::String
    func::Function
    param_names::Vector{String}
    guess::Union{Function, Nothing}
end

NLModel(name, func, param_names) = NLModel(name, func, param_names, nothing)

"""
    ODEModel(name, func, param_names, n_eq[, guess])

An ODE growth model in SciML in-place form `f!(du, u, p, t)`.

# Fields
- `name::String`: unique identifier, used in results.
- `func::Function`: ODE function with signature `f!(du, u, p, t)`.
- `param_names::Vector{String}`: human-readable name for each parameter.
- `n_eq::Int`: number of state equations (replaces the error-prone detection loop).
- `guess::Union{Function,Nothing}`: `(data::Matrix{Float64}) -> Vector{Float64}` returning
  a starting guess; `nothing` means the guess must be supplied in [`ModelSpec`](@ref).
"""
struct ODEModel <: AbstractGrowthModel
    name::String
    func::Function
    param_names::Vector{String}
    n_eq::Int
    guess::Union{Function, Nothing}
end

ODEModel(name, func, param_names, n_eq) = ODEModel(name, func, param_names, n_eq, nothing)

"""
    LogLinModel()

Sentinel type for log-linear (exponential phase) fitting. No parameters or
model function required — the fitting procedure is fully determined by
[`FitOptions`](@ref) fields `pt_smooth_derivative`, etc.
"""
struct LogLinModel <: AbstractGrowthModel end

"""
    DDDEModel(; max_degree=4, lambda_min=-5.0, lambda_max=-1.0, lambda_step=0.5)

Data-Driven Differential Equation (DDDE) discovery model. Uses sparse regression
(STLSQ) on a polynomial basis to identify the governing ODE directly from the
data, without assuming a fixed model form.

!!! note "Optional dependency"
    Fitting with `DDDEModel` requires `DataDrivenDiffEq`, `DataDrivenSparse`, and
    `ModelingToolkit` to be installed and loaded in the active Julia environment.
    These are not listed as Kinbiont dependencies — add them with `Pkg.add(...)`.

# Fields
- `max_degree::Int = 4`: maximum polynomial degree in the candidate basis.
- `lambda_min::Float64 = -5.0`: log₁₀ of the minimum STLSQ sparsity threshold.
- `lambda_max::Float64 = -1.0`: log₁₀ of the maximum STLSQ sparsity threshold.
- `lambda_step::Float64 = 0.5`: step in log₁₀ space between threshold values.
"""
struct DDDEModel <: AbstractGrowthModel
    max_degree::Int
    lambda_min::Float64
    lambda_max::Float64
    lambda_step::Float64
end

DDDEModel(;
    max_degree  = 4,
    lambda_min  = -5.0,
    lambda_max  = -1.0,
    lambda_step = 0.5,
) = DDDEModel(max_degree, lambda_min, lambda_max, lambda_step)

# ---------------------------------------------------------------------------
# 4. Model specification passed to fit()
# ---------------------------------------------------------------------------

"""
    ModelSpec(models, params[; lower, upper])

Specifies which models to try and their initial parameters.
[`fit`](@ref) evaluates every model and selects the best by AICc.

# Fields
- `models::Vector{<:AbstractGrowthModel}`: candidate models (from [`MODEL_REGISTRY`](@ref)
  or custom).
- `params::Vector{Vector{Float64}}`: initial parameter guess per model.
  Use an empty vector `Float64[]` for [`LogLinModel`](@ref).
- `lower::Union{Nothing,Vector{Union{Nothing,Vector{Float64}}}}`: per-model lower bounds;
  `nothing` or a slot set to `nothing` means no lower bound for that model.
- `upper::Union{Nothing,Vector{Union{Nothing,Vector{Float64}}}}`: per-model upper bounds.

# Example
```julia
spec = ModelSpec(
    [MODEL_REGISTRY["logistic"], MODEL_REGISTRY["gompertz"]],
    [[0.01, 0.5, 1.2],           [0.01, 0.5, 1.2, 0.5]],
)
```
"""
struct ModelSpec
    models::Vector{<:AbstractGrowthModel}
    params::Vector{Vector{Float64}}
    lower::Union{Nothing, Vector{Union{Nothing, Vector{Float64}}}}
    upper::Union{Nothing, Vector{Union{Nothing, Vector{Float64}}}}

    function ModelSpec(models, params; lower=nothing, upper=nothing)
        length(models) == length(params) ||
            error("models and params must have the same length")
        lower === nothing || length(lower) == length(models) ||
            error("lower bounds must have the same length as models")
        upper === nothing || length(upper) == length(models) ||
            error("upper bounds must have the same length as models")
        new(models, params, lower, upper)
    end
end

# Positional convenience constructor (no keyword args)
ModelSpec(models, params, lower, upper) = ModelSpec(models, params; lower=lower, upper=upper)

# ---------------------------------------------------------------------------
# 5. Result types
# ---------------------------------------------------------------------------

"""
    CurveFitResult

Fitting result for a single growth curve.

# Fields
- `label::String`: curve identifier from [`GrowthData`](@ref).
- `best_model::AbstractGrowthModel`: the model selected by AICc.
- `best_params::Vector{Any}`: fitted parameter vector for the best model.
- `best_aic::Float64`: AICc of the best model.
- `fitted_curve::Vector{Float64}`: model values at each time point.
- `times::Vector{Float64}`: time points (same as input data).
- `loss::Float64`: final objective value.
- `all_results::Vector{NamedTuple}`: raw result for every candidate model tried.
"""
struct CurveFitResult
    label::String
    best_model::AbstractGrowthModel
    best_params::Vector{Any}
    best_aic::Float64
    fitted_curve::Vector{Float64}
    times::Vector{Float64}
    loss::Float64
    all_results::Vector{<:NamedTuple}
end

"""
    GrowthFitResults

Top-level result returned by [`fit`](@ref).

# Fields
- `data::GrowthData`: the (possibly preprocessed) data that was fitted.
- `results::Vector{CurveFitResult}`: one entry per curve.
"""
struct GrowthFitResults
    data::GrowthData
    results::Vector{CurveFitResult}
end

Base.length(r::GrowthFitResults) = length(r.results)
Base.iterate(r::GrowthFitResults, args...) = iterate(r.results, args...)
Base.getindex(r::GrowthFitResults, i) = r.results[i]

export GrowthData
export FitOptions
export AbstractGrowthModel, NLModel, ODEModel, LogLinModel, DDDEModel
export ModelSpec
export CurveFitResult, GrowthFitResults
