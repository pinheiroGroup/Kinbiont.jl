# =============================================================================
# Model Recommendation via Feature Fingerprinting
# =============================================================================
# Public API:
#   build_model_fingerprint_db(; ...) -> ModelFingerprintDB
#   recommend_models(curve, times, db; top_k, method) -> Vector{String}
#   suggest_p0(curve, times, db, model_name) -> Vector{Float64}
#   smart_fit(data, db[, opts]; top_k, method) -> GrowthFitResults
# =============================================================================

using Random: AbstractRNG, GLOBAL_RNG, randperm
using Statistics: var, mean

_norm2(v) = sqrt(sum(x^2 for x in v))

# Number of time points used to store each normalised curve (`:curves` method).
const N_CURVE_PTS = 50

# ---------------------------------------------------------------------------
# Result struct
# ---------------------------------------------------------------------------

"""
    ModelFingerprintDB

Database of simulated growth-curve features used for model recommendation.
Built by [`build_model_fingerprint_db`](@ref) and consumed by
[`recommend_models`](@ref) / [`smart_fit`](@ref).

# Fields
- `model_names::Vector{String}`: model name for each row (one row per successful simulation).
- `model_objects::Vector{AbstractGrowthModel}`: corresponding model objects.
- `features::Matrix{Float64}`: `n_rows × 15` L2-normalised feature matrix (`:features` method).
- `feature_names::Vector{String}`: human-readable name for each of the 15 features.
- `params::Vector{Vector{Float64}}`: parameter vector used to generate each row.
  Used by [`suggest_p0`](@ref) to seed the optimiser with a nearest-neighbour initial guess.
- `curve_matrix::Matrix{Float64}`: `n_rows × N_CURVE_PTS` matrix of Pearson-ready curves
  (`:curves` method). Each row is the simulated curve min-max normalised to [0, 1],
  interpolated to a common unit time grid, then mean-centred and L2-normalised so that
  a dot product with a query vector equals the Pearson correlation.
"""
struct ModelFingerprintDB
    model_names  ::Vector{String}
    model_objects::Vector{AbstractGrowthModel}
    features     ::Matrix{Float64}   # n_rows × 15, L2-normalised
    feature_names::Vector{String}
    params       ::Vector{Vector{Float64}}   # parameter vector per row
    curve_matrix ::Matrix{Float64}   # n_rows × N_CURVE_PTS, Pearson-ready
end

const _FEATURE_NAMES = [
    "log_dynamic_range",
    "log_od0",
    "max_spec_gr",
    "t_max_gr_norm",
    "lag_est_norm",
    "stationary_est_norm",
    "auc_norm",
    "gr_variance",
    "gr_skewness",
    "late_slope",
    "od_end_ratio",
    "curvature_at_inflection",
    "early_auc_frac",
    "n_sign_changes_gr",
    "exp_phase_fraction",
]

# ---------------------------------------------------------------------------
# Feature extraction  (`:features` method)
# ---------------------------------------------------------------------------

function _trapz(t::Vector{Float64}, y::Vector{Float64})::Float64
    s = 0.0
    for i in 2:length(t)
        s += (t[i] - t[i-1]) * (y[i] + y[i-1])
    end
    return 0.5 * s
end

"""
    _extract_features(t, y) -> Vector{Float64}

Compute the 15-element shape-feature vector from a growth curve.
Undefined features (e.g. lag when growth never reaches 10%) are set to 0.
The returned vector is **not** L2-normalised (normalisation happens in the DB builder).
"""
function _extract_features(t::Vector{Float64}, y::Vector{Float64})::Vector{Float64}
    n = length(t)
    if n < 2
        return zeros(15)
    end

    tmax  = t[end] - t[1]
    ymax  = maximum(y)
    ymin  = minimum(y)

    logy = log.(max.(y, 1e-12))
    spec_gr = Vector{Float64}(undef, n - 1)
    for i in 1:(n-1)
        dt = t[i+1] - t[i]
        spec_gr[i] = dt > 0 ? (logy[i+1] - logy[i]) / dt : 0.0
    end

    if isempty(spec_gr)
        return zeros(15)
    end

    max_sgr = maximum(spec_gr)
    igr     = argmax(spec_gr)
    igr_t   = t[igr]

    f1 = log10(max(ymax / max(ymin, 1e-9), 1.0))
    f2 = log10(max(y[1], 1e-9))
    f3 = max_sgr
    f4 = tmax > 0 ? (igr_t - t[1]) / tmax : 0.0

    thresh_lag = y[1] + 0.1 * (ymax - y[1])
    idx_lag    = findfirst(yi -> yi >= thresh_lag, y)
    f5 = (idx_lag !== nothing && tmax > 0) ? (t[idx_lag] - t[1]) / tmax : 0.0

    thresh_stat = y[1] + 0.9 * (ymax - y[1])
    idx_stat    = findfirst(yi -> yi >= thresh_stat, y)
    f6 = (idx_stat !== nothing && tmax > 0) ? (t[idx_stat] - t[1]) / tmax : 0.0

    total_auc = _trapz(t, y)
    f7 = (tmax > 0 && ymax > 0) ? total_auc / (tmax * ymax) : 0.0

    f8 = length(spec_gr) > 1 ? var(spec_gr) : 0.0

    f9 = if length(spec_gr) >= 3
        mu  = mean(spec_gr)
        sig = sqrt(max(var(spec_gr), 1e-30))
        mean((s - mu)^3 for s in spec_gr) / sig^3
    else
        0.0
    end

    n_late = max(2, round(Int, 0.2 * n))
    t_late = t[(end - n_late + 1):end]
    y_late = y[(end - n_late + 1):end]
    tm     = mean(t_late)
    ym     = mean(y_late)
    denom  = sum((ti - tm)^2 for ti in t_late)
    f10    = denom > 0 ? sum((t_late[i] - tm) * (y_late[i] - ym) for i in eachindex(t_late)) / denom : 0.0

    f11 = ymax > 0 ? y[end] / ymax : 0.0

    yi  = min(igr + 1, n)
    f12 = if 2 <= yi <= n-1
        (y[yi+1] - 2*y[yi] + y[yi-1]) / max((t[yi] - t[yi-1])^2, 1e-12)
    else
        0.0
    end

    mid_t     = t[1] + 0.5 * tmax
    idx_mid   = searchsortedlast(t, mid_t)
    idx_mid   = clamp(idx_mid, 2, n)
    early_auc = _trapz(t[1:idx_mid], y[1:idx_mid])
    f13 = total_auc > 1e-12 ? early_auc / total_auc : 0.0

    f14 = Float64(sum(spec_gr[i] * spec_gr[i+1] < 0 for i in 1:(length(spec_gr)-1)))
    f15 = max_sgr > 0 ? count(s -> s > 0.5 * max_sgr, spec_gr) / length(spec_gr) : 0.0

    feat = [f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, f14, f15]
    return map(x -> isfinite(x) ? x : 0.0, feat)
end

# ---------------------------------------------------------------------------
# Direct-curve helpers  (`:curves` method)
# ---------------------------------------------------------------------------

# Linear interpolation on a sorted source grid.
function _interp_linear(
    t_src::Vector{Float64},
    y_src::Vector{Float64},
    t_new::Vector{Float64},
)::Vector{Float64}
    y_new = Vector{Float64}(undef, length(t_new))
    for (k, t) in enumerate(t_new)
        i = searchsortedlast(t_src, t)
        if i == 0
            y_new[k] = y_src[1]
        elseif i >= length(t_src)
            y_new[k] = y_src[end]
        else
            α = (t - t_src[i]) / (t_src[i+1] - t_src[i])
            y_new[k] = y_src[i] + α * (y_src[i+1] - y_src[i])
        end
    end
    return y_new
end

"""
    _pearson_ready_curve(t_sim, y, n_pts) -> Vector{Float64}

Convert a simulated (or observed) curve to a Pearson-correlation-ready vector:

1. Normalise time to [0, 1].
2. Interpolate `y` to `n_pts` evenly-spaced points on the unit grid.
3. Min-max normalise the interpolated values to [0, 1]  →  **scale-invariant**.
4. Mean-centre and L2-normalise.

The resulting unit vector is scale-invariant (same for OD or fluorescence).
A dot product of two such vectors equals their Pearson correlation.
Constant curves return a zero vector.
"""
function _pearson_ready_curve(
    t_sim::Vector{Float64},
    y::Vector{Float64},
    n_pts::Int,
)::Vector{Float64}
    t_span  = max(t_sim[end] - t_sim[1], 1e-12)
    t_norm  = (t_sim .- t_sim[1]) ./ t_span
    t_grid  = collect(range(0.0, 1.0; length=n_pts))
    y_interp = _interp_linear(t_norm, y, t_grid)

    lo, hi = minimum(y_interp), maximum(y_interp)
    y_mm   = (y_interp .- lo) ./ max(hi - lo, 1e-12)

    mu = mean(y_mm)
    z  = y_mm .- mu
    zn = _norm2(z)
    return zn > 1e-12 ? z ./ zn : z
end

# ---------------------------------------------------------------------------
# Latin Hypercube Sampling
# ---------------------------------------------------------------------------

function _lhs_sample(
    lb::Vector{Float64},
    ub::Vector{Float64},
    n::Int;
    rng::AbstractRNG = GLOBAL_RNG,
)::Matrix{Float64}
    d = length(lb)
    X = Matrix{Float64}(undef, n, d)
    for j in 1:d
        perm    = randperm(rng, n)
        X[:, j] = lb[j] .+ (ub[j] - lb[j]) .* ((perm .- rand(rng, n)) ./ n)
    end
    return X
end

# ---------------------------------------------------------------------------
# Parameter bounds
# ---------------------------------------------------------------------------

# Synthetic logistic curve for NL model guess, scaled to measurement units.
function _synthetic_logistic_data(tmax::Float64, n_points::Int, od_scale::Float64 = 1.0)::Matrix{Float64}
    t  = range(0.0, tmax; length=n_points)
    N0 = 0.03 * od_scale; K = 1.2 * od_scale; r = 0.3
    y  = @. K / (1 + (K/N0 - 1) * exp(-r * t))
    return Matrix(transpose(hcat(collect(t), y)))
end

const _ODE_PARAM_BOUNDS = [
    # (substring, lb, ub) — case-insensitive; carrying-capacity entries scaled separately
    ("gr",                   0.01, 3.0),
    ("n_max",                0.1,  5.0),
    ("n_res",                0.1,  5.0),
    ("lag",                  0.0, 12.0),
    ("t_lag",                0.0, 12.0),
    ("t_shift",              0.0, 12.0),
    ("t_stationary",         0.0, 12.0),
    ("shape",                0.2,  5.0),
    ("exit_lag_rate",        0.01, 2.0),
    ("inhibition_rate",      0.01, 2.0),
    ("inactivation_rate",    0.01, 2.0),
    ("death_rate",           0.001,0.5),
    ("linear_const",        -0.1,  0.5),
    ("linear_rate",         -0.1,  0.5),
    ("arbitrary_const",     -0.1,  0.5),
    ("doubling_time",        0.5, 24.0),
]

function _ode_param_bounds_for_name(pname::String)::Tuple{Float64,Float64}
    pl = lowercase(pname)
    for (pat, lo, hi) in _ODE_PARAM_BOUNDS
        if occursin(pat, pl)
            return (lo, hi)
        end
    end
    return (0.01, 2.0)
end

function _model_param_bounds(
    model::NLModel,
    tmax::Float64,
    n_points::Int,
    od_scale::Float64 = 1.0,
)::Tuple{Vector{Float64}, Vector{Float64}}
    data_mat = _synthetic_logistic_data(tmax, n_points, od_scale)
    p_hat    = try
        model.guess(data_mat)
    catch
        ones(length(model.param_names))
    end
    lb = max.(0.1 .* abs.(p_hat), 1e-4)
    ub = 10.0 .* abs.(p_hat)
    ub = max.(ub, lb .+ 1e-3)
    return lb, ub
end

function _model_param_bounds(
    model::ODEModel,
    ::Float64,
    ::Int,
    od_scale::Float64 = 1.0,
)::Tuple{Vector{Float64}, Vector{Float64}}
    lb = Vector{Float64}(undef, length(model.param_names))
    ub = Vector{Float64}(undef, length(model.param_names))
    for (j, pname) in enumerate(model.param_names)
        lo, hi = _ode_param_bounds_for_name(pname)
        # Carrying-capacity parameters are in measurement units → scale with od_scale.
        # Growth rates, lag times, shape parameters are dimensionless or time-units — unchanged.
        pl = lowercase(pname)
        if occursin("n_max", pl) || occursin("n_res", pl)
            lo *= od_scale
            hi *= od_scale
        end
        lb[j] = lo
        ub[j] = hi
    end
    return lb, ub
end

# ---------------------------------------------------------------------------
# Simulation helpers
# ---------------------------------------------------------------------------

function _simulate_model(
    model::NLModel,
    p::Vector{Float64},
    t_grid::Vector{Float64},
    ::Float64,
    ::Float64,
    ::Int;
    max_y::Float64 = 1000.0,
)::Union{Tuple{Vector{Float64},Vector{Float64}}, Nothing}
    y = try
        result = model.func(p, t_grid)
        Vector{Float64}(result)
    catch
        try
            [Float64(model.func(p, ti)) for ti in t_grid]
        catch
            return nothing
        end
    end
    return (all(isfinite, y) && maximum(abs, y) < max_y) ? (t_grid, y) : nothing
end

function _simulate_model(
    model::ODEModel,
    p::Vector{Float64},
    t_grid::Vector{Float64},
    od0::Float64,
    tmax::Float64,
    n_points::Int;
    max_y::Float64 = 1000.0,
)::Union{Tuple{Vector{Float64},Vector{Float64}}, Nothing}
    n_eq = model.n_eq
    # Biologically correct initial condition: all biomass starts in the first
    # compartment (lag state for HPM-type models), remaining compartments = 0.
    # This matches the legacy generating_IC convention used during fitting.
    u0 = vcat([od0], zeros(n_eq - 1))
    delta_t = tmax / n_points

    sim = try
        ODE_sim(model.name, u0, 0.0, tmax, delta_t, p)
    catch
        return nothing
    end

    length(sim.t) < 3 && return nothing
    t_sim = Vector{Float64}(sim.t)
    # Observable = total biomass (sum of all compartments).
    # For 1-state models this is just u[1]; for multi-state (HPM family) OD
    # measures lag + active cells so we must sum all states.
    y = [sum(Float64.(u)) for u in sim.u]
    return (all(isfinite, y) && maximum(abs, y) < max_y) ? (t_sim, y) : nothing
end

# ---------------------------------------------------------------------------
# Public: build_model_fingerprint_db
# ---------------------------------------------------------------------------

"""
    build_model_fingerprint_db(; models, n_samples, tmax, n_points, od_range, rng)
        -> ModelFingerprintDB

Build a fingerprint database by Latin-Hypercube-Sampling each model's parameter
space, simulating growth curves, and extracting shape features.

Automatically scales with `od_range`:
- Carrying-capacity parameter bounds (`n_max`, `n_res`) scale proportionally with
  `od_range[2]`, so the LHS covers the right absolute range for the measurement units.
- The numerical-stability discard threshold scales with `od_range[2]`, preventing
  valid fluorescence simulations (e.g. 2 000–20 000 AU) from being incorrectly discarded.

# Keyword arguments
- `models`: models to include (default: all of `MODEL_REGISTRY`).
- `n_samples::Int = 200`: LHS draws per model.
- `tmax::Float64 = 24.0`: simulation end time.
- `n_points::Int = 50`: number of time-grid points.
- `od_range::Tuple{Float64,Float64} = (0.02, 0.05)`: initial-condition range.
  For fluorescence data pass e.g. `od_range = (500.0, 2000.0)`.
- `rng`: random-number generator (default: `Random.GLOBAL_RNG`).
"""
function build_model_fingerprint_db(;
    models   = MODEL_REGISTRY,
    n_samples::Int              = 200,
    tmax::Float64               = 24.0,
    n_points::Int               = 50,
    od_range::Tuple{Float64,Float64} = (0.02, 0.05),
    rng::AbstractRNG            = GLOBAL_RNG,
)::ModelFingerprintDB

    t_grid = collect(range(0.0, tmax; length=n_points))

    # od_scale: multiplier for carrying-capacity bounds (n_max, n_res).
    # Reference is the default od_range[2] = 0.05 (OD scale).
    # For fluorescence od_range[2] = 2000 → od_scale = 40 000.
    od_scale = od_range[2] / 0.05

    # Discard threshold: keeps the default 1000 for OD, scales up for fluorescence.
    # For od_range[2] = 0.05 (OD):        max(1000, 100*0.05)     = 1 000
    # For od_range[2] = 2 000 (fluor.):   max(1000, 100*2000)     = 200 000
    max_y = max(1000.0, 100.0 * od_range[2])

    all_names    = String[]
    all_objects  = AbstractGrowthModel[]
    all_features = Vector{Float64}[]
    all_params   = Vector{Float64}[]
    all_curves   = Vector{Float64}[]   # N_CURVE_PTS Pearson-ready vectors

    model_list = models isa Dict ? values(models) : models

    for model in model_list
        model isa LogLinModel && continue
        model isa DDDEModel   && continue

        lb, ub = try
            _model_param_bounds(model, tmax, n_points, od_scale)
        catch
            continue
        end

        isempty(lb) && continue

        samples = _lhs_sample(lb, ub, n_samples; rng)

        for s in 1:n_samples
            p   = samples[s, :]
            od0 = od_range[1] + rand(rng) * (od_range[2] - od_range[1])

            result = _simulate_model(model, p, t_grid, od0, tmax, n_points; max_y)
            result === nothing && continue
            t_sim, y = result

            # `:features` method — 15-element L2-normalised shape vector
            feat  = _extract_features(t_sim, y)
            fnorm = _norm2(feat)
            fnorm > 1e-12 && (feat ./= fnorm)

            # `:curves` method — Pearson-ready normalised curve (scale-invariant)
            curve_vec = _pearson_ready_curve(t_sim, y, N_CURVE_PTS)

            push!(all_names,    model.name)
            push!(all_objects,  model)
            push!(all_features, feat)
            push!(all_params,   p)
            push!(all_curves,   curve_vec)
        end
    end

    if isempty(all_features)
        return ModelFingerprintDB(
            String[], AbstractGrowthModel[],
            Matrix{Float64}(undef, 0, 15),
            _FEATURE_NAMES,
            Vector{Float64}[],
            Matrix{Float64}(undef, 0, N_CURVE_PTS),
        )
    end

    feat_matrix  = Matrix{Float64}(reduce(hcat, all_features)')   # n_rows × 15
    curve_matrix = Matrix{Float64}(reduce(hcat, all_curves)')     # n_rows × N_CURVE_PTS

    return ModelFingerprintDB(all_names, all_objects, feat_matrix, _FEATURE_NAMES,
                              all_params, curve_matrix)
end

# ---------------------------------------------------------------------------
# Public: recommend_models
# ---------------------------------------------------------------------------

"""
    recommend_models(curve, times, db; top_k, method) -> Vector{String}

Return the names of up to `top_k` models from `db` ranked by similarity to the
query curve.

# Keyword arguments
- `top_k::Int = 3`: number of models to return.
- `method::Symbol = :features`: similarity backend.
  - `:features` — compare 15-element shape-feature vectors via cosine similarity.
    Fast, interpretable, but features are hand-designed and scale-dependent
    (`log_od0`, `log_dynamic_range` encode the measurement scale).
  - `:curves` — compare min-max normalised full curves via Pearson correlation.
    Scale-invariant (works identically for OD and fluorescence), no feature
    engineering, retains full shape information at the cost of 3× more storage.
"""
function recommend_models(
    curve::Vector{Float64},
    times::Vector{Float64},
    db::ModelFingerprintDB;
    top_k::Int     = 3,
    method::Symbol = :features,
)::Vector{String}
    isempty(db.model_names) && return String[]

    sims = if method == :curves
        q = _pearson_ready_curve(times, curve, size(db.curve_matrix, 2))
        db.curve_matrix * q
    else  # :features (default)
        q  = _extract_features(times, curve)
        qn = _norm2(q)
        q_norm = qn > 1e-12 ? q ./ qn : q
        db.features * q_norm
    end

    # Aggregate per unique model: take max similarity across all rows for that model.
    model_sim = Dict{String, Float64}()
    for (i, name) in enumerate(db.model_names)
        model_sim[name] = max(get(model_sim, name, -Inf), sims[i])
    end

    ranked = sort(collect(model_sim); by = x -> -x[2])
    return [name for (name, _) in ranked[1:min(top_k, length(ranked))]]
end

# ---------------------------------------------------------------------------
# Public: suggest_p0
# ---------------------------------------------------------------------------

"""
    suggest_p0(curve, times, db, model_name; method=:features) -> Vector{Float64}

Return a nearest-neighbour initial parameter guess for `model_name` by finding
the DB row for that model with the highest similarity to the query curve, then
returning the stored parameter vector.

# Keyword arguments
- `method::Symbol = :features`: similarity backend (`:features` or `:curves`).
  Use the same method as your `recommend_models` call for consistency.
"""
function suggest_p0(
    curve::Vector{Float64},
    times::Vector{Float64},
    db::ModelFingerprintDB,
    model_name::String;
    method::Symbol = :features,
)::Vector{Float64}
    isempty(db.params) && return Float64[]

    sims = if method == :curves
        q = _pearson_ready_curve(times, curve, size(db.curve_matrix, 2))
        db.curve_matrix * q
    else
        q  = _extract_features(times, curve)
        qn = _norm2(q)
        q_norm = qn > 1e-12 ? q ./ qn : q
        db.features * q_norm
    end

    best_i = 0
    best_s = -Inf
    for (i, name) in enumerate(db.model_names)
        if name == model_name && sims[i] > best_s
            best_s = sims[i]
            best_i = i
        end
    end

    best_i == 0 && return Float64[]
    return copy(db.params[best_i])
end

# ---------------------------------------------------------------------------
# Public: score_models  (diagnostic)
# ---------------------------------------------------------------------------

"""
    score_models(curve, times, db; method=:features) -> Vector{Pair{String,Float64}}

Return every model in the DB ranked by its maximum similarity score to the
query curve. Useful for diagnosing why a specific model is or is not being
recommended — compare its rank and score against the top models.

# Example
```julia
scores = score_models(y, t, db; method=:features)
for (name, score) in scores[1:10]
    println(rpad(name, 32), round(score; digits=4))
end
```
"""
function score_models(
    curve::Vector{Float64},
    times::Vector{Float64},
    db::ModelFingerprintDB;
    method::Symbol = :features,
)::Vector{Pair{String,Float64}}
    isempty(db.model_names) && return Pair{String,Float64}[]

    sims = if method == :curves
        q = _pearson_ready_curve(times, curve, size(db.curve_matrix, 2))
        db.curve_matrix * q
    else
        q  = _extract_features(times, curve)
        qn = _norm2(q)
        q_norm = qn > 1e-12 ? q ./ qn : q
        db.features * q_norm
    end

    model_sim = Dict{String, Float64}()
    for (i, name) in enumerate(db.model_names)
        model_sim[name] = max(get(model_sim, name, -Inf), sims[i])
    end

    sorted = sort(collect(model_sim); by = x -> -x[2])
    return [k => v for (k, v) in sorted]
end

# ---------------------------------------------------------------------------
# Internal: _default_p0  (fallback when no DB is available)
# ---------------------------------------------------------------------------

function _default_p0(model::AbstractGrowthModel, data::GrowthData)::Vector{Float64}
    if model isa NLModel
        model.guess !== nothing || return Float64[]
        data_mat = Matrix(transpose(hcat(data.times, data.curves[1, :])))
        return try model.guess(data_mat) catch; Float64[] end
    elseif model isa ODEModel
        tmax_d = isempty(data.times) ? 24.0 : Float64(maximum(data.times))
        lb, ub = try
            _model_param_bounds(model, tmax_d, length(data.times))
        catch
            return ones(length(model.param_names))
        end
        return (lb .+ ub) ./ 2.0
    end
    return Float64[]
end

# ---------------------------------------------------------------------------
# Public: smart_fit
# ---------------------------------------------------------------------------

"""
    smart_fit(data, db[, opts]; top_k, method) -> GrowthFitResults

Recommend models for each curve in `data` using `db`, take the union of
recommendations, and fit them via [`kinbiont_fit`](@ref).

Initial parameters are seeded from the nearest-neighbour fingerprint row
via [`suggest_p0`](@ref).

# Keyword arguments
- `top_k::Int = 3`: number of models to recommend per curve.
- `method::Symbol = :features`: passed to [`recommend_models`](@ref).
  Use `:curves` for scale-invariant recommendations (fluorescence data).
"""
function smart_fit(
    data::GrowthData,
    db::ModelFingerprintDB,
    opts::FitOptions = FitOptions();
    top_k::Int     = 3,
    method::Symbol = :features,
)::GrowthFitResults
    all_names   = String[]
    query_curve = data.curves[1, :]
    for i in axes(data.curves, 1)
        recs = recommend_models(data.curves[i, :], data.times, db; top_k, method)
        append!(all_names, recs)
    end
    unique!(all_names)

    isempty(all_names) && error("smart_fit: no models recommended (DB may be empty)")

    models = [MODEL_REGISTRY[n] for n in all_names]
    params = map(all_names) do name
        p = suggest_p0(query_curve, data.times, db, name)
        isempty(p) ? _default_p0(MODEL_REGISTRY[name], data) : p
    end
    spec = ModelSpec(models, params)
    return kinbiont_fit(data, spec, opts)
end

# ---------------------------------------------------------------------------
# Exports
# ---------------------------------------------------------------------------

export ModelFingerprintDB
export build_model_fingerprint_db
export recommend_models
export suggest_p0
export score_models
export smart_fit
