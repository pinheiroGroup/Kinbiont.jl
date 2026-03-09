# =============================================================================
# Model Recommendation via Feature Fingerprinting
# =============================================================================
# Public API:
#   build_model_fingerprint_db(; ...) -> ModelFingerprintDB
#   recommend_models(curve, times, db; top_k=3) -> Vector{String}
#   smart_fit(data, db[, opts]; top_k=3) -> GrowthFitResults
# =============================================================================

using Random: AbstractRNG, GLOBAL_RNG, randperm
using Statistics: var, mean

_norm2(v) = sqrt(sum(x^2 for x in v))

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
- `features::Matrix{Float64}`: `n_rows × 15` L2-normalised feature matrix.
- `feature_names::Vector{String}`: human-readable name for each of the 15 features.
- `params::Vector{Vector{Float64}}`: parameter vector used to generate each row's simulation.
  Used by [`suggest_p0`](@ref) to seed the optimiser with a nearest-neighbour initial guess.
"""
struct ModelFingerprintDB
    model_names  ::Vector{String}
    model_objects::Vector{AbstractGrowthModel}
    features     ::Matrix{Float64}   # n_rows × N_FEATURES, L2-normalised
    feature_names::Vector{String}
    params       ::Vector{Vector{Float64}}   # parameter vector per row
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
# Feature extraction
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
    # Need at least 2 points to compute any meaningful feature
    if n < 2
        return zeros(15)
    end

    tmax  = t[end] - t[1]
    ymax  = maximum(y)
    ymin  = minimum(y)

    # --- specific growth rates (on log-scale) ---
    logy = log.(max.(y, 1e-12))
    spec_gr = Vector{Float64}(undef, n - 1)
    for i in 1:(n-1)
        dt = t[i+1] - t[i]
        spec_gr[i] = dt > 0 ? (logy[i+1] - logy[i]) / dt : 0.0
    end

    if isempty(spec_gr)
        return zeros(15)
    end

    max_sgr     = maximum(spec_gr)
    igr         = argmax(spec_gr)   # index in spec_gr (length n-1)
    igr_t       = t[igr]            # time of max spec gr

    # 1. log_dynamic_range
    f1 = log10(max(ymax / max(ymin, 1e-9), 1.0))

    # 2. log_od0
    f2 = log10(max(y[1], 1e-9))

    # 3. max_spec_gr
    f3 = max_sgr

    # 4. t_max_gr_norm
    f4 = tmax > 0 ? (igr_t - t[1]) / tmax : 0.0

    # 5. lag_est_norm — time when y >= y[1] + 0.1*(ymax-y[1])
    thresh_lag = y[1] + 0.1 * (ymax - y[1])
    idx_lag    = findfirst(yi -> yi >= thresh_lag, y)
    f5 = (idx_lag !== nothing && tmax > 0) ? (t[idx_lag] - t[1]) / tmax : 0.0

    # 6. stationary_est_norm — time when y >= y[1] + 0.9*(ymax-y[1])
    thresh_stat = y[1] + 0.9 * (ymax - y[1])
    idx_stat    = findfirst(yi -> yi >= thresh_stat, y)
    f6 = (idx_stat !== nothing && tmax > 0) ? (t[idx_stat] - t[1]) / tmax : 0.0

    # 7. auc_norm
    total_auc = _trapz(t, y)
    f7 = (tmax > 0 && ymax > 0) ? total_auc / (tmax * ymax) : 0.0

    # 8. gr_variance
    f8 = length(spec_gr) > 1 ? var(spec_gr) : 0.0

    # 9. gr_skewness
    f9 = if length(spec_gr) >= 3
        mu  = mean(spec_gr)
        sig = sqrt(max(var(spec_gr), 1e-30))
        mean((s - mu)^3 for s in spec_gr) / sig^3
    else
        0.0
    end

    # 10. late_slope — OLS slope on last 20% of time points
    n_late = max(2, round(Int, 0.2 * n))
    t_late = t[(end - n_late + 1):end]
    y_late = y[(end - n_late + 1):end]
    tm     = mean(t_late)
    ym     = mean(y_late)
    denom  = sum((ti - tm)^2 for ti in t_late)
    f10    = denom > 0 ? sum((t_late[i] - tm) * (y_late[i] - ym) for i in eachindex(t_late)) / denom : 0.0

    # 11. od_end_ratio
    f11 = ymax > 0 ? y[end] / ymax : 0.0

    # 12. curvature_at_inflection — second finite diff at max-gr index
    # igr is index in spec_gr; map back to y index = igr (since spec_gr[i] = logy[i+1]-logy[i])
    yi = min(igr + 1, n)  # corresponding y index
    f12 = if 2 <= yi <= n-1
        (y[yi+1] - 2*y[yi] + y[yi-1]) / max((t[yi] - t[yi-1])^2, 1e-12)
    else
        0.0
    end

    # 13. early_auc_frac
    mid_t   = t[1] + 0.5 * tmax
    idx_mid = searchsortedlast(t, mid_t)
    idx_mid = clamp(idx_mid, 2, n)
    early_auc = _trapz(t[1:idx_mid], y[1:idx_mid])
    f13 = total_auc > 1e-12 ? early_auc / total_auc : 0.0

    # 14. n_sign_changes_gr
    f14 = Float64(sum(spec_gr[i] * spec_gr[i+1] < 0 for i in 1:(length(spec_gr)-1)))

    # 15. exp_phase_fraction
    f15 = max_sgr > 0 ? count(s -> s > 0.5 * max_sgr, spec_gr) / length(spec_gr) : 0.0

    feat = [f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, f14, f15]
    return map(x -> isfinite(x) ? x : 0.0, feat)
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
        perm      = randperm(rng, n)
        X[:, j]   = lb[j] .+ (ub[j] - lb[j]) .* ((perm .- rand(rng, n)) ./ n)
    end
    return X   # n × d
end

# ---------------------------------------------------------------------------
# Parameter bounds
# ---------------------------------------------------------------------------

# Synthetic logistic curve for NL model guess
function _synthetic_logistic_data(tmax::Float64, n_points::Int)::Matrix{Float64}
    t  = range(0.0, tmax; length=n_points)
    N0 = 0.03; K = 1.2; r = 0.3
    y  = @. K / (1 + (K/N0 - 1) * exp(-r * t))
    return Matrix(transpose(hcat(collect(t), y)))  # 2 × n_points
end

const _ODE_PARAM_BOUNDS = [
    # (substring, lb, ub) — case-insensitive
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
    return (0.01, 2.0)   # default
end

function _model_param_bounds(
    model::NLModel,
    tmax::Float64,
    n_points::Int,
)::Tuple{Vector{Float64}, Vector{Float64}}
    data_mat = _synthetic_logistic_data(tmax, n_points)
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
)::Tuple{Vector{Float64}, Vector{Float64}}
    lb = Vector{Float64}(undef, length(model.param_names))
    ub = Vector{Float64}(undef, length(model.param_names))
    for (j, pname) in enumerate(model.param_names)
        lo, hi  = _ode_param_bounds_for_name(pname)
        lb[j]   = lo
        ub[j]   = hi
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
    ::Float64,          # od0 (unused for NL)
    ::Float64,          # tmax (unused)
    ::Int,              # n_points (unused)
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
    return (all(isfinite, y) && maximum(abs, y) < 1000) ? (t_grid, y) : nothing
end

function _simulate_model(
    model::ODEModel,
    p::Vector{Float64},
    t_grid::Vector{Float64},
    od0::Float64,
    tmax::Float64,
    n_points::Int;
    ic_alpha::Float64 = 0.9,   # fraction of od0 in first compartment (sampled per draw)
)::Union{Tuple{Vector{Float64},Vector{Float64}}, Nothing}
    n_eq  = model.n_eq
    u0    = if n_eq == 1
        [od0]
    elseif n_eq == 2
        # ic_alpha is sampled in [0.05, 0.95] so the DB covers both lag-heavy
        # (small ic_alpha) and active-heavy (large ic_alpha) initial conditions,
        # giving 2-state models like aHPM a distinctive fingerprint region.
        [ic_alpha * od0, (1.0 - ic_alpha) * od0]
    else
        alpha_rest = (1.0 - ic_alpha) / (n_eq - 1)
        vcat([ic_alpha * od0], fill(alpha_rest * od0, n_eq - 1))
    end
    delta_t = tmax / n_points

    sim = try
        ODE_sim(model.name, u0, 0.0, tmax, delta_t, p)
    catch
        return nothing
    end

    # Extract first state variable; use actual sim.t to avoid length mismatch
    length(sim.t) < 3 && return nothing   # aborted simulation has too few points
    t_sim = Vector{Float64}(sim.t)
    y     = [Float64(u[1]) for u in sim.u]
    return (all(isfinite, y) && maximum(abs, y) < 1000) ? (t_sim, y) : nothing
end

# ---------------------------------------------------------------------------
# Public: build_model_fingerprint_db
# ---------------------------------------------------------------------------

"""
    build_model_fingerprint_db(; models, n_samples, tmax, n_points, od_range, rng)
        -> ModelFingerprintDB

Build a fingerprint database by Latin-Hypercube-Sampling each model's parameter
space, simulating growth curves, and extracting shape features.

# Keyword arguments
- `models`: collection of models to include (default: all of `MODEL_REGISTRY`).
- `n_samples::Int = 200`: LHS draws per model.
- `tmax::Float64 = 24.0`: simulation end time.
- `n_points::Int = 50`: number of time-grid points.
- `od_range::Tuple{Float64,Float64} = (0.02, 0.05)`: range for initial OD sampling.
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

    all_names    = String[]
    all_objects  = AbstractGrowthModel[]
    all_features = Vector{Float64}[]
    all_params   = Vector{Float64}[]

    model_list = models isa Dict ? values(models) : models

    for model in model_list
        # Skip LogLinModel and DDDEModel — not simulatable with a parameter vector
        model isa LogLinModel  && continue
        model isa DDDEModel    && continue

        lb, ub = try
            _model_param_bounds(model, tmax, n_points)
        catch
            continue
        end

        # Guard: ensure lb/ub have at least 1 element
        isempty(lb) && continue

        samples = _lhs_sample(lb, ub, n_samples; rng)

        is_multi_state = model isa ODEModel && model.n_eq >= 2

        for s in 1:n_samples
            p   = samples[s, :]
            od0 = od_range[1] + rand(rng) * (od_range[2] - od_range[1])

            # Sample ic_alpha in [0.05, 0.95] for multi-state models so the DB
            # covers the full range of lag-heavy vs active-heavy initial conditions,
            # giving 2-state models like aHPM a distinctive fingerprint region.
            result = if is_multi_state
                ic_alpha = 0.05 + rand(rng) * 0.90
                _simulate_model(model, p, t_grid, od0, tmax, n_points; ic_alpha)
            else
                _simulate_model(model, p, t_grid, od0, tmax, n_points)
            end
            result === nothing && continue
            t_sim, y = result

            feat = _extract_features(t_sim, y)
            fnorm = _norm2(feat)
            fnorm > 1e-12 && (feat ./= fnorm)

            push!(all_names,    model.name)
            push!(all_objects,  model)
            push!(all_features, feat)
            push!(all_params,   p)
        end
    end

    if isempty(all_features)
        return ModelFingerprintDB(String[], AbstractGrowthModel[],
                                  Matrix{Float64}(undef, 0, 15), _FEATURE_NAMES,
                                  Vector{Float64}[])
    end

    feat_matrix = Matrix{Float64}(reduce(hcat, all_features)')  # n_rows × 15

    return ModelFingerprintDB(all_names, all_objects, feat_matrix, _FEATURE_NAMES, all_params)
end

# ---------------------------------------------------------------------------
# Public: recommend_models
# ---------------------------------------------------------------------------

"""
    recommend_models(curve, times, db; top_k=3) -> Vector{String}

Given an observed growth curve, return the names of up to `top_k` models
from `db` ranked by cosine similarity to the query's feature vector.
"""
function recommend_models(
    curve::Vector{Float64},
    times::Vector{Float64},
    db::ModelFingerprintDB;
    top_k::Int = 3,
)::Vector{String}
    isempty(db.model_names) && return String[]

    q      = _extract_features(times, curve)
    qn     = _norm2(q)
    q_norm = qn > 1e-12 ? q ./ qn : q

    # Cosine similarities (db rows already L2-normalised)
    sims = db.features * q_norm   # length n_rows

    # Aggregate per unique model: take max similarity
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
    suggest_p0(curve, times, db, model_name) -> Vector{Float64}

Return a nearest-neighbour initial parameter guess for `model_name` by finding
the DB row for that model whose feature vector is most similar (cosine similarity)
to the query curve's feature vector, then returning the parameter vector that
was used to simulate that row.

Falls back to the midpoint of parameter bounds if the DB has no rows for the
requested model or if the DB was built without parameter storage.
"""
function suggest_p0(
    curve::Vector{Float64},
    times::Vector{Float64},
    db::ModelFingerprintDB,
    model_name::String,
)::Vector{Float64}
    isempty(db.params) && return Float64[]

    q      = _extract_features(times, curve)
    qn     = _norm2(q)
    q_norm = qn > 1e-12 ? q ./ qn : q

    sims = db.features * q_norm   # cosine similarity for every row

    # Find the best-matching row that belongs to this model
    best_i = 0
    best_s = -Inf
    for (i, name) in enumerate(db.model_names)
        if name == model_name && sims[i] > best_s
            best_s = sims[i]
            best_i = i
        end
    end

    best_i == 0 && return Float64[]   # model not in DB
    return copy(db.params[best_i])
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
    smart_fit(data, db[, opts]; top_k=3) -> GrowthFitResults

Recommend models for each curve in `data` using `db`, take the union of
recommendations, and fit them via [`kinbiont_fit`](@ref).

Initial parameters are seeded from the nearest-neighbour fingerprint row
via [`suggest_p0`](@ref), giving the optimiser a warm start close to a
parameter set that already produced a similar curve shape.
"""
function smart_fit(
    data::GrowthData,
    db::ModelFingerprintDB,
    opts::FitOptions = FitOptions();
    top_k::Int = 3,
)::GrowthFitResults
    all_names = String[]
    # Use the first curve as the representative query for p0 suggestion
    query_curve = data.curves[1, :]
    for i in axes(data.curves, 1)
        recs = recommend_models(data.curves[i, :], data.times, db; top_k)
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
export smart_fit
