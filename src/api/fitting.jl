# =============================================================================
# Unified fit() Entry Point
# =============================================================================
# Single public function:
#
#   fit(data::GrowthData, spec::ModelSpec[, opts::FitOptions]) -> GrowthFitResults
#
# Dispatches internally to the existing low-level fitting functions:
#   NLModel  -> fit_NL_model (NL_fit_one_well.jl)
#   ODEModel -> fitting_one_well_custom_ODE (Fit_one_well_functions.jl)
#   LogLinModel -> fitting_one_well_Log_Lin (Fit_one_well_functions.jl)
#
# AICc-based model selection is performed for each curve across all candidate
# models in the ModelSpec.
# =============================================================================

# Internal result type for a single (curve, model) pair — not exported
const _CandidateResult = @NamedTuple begin
    model_name::String
    params::Vector{Any}
    fitted_curve::Vector{Float64}
    times::Vector{Float64}
    loss::Float64
    aic::Float64
end

# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

"""
    kinbiont_fit(data::GrowthData, spec::ModelSpec[, opts::FitOptions]) -> GrowthFitResults

Fit every curve in `data` to every model in `spec` and select the best fit per
curve using AICc. Returns a [`GrowthFitResults`](@ref) containing one
[`CurveFitResult`](@ref) per curve.

Preprocessing (smoothing, blank subtraction, clustering, …) is applied
according to `opts` before fitting. Pass `FitOptions()` for all defaults.

# Example
```julia
data  = GrowthData(my_matrix, times, labels)
spec  = ModelSpec([MODEL_REGISTRY["NL_logistic"]], [[1.2, 0.5, 0.01]])
opts  = FitOptions(smooth=true, loss="RE")

results = kinbiont_fit(data, spec, opts)

for r in results
    println(r.label, " → ", r.best_model.name, "  AICc=", round(r.best_aic; digits=2))
end
```
"""
function kinbiont_fit(
    data::GrowthData,
    spec::ModelSpec,
    opts::FitOptions = FitOptions(),
)::GrowthFitResults

    _validate_spec(data, spec)
    processed = preprocess(data, opts)

    results = map(axes(processed.curves, 1)) do i
        curve = processed.curves[i, :]
        label = processed.labels[i]
        _fit_curve(curve, processed.times, label, spec, opts)
    end

    return GrowthFitResults(processed, results)
end

# ---------------------------------------------------------------------------
# Per-curve fitting and model selection
# ---------------------------------------------------------------------------

function _fit_curve(
    curve::Vector{Float64},
    times::Vector{Float64},
    label::String,
    spec::ModelSpec,
    opts::FitOptions,
)::CurveFitResult

    # Legacy functions expect a 2×n matrix [times; values]
    data_mat = Matrix(transpose(hcat(times, curve)))

    candidates = map(eachindex(spec.models)) do k
        model  = spec.models[k]
        p0     = _resolve_p0(spec.params[k], model, data_mat)
        lb     = spec.lower === nothing ? nothing : spec.lower[k]
        ub     = spec.upper === nothing ? nothing : spec.upper[k]
        _fit_single(data_mat, label, model, p0, lb, ub, opts)
    end

    best_idx = argmin([r.aic for r in candidates])   # Int index
    best = candidates[best_idx]

    return CurveFitResult(
        label,
        spec.models[best_idx],
        best.params,
        best.aic,
        best.fitted_curve,
        best.times,
        best.loss,
        candidates,
    )
end

# LogLinModel takes no parameters — always return empty vector
_resolve_p0(::Vector{Float64}, ::LogLinModel, ::Matrix{Float64}) = Float64[]

# Use model's built-in guess if p0 is empty and a guess function exists
function _resolve_p0(
    p0::Vector{Float64},
    model::AbstractGrowthModel,
    data_mat::Matrix{Float64},
)::Vector{Float64}
    isempty(p0) || return p0
    m = model
    if hasproperty(m, :guess) && m.guess !== nothing
        return m.guess(data_mat)
    end
    error("$(m.name): empty params passed and no guess function available")
end

# ---------------------------------------------------------------------------
# Dispatch on model type — no if/elseif strings
# ---------------------------------------------------------------------------

function _fit_single(
    data_mat::Matrix{Float64},
    label::String,
    model::NLModel,
    p0::Vector{Float64},
    lb,
    ub,
    opts::FitOptions,
)::_CandidateResult

    opt_params = _build_opt_params(lb, ub, opts.opt_params)

    raw = fit_NL_model(
        data_mat,
        model.name,
        label,
        model.func,
        p0;
        smoothing              = false,   # preprocessing already done
        type_of_loss           = opts.loss,
        pt_smooth_derivative   = opts.pt_smooth_derivative,
        optimizer              = opts.optimizer,
        multistart             = opts.multistart,
        n_restart              = opts.n_restart,
        auto_diff_method       = opts.auto_diff_method,
        cons                   = opts.cons,
        opt_params...,
    )

    # raw = ("NL", res_param, fitted_model, times)
    # res_param layout: [label, well, model_name, p1..pn, th_gr, em_gr, loss]
    res_param    = raw[2]
    fitted_curve = Vector{Float64}(raw[3])
    fit_times    = Vector{Float64}(raw[4])
    loss_val     = Float64(res_param[end])
    params_raw   = res_param[4:end-3]     # strip metadata columns
    n_params     = length(p0)
    aic          = AICc_evaluation2(n_params, 2.0, fit_times, loss_val;
                                    correction = opts.aic_correction)

    return (
        model_name   = model.name,
        params       = Vector{Any}(params_raw),
        fitted_curve = fitted_curve,
        times        = fit_times,
        loss         = loss_val,
        aic          = Float64(aic),
    )
end

function _fit_single(
    data_mat::Matrix{Float64},
    label::String,
    model::ODEModel,
    p0::Vector{Float64},
    lb,
    ub,
    opts::FitOptions,
)::_CandidateResult

    opt_params = _build_opt_params(lb, ub, opts.opt_params)

    raw = fitting_one_well_custom_ODE(
        data_mat,
        model.name,
        label,
        model.func,
        p0,
        model.n_eq;
        Integration_method     = opts.integrator,
        smoothing              = false,   # preprocessing already done
        type_of_loss           = opts.loss,
        pt_smooth_derivative   = opts.pt_smooth_derivative,
        optimizer              = opts.optimizer,
        multistart             = opts.multistart,
        n_restart              = opts.n_restart,
        auto_diff_method       = opts.auto_diff_method,
        cons                   = opts.cons,
        opt_params...,
    )

    # raw = ("custom_ODE", res_param, fitted_curve_values, fitted_times)
    # res_param: [name_well, "custom_model", params_vector, th_gr, em_gr, loss]
    res_param    = raw[2]
    fitted_curve = Vector{Float64}(raw[3])
    fit_times    = Vector{Float64}(raw[4])
    loss_val     = Float64(res_param[end])
    params_raw   = res_param[3]           # Vector of fitted parameters
    n_params     = length(p0)
    aic          = AICc_evaluation2(n_params, 2.0, fit_times, loss_val;
                                    correction = opts.aic_correction)

    return (
        model_name   = model.name,
        params       = Vector{Any}(params_raw),
        fitted_curve = fitted_curve,
        times        = fit_times,
        loss         = loss_val,
        aic          = Float64(aic),
    )
end

function _fit_single(
    data_mat::Matrix{Float64},
    label::String,
    ::LogLinModel,
    _p0::Vector{Float64},
    _lb,
    _ub,
    opts::FitOptions,
)::_CandidateResult

    raw = fitting_one_well_Log_Lin(
        data_mat,
        "log_lin",
        label;
        type_of_smoothing      = "NO",    # preprocessing already done
        pt_smoothing_derivative = opts.pt_smooth_derivative,
    )

    # raw is a Tuple: (method, params, fit_matrix, smoothed_data, confidence_band)
    # fit_matrix is n×2: columns are [times, log_fitted_values]
    params_raw   = raw[2]
    fit_matrix   = raw[3]
    fit_times    = Vector{Float64}(fit_matrix[:, 1])
    fitted_curve = Vector{Float64}(fit_matrix[:, 2])
    # Log-lin loss is the residual of the linear regression
    loss_val     = 1.0 - Float64(params_raw[end])   # 1 - R²
    n_params     = 2   # slope + intercept
    aic          = AICc_evaluation2(n_params, 2.0, fit_times, max(loss_val, 1e-12);
                                    correction = opts.aic_correction)

    return (
        model_name   = "log_lin",
        params       = Vector{Any}(params_raw),
        fitted_curve = fitted_curve,
        times        = fit_times,
        loss         = loss_val,
        aic          = Float64(aic),
    )
end

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

function _build_opt_params(lb, ub, extra::NamedTuple = (;))
    bounds = if lb === nothing && ub === nothing
        (;)
    elseif lb === nothing
        (ub = ub,)
    elseif ub === nothing
        (lb = lb,)
    else
        (lb = lb, ub = ub)
    end
    return merge(bounds, extra)
end

function _validate_spec(data::GrowthData, spec::ModelSpec)
    isempty(spec.models) && error("ModelSpec: models list is empty")
    n = length(spec.models)
    length(spec.params) == n || error(
        "ModelSpec: params length $(length(spec.params)) ≠ models length $n"
    )
end

export kinbiont_fit
