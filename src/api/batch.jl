# =============================================================================
# GUI-compatible batch fitting helpers
# =============================================================================

import OptimizationNLopt

const BATCH_OPTIMIZER_MAP = Dict{String, Any}(
    "LN_BOBYQA" => OptimizationNLopt.NLopt.LN_BOBYQA,
    "LN_COBYLA" => OptimizationNLopt.NLopt.LN_COBYLA,
    "GN_ISRES" => OptimizationNLopt.NLopt.GN_ISRES,
    "GN_DIRECT_L" => OptimizationNLopt.NLopt.GN_DIRECT_L,
    "BBO_adaptive_de_rand_1_bin_radiuslimited" => BBO_adaptive_de_rand_1_bin_radiuslimited(),
)

_batch_optimizer(name::String) = get(BATCH_OPTIMIZER_MAP, name, BATCH_OPTIMIZER_MAP["LN_BOBYQA"])

function _batch_model_name(m)
    m isa LogLinModel && return "log_lin"
    return m.name
end

function _batch_csv_cell(x)
    if x === missing || x === nothing
        return "\"\""
    elseif x isa AbstractFloat && !isfinite(x)
        return "\"\""
    end
    return "\"" * replace(string(x), "\"" => "\"\"") * "\""
end

function _batch_write_csv(path::String, header, rows)
    open(path, "w") do io
        println(io, join(_batch_csv_cell.(header), ","))
        for row in rows
            println(io, join(_batch_csv_cell.(row), ","))
        end
    end
end

function _batch_format_time(t)
    !isfinite(t) && return string(t)
    isapprox(t, round(t); atol=1e-10) && return string(Int(round(t)))
    return string(Float64(round(t; sigdigits=12)))
end

function _batch_quantile(values, p)
    xs = sort(Float64[v for v in values if isfinite(v)])
    isempty(xs) && return NaN
    length(xs) == 1 && return xs[1]
    pos = 1 + (length(xs) - 1) * clamp(Float64(p), 0.0, 1.0)
    lo = floor(Int, pos)
    hi = ceil(Int, pos)
    lo == hi && return xs[lo]
    return xs[lo] + (pos - lo) * (xs[hi] - xs[lo])
end

_batch_positive_or(x, fallback) = isfinite(x) && x > 0 ? Float64(x) : Float64(fallback)
_batch_clip_initial(x; lower=1e-6, upper=49.0) =
    isfinite(x) ? clamp(Float64(x), Float64(lower), Float64(upper)) : 1.0

function _batch_growth_features(time_numeric::Vector{Float64}, od_for_fit::Vector{Float64})
    valid = findall(i -> isfinite(time_numeric[i]) && isfinite(od_for_fit[i]), eachindex(time_numeric))
    if length(valid) < 2
        return Dict{Symbol, Float64}(
            :baseline => 0.01, :plateau => 1.0, :amplitude => 1.0,
            :growth_rate => 1.0, :max_slope => 1.0, :lag_time => 0.0,
            :mid_time => 0.0, :inflection_time => 0.0,
            :doubling_time => log(2), :duration => 1.0,
            :terminal_slope => 0.0, :q0 => 1.0,
        )
    end

    ord = sortperm(time_numeric[valid])
    t = time_numeric[valid][ord]
    y = max.(od_for_fit[valid][ord], 1e-6)
    duration = _batch_positive_or(t[end] - t[1], 1.0)
    n_edge = clamp(ceil(Int, length(y) * 0.1), 2, min(8, length(y)))
    baseline = max(_batch_quantile(y[1:n_edge], 0.5), minimum(y), 1e-4)
    late_level = _batch_quantile(y[end-n_edge+1:end], 0.5)
    plateau = max(_batch_quantile(y, 0.95), late_level, baseline + 1e-4)
    amplitude = max(plateau - baseline, maximum(y) - minimum(y), 1e-4)

    slopes = Float64[]
    slope_times = Float64[]
    log_slopes = Float64[]
    for i in 1:(length(t)-1)
        dt = t[i+1] - t[i]
        dt > 0 || continue
        push!(slopes, (y[i+1] - y[i]) / dt)
        push!(slope_times, (t[i] + t[i+1]) / 2)
        push!(log_slopes, (log(y[i+1]) - log(y[i])) / dt)
    end

    pos_slopes = filter(x -> isfinite(x) && x > 0, slopes)
    max_slope = isempty(pos_slopes) ? amplitude / duration : maximum(pos_slopes)
    pos_log_slopes = filter(x -> isfinite(x) && x > 0, log_slopes)
    max_log_slope = isempty(pos_log_slopes) ? 0.0 : maximum(pos_log_slopes)
    growth_rate = _batch_clip_initial(max(max_log_slope, 4 * max_slope / max(plateau, 1e-4), 1 / duration); upper=20.0)

    function first_cross(thr)
        idx = findfirst(v -> v >= thr, y)
        idx === nothing ? t[1] : t[idx]
    end
    lag_time = clamp(first_cross(baseline + 0.10 * amplitude), t[1], t[end])
    mid_time = clamp(first_cross(baseline + 0.50 * amplitude), t[1], t[end])
    inflection = isempty(slope_times) ? t[1] : slope_times[argmax(slopes)]
    tail_start = max(1, length(slopes) - max(1, ceil(Int, 0.2 * length(slopes))) + 1)
    terminal_slope = isempty(slopes) ? 0.0 : _batch_quantile(slopes[tail_start:end], 0.5)

    return Dict{Symbol, Float64}(
        :baseline => _batch_clip_initial(baseline),
        :plateau => _batch_clip_initial(plateau),
        :amplitude => _batch_clip_initial(amplitude),
        :growth_rate => growth_rate,
        :max_slope => _batch_clip_initial(max_slope),
        :lag_time => _batch_clip_initial(lag_time - t[1]; upper=duration),
        :mid_time => _batch_clip_initial(mid_time - t[1]; upper=duration),
        :inflection_time => _batch_clip_initial(inflection - t[1]; upper=duration),
        :doubling_time => _batch_clip_initial(log(2) / growth_rate; upper=duration),
        :duration => duration,
        :terminal_slope => _batch_clip_initial(abs(terminal_slope); upper=20.0),
        :q0 => _batch_clip_initial(exp(-growth_rate * max(lag_time - t[1], 0.0)); upper=50.0),
    )
end

const BATCH_MODEL_INIT = Dict{String, Function}(
    "aHPM" => f -> [f[:growth_rate], 1.0/max(f[:lag_time], 0.1), f[:plateau], 1.0],
    "HPM" => f -> [f[:growth_rate], 1.0/max(f[:lag_time], 0.1), f[:plateau]],
    "HPM_exp" => f -> [f[:growth_rate], 1.0/max(f[:lag_time], 0.1)],
    "HPM_3_death" => f -> [f[:growth_rate], 1.0/max(f[:lag_time], 0.1), 0.05*f[:growth_rate], 0.05*f[:growth_rate]],
    "HPM_3_inhibition" => f -> [f[:growth_rate], 1.0/max(f[:lag_time], 0.1), 0.05*f[:growth_rate]],
    "HPM_inhibition" => f -> [f[:growth_rate], 0.05*f[:growth_rate], 0.5*f[:growth_rate], f[:plateau]],
    "aHPM_3_death_resistance" => f -> [f[:growth_rate], 1.0/max(f[:lag_time], 0.1), 0.05*f[:growth_rate], 0.05*f[:growth_rate], 1.0, 1.0],
    "aHPM_inhibition" => f -> [f[:growth_rate], 0.05*f[:growth_rate], 0.5*f[:growth_rate], f[:plateau], 1.0],
    "logistic" => f -> [f[:growth_rate], f[:plateau]],
    "gompertz" => f -> [f[:growth_rate], f[:plateau]],
    "exponential" => f -> [f[:growth_rate]],
    "alogistic" => f -> [f[:growth_rate], f[:plateau], 1.0],
    "hyper_gompertz" => f -> [f[:growth_rate], f[:plateau], 1.0],
    "hyper_logistic" => f -> [f[:doubling_time], f[:growth_rate], f[:plateau], 1.0],
    "bertalanffy_richards" => f -> [f[:growth_rate], f[:plateau], 1.0],
    "baranyi_exp" => f -> [f[:growth_rate], f[:lag_time], 1.0],
    "baranyi_richards" => f -> [f[:growth_rate], f[:plateau], f[:lag_time], 1.0],
    "baranyi_roberts" => f -> [f[:growth_rate], f[:plateau], f[:lag_time], 1.0, 1.0],
    "NL_exponential" => f -> [f[:baseline], f[:growth_rate]],
    "NL_logistic" => f -> [f[:plateau], f[:growth_rate], f[:lag_time]],
    "NL_Gompertz" => f -> [f[:plateau], f[:growth_rate], f[:lag_time]],
    "NL_Richards" => f -> [f[:plateau], 1.0, f[:growth_rate], f[:lag_time]],
    "NL_Weibull" => f -> [f[:plateau], f[:baseline], f[:growth_rate], 1.0],
    "NL_Bertalanffy" => f -> [f[:baseline], f[:plateau], f[:growth_rate], 1.0],
)

function _batch_initial_value(param_name, features)
    # Priority-ordered substring match: more specific suffix/prefix tokens
    # (e.g. "nlag", "death", "doubl") are checked before broader ones
    # ("lag", "growth"). Reordering this block changes which initial value
    # ambiguous names like `growth_rate_max` receive.
    compact = replace(lowercase(String(param_name)), r"[^a-z0-9]" => "")
    if compact in ("nlag", "xlag", "ylag", "odlag") ||
       ((startswith(compact, "n") || startswith(compact, "x") ||
         startswith(compact, "y") || startswith(compact, "od")) && occursin("lag", compact))
        return features[:baseline]
    elseif occursin("death", compact) || occursin("decay", compact) ||
           occursin("decline", compact) || occursin("mort", compact) ||
           occursin("inhib", compact) || occursin("inactiv", compact) ||
           occursin("resist", compact)
        return _batch_clip_initial(0.05 * features[:growth_rate]; upper=20.0)
    elseif occursin("doubl", compact) || compact in ("dt", "td")
        return features[:doubling_time]
    elseif occursin("slope", compact) || occursin("linear", compact)
        return max(features[:terminal_slope], 1e-6)
    elseif compact in ("n0", "x0", "y0", "od0") || occursin("initial", compact) ||
           occursin("inoc", compact) || occursin("baseline", compact)
        return features[:baseline]
    elseif compact in ("lag", "tl", "tlag", "lagtime") || startswith(compact, "tlag") ||
           occursin("lambda", compact) || occursin("delay", compact) || occursin("lag", compact)
        return features[:lag_time]
    elseif occursin("mu", compact) || compact in ("r", "gr") ||
           occursin("growth", compact) || occursin("rate", compact) || occursin("gr", compact)
        return features[:growth_rate]
    elseif compact == "k" || occursin("nmax", compact) || occursin("ymax", compact) ||
           occursin("xmax", compact) || occursin("maxod", compact) ||
           occursin("carrying", compact) || occursin("capacity", compact) ||
           occursin("plateau", compact) || occursin("asymptote", compact)
        return features[:plateau]
    elseif occursin("amplitude", compact) || compact == "amp"
        return features[:amplitude]
    elseif compact in ("t0", "tmid", "tmax", "tinf", "tinflection") ||
           occursin("inflection", compact) || occursin("midtime", compact)
        return features[:mid_time]
    elseif compact in ("q0", "h0")
        return features[:q0]
    end
    return 1.0
end

function _batch_initial_params(model_name::String, model, t, y)
    features = _batch_growth_features(t, y)
    if haskey(BATCH_MODEL_INIT, model_name)
        init = BATCH_MODEL_INIT[model_name](features)
        length(init) == length(model.param_names) &&
            return [_batch_clip_initial(Float64(v)) for v in init]
    end
    return [_batch_clip_initial(_batch_initial_value(name, features)) for name in model.param_names]
end

function _batch_param_bounds(model, t, y)
    features = _batch_growth_features(t, y)
    lower = Float64[]
    upper = Float64[]
    plateau_cap = max(features[:plateau] * 3, 1.0)
    time_cap = max(features[:duration] * 2, 10.0)
    for name in model.param_names
        compact = replace(lowercase(String(name)), r"[^a-z0-9]" => "")
        if occursin("max", compact) || occursin("capacity", compact) || compact == "k" ||
           occursin("plateau", compact) || occursin("asymptote", compact)
            push!(lower, 1e-6); push!(upper, plateau_cap)
        elseif occursin("lag", compact) || occursin("lambda", compact) || occursin("delay", compact)
            push!(lower, 0.0); push!(upper, max(time_cap, 50.0))
        elseif occursin("growth", compact) || occursin("rate", compact) ||
               occursin("mu", compact) || compact in ("r", "gr")
            push!(lower, 1e-6); push!(upper, 20.0)
        elseif occursin("shape", compact) || compact in ("nu", "theta", "beta", "gamma", "m", "v")
            push!(lower, 1e-6); push!(upper, 10.0)
        elseif occursin("slope", compact) || occursin("linear", compact)
            push!(lower, -10.0); push!(upper, 10.0)
        else
            push!(lower, 0.0); push!(upper, 50.0)
        end
    end
    return lower, upper
end

function _batch_linear_interp(x, xs, ys)
    isempty(xs) && return NaN
    x <= first(xs) && return Float64(first(ys))
    x >= last(xs) && return Float64(last(ys))
    for i in 1:(length(xs)-1)
        if xs[i] <= x <= xs[i+1]
            dx = xs[i+1] - xs[i]
            dx == 0 && return Float64(ys[i])
            return Float64(ys[i] + (x - xs[i]) / dx * (ys[i+1] - ys[i]))
        end
    end
    return Float64(last(ys))
end

function _batch_loss_rmse(t, y, fit_t, fit_y, score_end)
    isempty(fit_t) && return Inf
    ss = 0.0
    n = 0
    for (i, ti) in enumerate(t)
        (ti < first(fit_t) || ti > min(score_end, last(fit_t))) && continue
        pred = _batch_linear_interp(ti, fit_t, fit_y)
        isfinite(pred) || continue
        ss += (y[i] - pred)^2
        n += 1
    end
    return n == 0 ? Inf : sqrt(ss / n)
end

function _batch_prepare_curve(
    od_raw::Vector{Float64};
    blank_value::Real=0.0,
    subtract_blank::Bool=false,
    blank_method::String="pointbypoint",
    blank_timeseries::Vector{Float64}=Float64[],
    unblanked_floor::Float64=0.01,
)
    if !(subtract_blank && blank_value > 0.0)
        return (
            od_for_fit=max.(od_raw, unblanked_floor),
            od_subtracted_display=nothing,
            anchor=od_raw[1],
            shift=0.0,
        )
    end

    pointwise = blank_method == "pointbypoint" && length(blank_timeseries) == length(od_raw)
    blank_trace = pointwise ? blank_timeseries : fill(Float64(blank_value), length(od_raw))
    method = blank_method == "clip" ? :clip : blank_method == "shift" ? :shift : :pointbypoint
    corrected = pointwise ? od_raw .- blank_trace : od_raw .- blank_value
    od_for_fit = vec(apply_blank_timeseries(
        reshape(od_raw, 1, :), blank_trace; method, floor=1e-4,
    ))

    return (
        od_for_fit=od_for_fit,
        od_subtracted_display=max.(corrected, 0.0),
        anchor=max(corrected[1], 0.0),
        shift=method == :clip ? 0.0 : od_for_fit[1] - corrected[1],
    )
end

function _batch_run_attempt(
    optimizer::String,
    time_numeric::Vector{Float64},
    od_for_fit::Vector{Float64},
    shift::Float64,
    subtract_blank::Bool,
    blank_value::Real,
    label::String,
    model_name::String,
    model_names::Vector{String},
    maxiters::Int,
    abstol::Float64,
    smooth::Bool,
    smooth_window::Int,
    optimizer_seed::Int,
)
    time_fit = time_numeric
    od_fit = od_for_fit

    model_keys = isempty(model_names) ? [model_name] : filter(!=("log_lin"), model_names)
    isempty(model_keys) && error("No parametric models selected")
    models = [MODEL_REGISTRY[k] for k in model_keys]
    initial_params = [_batch_initial_params(model_keys[i], models[i], time_fit, od_fit) for i in eachindex(models)]
    bounds = [_batch_param_bounds(m, time_fit, od_fit) for m in models]
    lower = Union{Nothing, Vector{Float64}}[b[1] for b in bounds]
    upper = Union{Nothing, Vector{Float64}}[b[2] for b in bounds]
    spec = ModelSpec(models, initial_params; lower=lower, upper=upper)
    opt_params = abstol > 0.0 ? (maxiters=maxiters, abstol=abstol) : (maxiters=maxiters,)
    opts = FitOptions(
        scattering_correction=false,
        smooth=smooth,
        smooth_method=:boxcar,
        boxcar_window=smooth_window,
        cut_stationary_phase=true,
        stationary_percentile_thr=0.05,
        stationary_pt_smooth_derivative=10,
        stationary_win_size=5,
        loss="RE",
        optimizer=_batch_optimizer(optimizer),
        optimizer_seed=optimizer_seed,
        opt_params=opt_params,
    )

    fit_results = kinbiont_fit(GrowthData(reshape(od_fit, 1, :), time_fit, [label]), spec, opts)
    r = fit_results[1]
    preprocessed_time = Float64.(fit_results.data.times)
    preprocessed_od = Float64.(vec(fit_results.data.curves[1, :]))
    fit_od_curve = subtract_blank && blank_value > 0.0 ? r.fitted_curve .- shift : r.fitted_curve
    fit_time_out = collect(r.times)
    fit_od_out = collect(fit_od_curve)

    stationary_phase_start = Float64(last(r.times))
    preprocessed_od_out = subtract_blank && blank_value > 0.0 ? preprocessed_od .- shift : preprocessed_od
    return (
        best_params=r.best_params,
        param_names=r.best_model.param_names,
        model_name=_batch_model_name(r.best_model),
        fit_time_out=fit_time_out,
        fit_od_out=fit_od_out,
        preprocessed_time=preprocessed_time,
        preprocessed_od=preprocessed_od_out,
        stationary_phase_start=stationary_phase_start,
        aic=r.best_aic,
        loss_rmse=_batch_loss_rmse(
            preprocessed_time,
            preprocessed_od,
            Float64.(r.times),
            Float64.(r.fitted_curve),
            stationary_phase_start,
        ),
        loss_re=r.loss,
    )
end

"""
    loglin_stationary_nmax(raw; pt_smoothing_derivative=7) -> Float64

Return the empirical maximum OD at the stationary cutoff associated with the
log-linear `mu_max` stored in `raw`. The legacy q95 value in `raw[2][16]` is
not used.
"""
function loglin_stationary_nmax(
    raw;
    pt_smoothing_derivative::Int=7,
)::Float64
    params = raw[2]
    length(params) >= 7 && params[7] !== missing || return NaN
    smoothed = raw[4]
    ismissing(smoothed) && return NaN

    opts = FitOptions(
        stationary_percentile_thr=0.05,
        stationary_pt_smooth_derivative=pt_smoothing_derivative,
        stationary_win_size=5,
        stationary_thr_od=0.02,
    )
    cutoff = find_stationary_cutoff_from_mu(
        Matrix{Float64}(smoothed), Float64(params[7]), opts,
    )
    nmax = Float64(smoothed[2, cutoff])
    return isfinite(nmax) ? nmax : NaN
end

function _batch_loglin_fields(t, y, label, experiment; pt_avg=7, pt_deriv=7, pt_min_win=7, threshold=0.9)
    out = Dict{String, Any}(
        "gr_loglin" => NaN,
        "gr_loglin_se" => NaN,
        "gr_max_sliding" => NaN,
        "t_exp_start_loglin" => NaN,
        "t_exp_end_loglin" => NaN,
        "doubling_time_loglin" => NaN,
        "R_squared_loglin" => NaN,
        "lag_loglin" => NaN,
        "N_max_emp" => NaN,
        "loglin_converged" => false,
    )
    length(y) >= max(10, pt_deriv + pt_min_win + 2) || return out
    try
        raw = fitting_one_well_Log_Lin(
            Matrix(transpose(hcat(t, y))), label, experiment;
            type_of_smoothing="rolling_avg",
            pt_avg=pt_avg,
            pt_smoothing_derivative=pt_deriv,
            pt_min_size_of_win=pt_min_win,
            type_of_win="maximum",
            threshold_of_exp=threshold,
        )
        params = raw[2]
        if length(params) >= 14 && params[7] !== missing
            out["gr_loglin"] = Float64(params[7])
            out["gr_loglin_se"] = Float64(params[8])
            out["gr_max_sliding"] = Float64(params[6])
            out["t_exp_start_loglin"] = Float64(params[3])
            out["t_exp_end_loglin"] = Float64(params[4])
            out["doubling_time_loglin"] = Float64(params[9])
            out["R_squared_loglin"] = Float64(params[14])^2
            if length(params) >= 15
                out["lag_loglin"] = params[15] === missing ? NaN : Float64(params[15])
            end
            out["N_max_emp"] = loglin_stationary_nmax(
                raw; pt_smoothing_derivative=pt_deriv,
            )
            out["loglin_converged"] = true
        end
    catch
    end
    return out
end

function _batch_fit_one(
    time_numeric::Vector{Float64},
    od_raw::Vector{Float64},
    label::String,
    experiment::String;
    blank_value::Real=0.0,
    subtract_blank::Bool=false,
    blank_method::String="pointbypoint",
    blank_timeseries::Vector{Float64}=Float64[],
    model_name::String="aHPM",
    model_names::Vector{String}=String[],
    optimizer::String="LN_BOBYQA",
    deterministic_optimizers::Vector{String}=String[],
    stochastic_optimizers::Vector{String}=String[],
    stochastic_runs::Int=1,
    optimizer_seed::Int=42,
    maxiters::Int=100000,
    abstol::Float64=1e-15,
    smooth::Bool=false,
    smooth_window::Int=3,
    compute_loglin::Bool=false,
    loglin_pt_avg::Int=7,
    loglin_pt_smoothing_derivative::Int=7,
    loglin_pt_min_size_of_win::Int=7,
    loglin_threshold_of_exp::Float64=0.9,
)
    if smooth && (smooth_window < 3 || iseven(smooth_window))
        throw(ArgumentError("smooth_window must be an odd integer greater than or equal to 3"))
    end

    prepared = _batch_prepare_curve(
        od_raw;
        blank_value,
        subtract_blank,
        blank_method,
        blank_timeseries,
        unblanked_floor=0.01,
    )
    od_for_fit = prepared.od_for_fit
    od_subtracted_display = prepared.od_subtracted_display
    shift = prepared.shift

    attempts = Tuple{String, Int}[]
    if !isempty(deterministic_optimizers) || !isempty(stochastic_optimizers)
        append!(attempts, [(opt, 1) for opt in deterministic_optimizers])
        for opt in stochastic_optimizers, run_idx in 1:max(1, stochastic_runs)
            push!(attempts, (opt, run_idx))
        end
    else
        push!(attempts, (optimizer, 1))
    end
    isempty(attempts) && error("No optimizers selected")

    outcomes = NamedTuple[]
    for (opt, run_idx) in attempts
        attempt_seed = opt in stochastic_optimizers || (isempty(stochastic_optimizers) &&
            opt in ("GN_ISRES", "BBO_adaptive_de_rand_1_bin_radiuslimited")) ?
            optimizer_seed + run_idx - 1 : optimizer_seed
        try
            res = _batch_run_attempt(
                opt, time_numeric, od_for_fit, shift,
                subtract_blank, blank_value, label, model_name, model_names,
                maxiters, abstol, smooth, smooth_window, attempt_seed,
            )
            push!(outcomes, (optimizer=opt, run=run_idx, seed=attempt_seed, status="ok", loss=res.loss_rmse,
                loss_rmse=res.loss_rmse, loss_re=res.loss_re, aic=res.aic, result=res))
        catch e
            push!(outcomes, (optimizer=opt, run=run_idx, seed=attempt_seed, status="error: $(string(e))",
                loss=Inf, loss_rmse=Inf, loss_re=NaN, aic=NaN, result=nothing))
        end
    end
    successful = filter(o -> o.result !== nothing && isfinite(o.loss_rmse), outcomes)
    isempty(successful) && error("All optimizer attempts failed")
    best = successful[argmin([o.loss_rmse for o in successful])]
    win = best.result

    result = Dict{String, Any}(
        "experiment" => experiment,
        "well" => label,
        "experimental_time" => time_numeric,
        "experimental_od" => od_raw,
        "fit_time" => win.fit_time_out,
        "fit_od" => win.fit_od_out,
        "parameters" => win.best_params,
        "param_names" => win.param_names,
        "model" => win.model_name,
        "blank_value" => blank_value,
        "blank_subtraction" => subtract_blank,
        "blank_method" => blank_method,
        "stationary_phase_start" => win.stationary_phase_start,
        "maxiters" => maxiters,
        "abstol" => abstol,
        "preprocessing" => Dict(
            "smooth" => smooth,
            "smooth_method" => smooth ? "boxcar" : "none",
            "smooth_window" => smooth_window,
            "cut_stationary_phase" => true,
            "stationary_percentile_thr" => 0.05,
            "stationary_pt_smooth_derivative" => 10,
            "stationary_win_size" => 5,
        ),
        "aic" => win.aic,
        "loss" => win.loss_rmse,
        "loss_rmse" => win.loss_rmse,
        "loss_re" => win.loss_re,
        "optimizer_used" => best.optimizer,
        "optimizer_run" => best.run,
        "optimizer_seed" => best.seed,
        "all_attempts" => [Dict("optimizer" => o.optimizer, "run" => o.run, "seed" => o.seed,
            "status" => o.status, "loss" => o.loss, "loss_rmse" => o.loss_rmse,
            "loss_re" => o.loss_re, "aic" => o.aic) for o in outcomes],
    )
    od_subtracted_display !== nothing && (result["experimental_od_subtracted"] = od_subtracted_display)
    if smooth
        result["smoothed_time"] = win.preprocessed_time
        result["smoothed_od"] = win.preprocessed_od
    end
    if compute_loglin
        merge!(result, _batch_loglin_fields(
            time_numeric, od_for_fit, label, experiment;
            pt_avg=loglin_pt_avg,
            pt_deriv=loglin_pt_smoothing_derivative,
            pt_min_win=loglin_pt_min_size_of_win,
            threshold=loglin_threshold_of_exp,
        ))
    end
    return result
end

"""
    kinbiont_batch_fit(data::GrowthData; kwargs...)

Run a GUIbiont-compatible batch fit on a `GrowthData` object. The returned
named tuple contains `results`, `skipped`, and `errors` dictionaries shaped
like GUIbiont's `/api/batch-fit` response.

The returned `model` field is `"multi"` when more than one parametric model
is compared (and `model_names` holds the list); when a single model is
fitted, `model` is its name and `model_names == [model]`.

Each result carries two loss fields: `loss_re` is the relative-error objective
minimized by Kinbiont, whereas `loss_rmse` is recomputed against the same
preprocessed observations through the Kinbiont stationary-phase cutoff.
GUIbiont selects the optimizer attempt with the lowest `loss_rmse`.

Set `smooth=true` to apply the same centered moving average to every optimizer
attempt before stationary-phase detection. `smooth_window` is the odd window
width and defaults to 3 (one point on either side of the current point).
"""
function kinbiont_batch_fit(
    data::GrowthData;
    experiment::String="experiment",
    labels::Vector{String}=String[],
    model_name::String="aHPM",
    model_names::Vector{String}=String[],
    optimizer::String="LN_BOBYQA",
    deterministic_optimizers::Vector{String}=String[],
    stochastic_optimizers::Vector{String}=String[],
    stochastic_runs::Int=1,
    optimizer_seed::Int=42,
    maxiters::Int=100000,
    abstol::Float64=1e-15,
    skip_flat_threshold::Float64=0.02,
    smooth::Bool=false,
    smooth_window::Int=3,
    compute_loglin::Bool=false,
    loglin_pt_avg::Int=7,
    loglin_pt_smoothing_derivative::Int=7,
    loglin_pt_min_size_of_win::Int=7,
    loglin_threshold_of_exp::Float64=0.9,
    blank_subtraction::Bool=false,
    blank_method::String="pointbypoint",
    blank_value::Real=0.0,
    blank_timeseries::Vector{Float64}=Float64[],
)
    if smooth && (smooth_window < 3 || iseven(smooth_window))
        throw(ArgumentError("smooth_window must be an odd integer greater than or equal to 3"))
    end

    selected = isempty(labels) ? data.labels : labels
    results = Dict{String, Any}[]
    skipped = Dict{String, Any}[]
    errors = String[]

    for label in selected
        idx = findfirst(==(label), data.labels)
        if idx === nothing
            push!(errors, "Well '$label' not found")
            continue
        end
        t = data.times
        yraw = vec(data.curves[idx, :])
        valid = findall(i -> isfinite(t[i]) && isfinite(yraw[i]), eachindex(t))
        if length(valid) < 10
            push!(errors, "$label: insufficient data points")
            continue
        end
        tv = Float64.(t[valid])
        yv = Float64.(yraw[valid])
        blank_ts = isempty(blank_timeseries) ? Float64[] : Float64.(blank_timeseries[valid])
        amp = maximum(yv) - minimum(yv)
        if skip_flat_threshold > 0.0 && amp < skip_flat_threshold
            push!(skipped, Dict{String, Any}(
                "well" => label,
                "amplitude" => amp,
                "reason" => "flat curve (amplitude $(round(amp, digits=4)) < threshold $(skip_flat_threshold))",
            ))
            continue
        end
        try
            push!(results, _batch_fit_one(
                tv, yv, label, experiment;
                blank_value=blank_value,
                subtract_blank=blank_subtraction,
                blank_method=blank_method,
                blank_timeseries=blank_ts,
                model_name=model_name,
                model_names=model_names,
                optimizer=optimizer,
                deterministic_optimizers=deterministic_optimizers,
                stochastic_optimizers=stochastic_optimizers,
                stochastic_runs=stochastic_runs,
                optimizer_seed=optimizer_seed,
                maxiters=maxiters,
                abstol=abstol,
                smooth=smooth,
                smooth_window=smooth_window,
                compute_loglin=compute_loglin,
                loglin_pt_avg=loglin_pt_avg,
                loglin_pt_smoothing_derivative=loglin_pt_smoothing_derivative,
                loglin_pt_min_size_of_win=loglin_pt_min_size_of_win,
                loglin_threshold_of_exp=loglin_threshold_of_exp,
            ))
        catch e
            push!(errors, "$label: $(string(e))")
        end
    end

    model_out = isempty(model_names) ? model_name : "multi"
    model_names_out = isempty(model_names) ? [model_name] : model_names
    return (
        experiment=experiment,
        model=model_out,
        model_names=model_names_out,
        smooth=smooth,
        smooth_window=smooth_window,
        results=results,
        skipped=skipped,
        errors=errors,
    )
end

"""
    save_batch_results(batch, dir="."; prefix=nothing)

Write GUIbiont-compatible batch summary and fitted-curve CSV files. Returns
`(summary=..., fitted_curves=...)`.
"""
function save_batch_results(batch, dir::String="."; prefix::Union{Nothing, String}=nothing)
    mkpath(dir)
    prefix = prefix === nothing ? "$(batch.experiment)_batch_fit" : prefix
    results = batch.results

    param_headers = String[]
    for r in results
        for name in get(r, "param_names", String[])
            s = String(name)
            s in param_headers || push!(param_headers, s)
        end
    end

    any_loglin = any(r -> get(r, "loglin_converged", false) === true ||
                          (haskey(r, "gr_loglin") && isfinite(r["gr_loglin"])), results)
    loglin_headers = any_loglin ? ["gr_loglin", "gr_loglin_se", "gr_max_sliding",
        "t_exp_start_loglin", "t_exp_end_loglin", "doubling_time_loglin",
        "R_squared_loglin", "lag_loglin", "N_max_emp", "loglin_converged"] : String[]

    summary_header = Any["experiment", "well", "model", param_headers..., loglin_headers...,
        "stationary_phase_start", "aic", "loss_rmse", "loss_re", "optimizer_used"]
    summary_rows = Vector{Vector{Any}}()
    for r in results
        pmap = Dict(String(r["param_names"][i]) => r["parameters"][i] for i in eachindex(r["param_names"]))
        row = Any[batch.experiment, r["well"], r["model"]]
        append!(row, [get(pmap, h, "") for h in param_headers])
        append!(row, [get(r, h, "") for h in loglin_headers])
        append!(row, [get(r, "stationary_phase_start", ""), get(r, "aic", ""),
                      get(r, "loss_rmse", ""), get(r, "loss_re", ""),
                      get(r, "optimizer_used", "")])
        push!(summary_rows, row)
    end
    summary_path = joinpath(dir, prefix * ".csv")
    _batch_write_csv(summary_path, summary_header, summary_rows)

    curve_times = [_batch_format_time.(r["fit_time"]) for r in results if haskey(r, "fit_time")]
    time_labels = isempty(curve_times) ? String[] :
        sort(unique(vcat(curve_times...)); by=x -> parse(Float64, x))
    curve_header = Any["experiment", "well", "model", ["t_" * t for t in time_labels]...]
    curve_rows = Vector{Vector{Any}}()
    for r in results
        by_time = Dict(_batch_format_time(r["fit_time"][i]) => r["fit_od"][i] for i in eachindex(r["fit_time"]))
        push!(curve_rows, Any[batch.experiment, r["well"], r["model"], [get(by_time, t, "") for t in time_labels]...])
    end
    fitted_path = joinpath(dir, prefix * "_fitted_curves.csv")
    _batch_write_csv(fitted_path, curve_header, curve_rows)

    return (summary=summary_path, fitted_curves=fitted_path)
end

function _batch_fit_loglin_one(
    time_numeric::Vector{Float64},
    od_raw::Vector{Float64},
    label::String,
    experiment::String;
    blank_value::Real=0.0,
    subtract_blank::Bool=false,
    blank_method::String="pointbypoint",
    blank_timeseries::Vector{Float64}=Float64[],
    type_of_smoothing::String="rolling_avg",
    pt_avg::Int=7,
    pt_smoothing_derivative::Int=7,
    pt_min_size_of_win::Int=7,
    type_of_win::String="maximum",
    threshold_of_exp::Float64=0.9,
    start_exp_win_thr::Float64=0.05,
    thr_lowess::Float64=0.05,
    unblanked_floor::Float64=1e-4,
)
    prepared = _batch_prepare_curve(
        od_raw;
        blank_value,
        subtract_blank,
        blank_method,
        blank_timeseries,
        unblanked_floor,
    )
    od_for_fit = prepared.od_for_fit
    od_subtracted_display = prepared.od_subtracted_display

    raw = fitting_one_well_Log_Lin(
        Matrix(transpose(hcat(time_numeric, od_for_fit))),
        label,
        experiment;
        type_of_smoothing=type_of_smoothing,
        pt_avg=pt_avg,
        pt_smoothing_derivative=pt_smoothing_derivative,
        pt_min_size_of_win=pt_min_size_of_win,
        type_of_win=type_of_win,
        threshold_of_exp=threshold_of_exp,
        start_exp_win_thr=start_exp_win_thr,
        thr_lowess=thr_lowess,
    )
    params = raw[2]
    result = Dict{String, Any}(
        "experiment" => experiment,
        "well" => label,
        "method" => "Log-lin",
        "experimental_time" => time_numeric,
        "experimental_od" => od_raw,
        "blank_value" => blank_value,
        "blank_subtraction" => subtract_blank,
        "blank_method" => blank_method,
    )
    if length(params) >= 14 && params[7] !== missing
        result["gr_loglin"] = Float64(params[7])
        result["gr_loglin_se"] = Float64(params[8])
        result["gr_max_sliding"] = Float64(params[6])
        result["t_exp_start_loglin"] = Float64(params[3])
        result["t_exp_end_loglin"] = Float64(params[4])
        result["doubling_time_loglin"] = Float64(params[9])
        result["R_squared_loglin"] = Float64(params[14])^2
        result["lag_loglin"] = length(params) >= 15 && params[15] !== missing ? Float64(params[15]) : NaN
        # Deliberately ignore legacy params[16] (curve q95).
        result["N_max_emp"] = loglin_stationary_nmax(
            raw; pt_smoothing_derivative,
        )
        result["loglin_converged"] = true
    else
        for name in ["gr_loglin", "gr_loglin_se", "gr_max_sliding",
                     "t_exp_start_loglin", "t_exp_end_loglin",
                     "doubling_time_loglin", "R_squared_loglin",
                     "lag_loglin", "N_max_emp"]
            result[name] = NaN
        end
        result["loglin_converged"] = false
    end
    od_subtracted_display !== nothing && (result["experimental_od_subtracted"] = od_subtracted_display)
    return result
end

"""
    kinbiont_fit_loglin(data::GrowthData; kwargs...) -> Dict

Fit one growth curve with the log-linear method and return the same named
fields used by GUIbiont, including the stationary-cutoff `N_max_emp`.
"""
function kinbiont_fit_loglin(
    data::GrowthData;
    experiment::String="experiment",
    label::Union{Nothing, String}=nothing,
    blank_subtraction::Bool=false,
    blank_method::String="pointbypoint",
    blank_value::Real=0.0,
    blank_timeseries::Vector{Float64}=Float64[],
    unblanked_floor::Float64=1e-4,
    type_of_smoothing::String="rolling_avg",
    pt_avg::Int=7,
    pt_smoothing_derivative::Int=7,
    pt_min_size_of_win::Int=7,
    type_of_win::String="maximum",
    threshold_of_exp::Float64=0.9,
    start_exp_win_thr::Float64=0.05,
    thr_lowess::Float64=0.05,
)
    selected_label = isnothing(label) ? only(data.labels) : label
    idx = findfirst(==(selected_label), data.labels)
    idx === nothing && throw(ArgumentError("Well '$selected_label' not found"))

    t = data.times
    y = vec(data.curves[idx, :])
    valid = findall(i -> isfinite(t[i]) && isfinite(y[i]), eachindex(t))
    min_pts = max(10, pt_smoothing_derivative + pt_min_size_of_win + 2)
    length(valid) >= min_pts || throw(ArgumentError("Insufficient data points"))

    blank_ts = isempty(blank_timeseries) ? Float64[] : Float64.(blank_timeseries[valid])
    return _batch_fit_loglin_one(
        Float64.(t[valid]), Float64.(y[valid]), selected_label, experiment;
        blank_value,
        subtract_blank=blank_subtraction,
        blank_method,
        blank_timeseries=blank_ts,
        unblanked_floor,
        type_of_smoothing,
        pt_avg,
        pt_smoothing_derivative,
        pt_min_size_of_win,
        type_of_win,
        threshold_of_exp,
        start_exp_win_thr,
        thr_lowess,
    )
end

"""
    kinbiont_fit_loglin(data, preprocess_opts::FitOptions; kwargs...) -> Dict

Convenience overload that reuses the blank-correction settings from an
existing `FitOptions` object.
"""
function kinbiont_fit_loglin(
    data::GrowthData,
    preprocess_opts::FitOptions;
    kwargs...,
)
    blank_ts = isnothing(preprocess_opts.blank_timeseries) ?
        Float64[] : preprocess_opts.blank_timeseries
    return kinbiont_fit_loglin(
        data;
        blank_subtraction=preprocess_opts.blank_subtraction,
        blank_method=String(preprocess_opts.blank_method),
        blank_value=preprocess_opts.blank_value,
        blank_timeseries=blank_ts,
        unblanked_floor=preprocess_opts.negative_threshold,
        kwargs...,
    )
end

"""
    kinbiont_batch_loglin(data::GrowthData; kwargs...)

Run GUIbiont-compatible batch log-linear fitting on a `GrowthData` object.
"""
function kinbiont_batch_loglin(
    data::GrowthData;
    experiment::String="experiment",
    labels::Vector{String}=String[],
    blank_subtraction::Bool=false,
    blank_method::String="pointbypoint",
    blank_value::Real=0.0,
    blank_timeseries::Vector{Float64}=Float64[],
    type_of_smoothing::String="rolling_avg",
    pt_avg::Int=7,
    pt_smoothing_derivative::Int=7,
    pt_min_size_of_win::Int=7,
    type_of_win::String="maximum",
    threshold_of_exp::Float64=0.9,
    start_exp_win_thr::Float64=0.05,
    thr_lowess::Float64=0.05,
    skip_flat_threshold::Float64=0.02,
)
    selected = isempty(labels) ? data.labels : labels
    results = Dict{String, Any}[]
    skipped = Dict{String, Any}[]
    errors = String[]
    min_pts = max(10, pt_smoothing_derivative + pt_min_size_of_win + 2)

    for label in selected
        idx = findfirst(==(label), data.labels)
        if idx === nothing
            push!(errors, "Well '$label' not found")
            continue
        end
        t = data.times
        yraw = vec(data.curves[idx, :])
        valid = findall(i -> isfinite(t[i]) && isfinite(yraw[i]), eachindex(t))
        if length(valid) < min_pts
            push!(errors, "$label: insufficient data points")
            continue
        end
        tv = Float64.(t[valid])
        yv = Float64.(yraw[valid])
        amp = maximum(yv) - minimum(yv)
        if skip_flat_threshold > 0.0 && amp < skip_flat_threshold
            push!(skipped, Dict{String, Any}(
                "well" => label,
                "amplitude" => amp,
                "reason" => "flat curve (amplitude $(round(amp, digits=4)) < threshold $(skip_flat_threshold))",
            ))
            continue
        end
        blank_ts = isempty(blank_timeseries) ? Float64[] : Float64.(blank_timeseries[valid])
        try
            push!(results, _batch_fit_loglin_one(
                tv, yv, label, experiment;
                blank_value=blank_value,
                subtract_blank=blank_subtraction,
                blank_method=blank_method,
                blank_timeseries=blank_ts,
                type_of_smoothing=type_of_smoothing,
                pt_avg=pt_avg,
                pt_smoothing_derivative=pt_smoothing_derivative,
                pt_min_size_of_win=pt_min_size_of_win,
                type_of_win=type_of_win,
                threshold_of_exp=threshold_of_exp,
                start_exp_win_thr=start_exp_win_thr,
                thr_lowess=thr_lowess,
            ))
        catch e
            push!(errors, "$label: $(string(e))")
        end
    end

    return (
        experiment=experiment,
        model="log_lin",
        model_names=["log_lin"],
        results=results,
        skipped=skipped,
        errors=errors,
    )
end

"""
    save_batch_loglin_results(batch, dir="."; prefix=nothing)

Write GUIbiont-compatible `batch_fit_loglin` CSV output.
"""
function save_batch_loglin_results(batch, dir::String="."; prefix::Union{Nothing, String}=nothing)
    mkpath(dir)
    prefix = prefix === nothing ? "$(batch.experiment)_batch_fit_loglin" : prefix
    header = Any["experiment", "well", "method",
        "gr_loglin", "gr_loglin_se", "gr_max_sliding",
        "t_exp_start_loglin", "t_exp_end_loglin",
        "doubling_time_loglin", "R_squared_loglin",
        "lag_loglin", "N_max_emp", "loglin_converged"]
    rows = Vector{Vector{Any}}()
    for r in batch.results
        push!(rows, Any[
            batch.experiment, r["well"], "log_lin",
            get(r, "gr_loglin", ""), get(r, "gr_loglin_se", ""),
            get(r, "gr_max_sliding", ""), get(r, "t_exp_start_loglin", ""),
            get(r, "t_exp_end_loglin", ""), get(r, "doubling_time_loglin", ""),
            get(r, "R_squared_loglin", ""), get(r, "lag_loglin", ""),
            get(r, "N_max_emp", ""),
            get(r, "loglin_converged", false) ? "true" : "false",
        ])
    end
    path = joinpath(dir, prefix * ".csv")
    _batch_write_csv(path, header, rows)
    return (summary=path,)
end

export kinbiont_batch_fit
export save_batch_results
export loglin_stationary_nmax
export kinbiont_fit_loglin
export kinbiont_batch_loglin
export save_batch_loglin_results
