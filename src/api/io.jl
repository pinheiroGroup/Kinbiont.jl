# =============================================================================
# Results I/O
# =============================================================================
# save_results(results, dir) writes three CSV files:
#
#   <prefix>_summary.csv      — one row per curve: best model, AIC, loss, params
#   <prefix>_fitted_curves.csv — long-format: label, time, observed, fitted
#   <prefix>_all_models.csv   — one row per (curve × candidate model)
#
# All param columns are named param_1 … param_n and padded with `missing`
# when models have fewer parameters, so every file has a rectangular shape.
# =============================================================================

using CSV
using DataFrames

# LogLinModel has no .name / .param_names fields — provide fallbacks
_model_name(m::AbstractGrowthModel)        = m.name
_model_name(::LogLinModel)                 = "log_lin"
_model_param_names(m::AbstractGrowthModel) = m.param_names
_model_param_names(::LogLinModel)          = ["intercept", "slope"]

# Look up param names by model name string (for all_results which store name, not object).
const _EXTRA_PARAM_NAMES = Dict(
    "log_lin" => ["intercept", "slope"],
    "DDDE"    => String[],   # data-driven; parameter count varies
)
function _param_names_by_name(name::String)
    haskey(_EXTRA_PARAM_NAMES, name) && return _EXTRA_PARAM_NAMES[name]
    haskey(MODEL_REGISTRY, name)     && return MODEL_REGISTRY[name].param_names
    return String[]
end

"""
    save_results(results::GrowthFitResults, dir::String; prefix="kinbiont")

Write fitting results to CSV files inside `dir` (created if it does not exist).

Three files are produced:
- `<prefix>_summary.csv`: one row per curve — label, cluster id, best model
  name, AICc, loss, and one `param_k` column per parameter (padded with
  `missing` for models with fewer parameters than the widest result).
- `<prefix>_fitted_curves.csv`: long-format table with columns
  `label`, `time`, `observed`, `fitted`.
- `<prefix>_all_models.csv`: one row per (curve, candidate model) pair,
  useful for comparing AICc across all models tried.

# Example
```julia
results = kinbiont_fit(data, spec, opts)
save_results(results, "output/my_experiment")
# writes output/my_experiment_summary.csv, etc.
```
"""
function save_results(
    results::GrowthFitResults,
    dir::String;
    prefix::String = "kinbiont",
)
    mkpath(dir)
    base = joinpath(dir, prefix)
    _save_summary(results, base * "_summary.csv")
    _save_fitted_curves(results, base * "_fitted_curves.csv")
    _save_all_models(results, base * "_all_models.csv")
    return (
        summary       = base * "_summary.csv",
        fitted_curves = base * "_fitted_curves.csv",
        all_models    = base * "_all_models.csv",
    )
end

# ---------------------------------------------------------------------------
# 1. Summary table
# ---------------------------------------------------------------------------

function _save_summary(results::GrowthFitResults, path::String)
    n_max_params = maximum(length(r.best_params) for r in results.results)

    # Build DataFrame with explicit param columns so missing pads correctly
    df = DataFrame(
        label       = [r.label for r in results.results],
        cluster     = [
            results.data.clusters === nothing ? missing :
            results.data.clusters[findfirst(==(r.label), results.data.labels)]
            for r in results.results
        ],
        best_model  = [_model_name(r.best_model) for r in results.results],
        n_params    = [length(r.best_params) for r in results.results],
        param_names = [join(_model_param_names(r.best_model), ";") for r in results.results],
        aic         = [r.best_aic for r in results.results],
        loss        = [r.loss for r in results.results],
    )

    for k in 1:n_max_params
        col = [
            k <= length(r.best_params) ? r.best_params[k] : missing
            for r in results.results
        ]
        df[!, Symbol("param_$k")] = col
    end

    CSV.write(path, df)
    @info "Saved summary → $path"
end

# ---------------------------------------------------------------------------
# 2. Fitted curves (long format)
# ---------------------------------------------------------------------------

function _save_fitted_curves(results::GrowthFitResults, path::String)
    rows = DataFrame(
        label    = String[],
        time     = Float64[],
        observed = Float64[],
        fitted   = Float64[],
    )

    for r in results.results
        idx = findfirst(==(r.label), results.data.labels)
        obs = results.data.curves[idx, :]

        # fitted_curve may be shorter than times when negative values were removed
        n_fit = length(r.fitted_curve)
        for (j, t) in enumerate(r.times[1:n_fit])
            obs_j = j <= length(obs) ? obs[j] : missing
            push!(rows, (r.label, t, obs_j, r.fitted_curve[j]))
        end
    end

    CSV.write(path, rows)
    @info "Saved fitted curves → $path"
end

# ---------------------------------------------------------------------------
# 3. All candidate models per curve
# ---------------------------------------------------------------------------

function _save_all_models(results::GrowthFitResults, path::String)
    n_max_params = maximum(
        maximum(length(c.params) for c in r.all_results)
        for r in results.results
    )

    df = DataFrame(
        label       = String[],
        model_name  = String[],
        param_names = String[],
        aic         = Float64[],
        loss        = Float64[],
        is_best     = Bool[],
    )
    for k in 1:n_max_params
        df[!, Symbol("param_$k")] = Any[]
    end

    for r in results.results
        for c in r.all_results
            row = Dict{Symbol, Any}(
                :label       => r.label,
                :model_name  => c.model_name,
                :param_names => join(_param_names_by_name(c.model_name), ";"),
                :aic         => c.aic,
                :loss        => c.loss,
                :is_best     => c.model_name == _model_name(r.best_model),
            )
            for k in 1:n_max_params
                row[Symbol("param_$k")] = k <= length(c.params) ? c.params[k] : missing
            end
            push!(df, row)
        end
    end

    CSV.write(path, df)
    @info "Saved all-models comparison → $path"
end

export save_results
