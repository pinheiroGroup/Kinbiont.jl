# [Fitting](@id fitting)

```@contents
Pages = ["index.md"]
Depth = 3
```

To run all examples on this page:

```julia
using Kinbiont, CSV, DataFrames, Plots, Random
```

---

## The three model types

Kinbiont dispatches on the type of model object:

| Type | What it is | Entry in `MODEL_REGISTRY` |
|---|---|---|
| `ODEModel` | Mechanistic ODE in SciML in-place form `f!(du, u, p, t)` | Yes (`"logistic"`, `"aHPM"`, …) |
| `NLModel` | Closed-form algebraic equation fitted to OD directly | Yes (`"NL_Gompertz"`, …) |
| `LogLinModel` | Log-linear fit to exponential phase; no parameters needed | No (use `LogLinModel()`) |

---

## Discovering available models

```julia
# List all registered model names
sort(collect(keys(MODEL_REGISTRY)))

# Inspect one model
m = MODEL_REGISTRY["aHPM"]
println(m.name)         # "aHPM"
println(m.param_names)  # ["gr", "exit_lag_rate", "N_max", "shape"]

m2 = MODEL_REGISTRY["NL_Gompertz"]
println(m2.param_names) # ["N_max", "growth_rate", "lag"]

m3 = MODEL_REGISTRY["logistic"]
println(m3.param_names) # ["gr", "N_max"]
```

---

## Fitting a single curve

```julia
using Kinbiont, Random
Random.seed!(1)

# Simulate one curve
sim  = ODE_sim("aHPM", [0.05], 0.0, 48.0, 1.0, [0.4, 0.1, 1.2, 2.0])
times = Float64.(sim.t)
od    = Float64.(reduce(hcat, sim.u)[1, :]) .+ 0.005 .* randn(length(times))
data  = GrowthData(reshape(od, 1, length(od)), times, ["sim_A1"])

spec = ModelSpec(
    [MODEL_REGISTRY["aHPM"]],
    [[1.0, 1.0, 1.0, 1.0]];
    lower = [[0.0, 0.0, 0.0, 0.0]],
    upper = [[50.0, 50.0, 50.0, 50.0]],
)

results = kinbiont_fit(data, spec, FitOptions(smooth=true))
r = results[1]   # CurveFitResult

println("Model:  ", r.best_model.name)
println("Params: ", r.best_params)   # [gr, exit_lag_rate, N_max, shape]
println("AICc:   ", round(r.best_aic; digits=2))
println("Loss:   ", round(r.loss; sigdigits=3))

# Plot
scatter(times, od; label="data", markersize=2, alpha=0.6,
        xlabel="Time (h)", ylabel="OD")
plot!(r.times, r.fitted_curve; label="aHPM fit", linewidth=2)
```

---

## Model selection

Pass multiple models in one `ModelSpec`. Kinbiont fits all candidates and
selects the best by AICc automatically.

```julia
spec = ModelSpec(
    [MODEL_REGISTRY["aHPM"], MODEL_REGISTRY["logistic"], MODEL_REGISTRY["NL_Gompertz"]],
    [
        [1.0, 1.0, 1.0, 1.0],   # aHPM: [gr, exit_lag_rate, N_max, shape]
        [0.5, 1.0],              # logistic: [gr, N_max]
        [1.0, 0.5, 5.0],        # NL_Gompertz: [N_max, growth_rate, lag]
    ];
    lower = [[0.0,0.0,0.0,0.0], [0.0,0.0], [0.0,0.0,0.0]],
    upper = [[50.,50.,50.,50.], [10.,10.], [10.,10.,100.]],
)

results = kinbiont_fit(data, spec, FitOptions(smooth=true))
r = results[1]
println("Best model: ", r.best_model.name, "  AICc=", round(r.best_aic; digits=2))

# Inspect all candidates
for c in r.all_results
    println("  ", rpad(c.model_name, 16), " AICc=", round(c.aic; digits=2))
end
```

---

## Bounds and initial guesses

The default optimizer (`BBO_adaptive_de_rand_1_bin_radiuslimited`) requires box
bounds. If you omit them, Kinbiont warns and uses `param_guess/10` to `param_guess*10`.

Rules of thumb:
- Lower bounds for rates and concentrations: `0.0`
- Upper bound for growth rates: 10–50 (per hour; depends on organism)
- Upper bound for N_max: 2–5× the maximum OD you expect

```julia
spec = ModelSpec(
    [MODEL_REGISTRY["aHPM"]],
    [[0.5, 0.5, 1.0, 2.0]];   # [gr, exit_lag_rate, N_max, shape]
    lower = [[0.001, 0.0, 0.01, 0.1]],
    upper = [[5.0,   5.0, 5.0,  10.0]],
)
```

To use an optimizer that does not require bounds (e.g. `NOMADOpt()`):

```julia
using OptimizationNOMAD
opts = FitOptions(optimizer=NOMADOpt())
spec = ModelSpec([MODEL_REGISTRY["aHPM"]], [[0.5, 0.5, 1.0, 2.0]])   # no bounds
```

---

## Fitting a whole plate

```julia
using Kinbiont, CSV

data = GrowthData("data_examples/ecoli_sample.csv")

spec = ModelSpec(
    [MODEL_REGISTRY["aHPM"]],
    [[1.0, 1.0, 1.0, 1.0]];
    lower = [[0.0, 0.0, 0.0, 0.0]],
    upper = [[50.0, 50.0, 50.0, 50.0]],
)

opts    = FitOptions(smooth=true, smooth_method=:rolling_avg, smooth_pt_avg=7)
results = kinbiont_fit(data, spec, opts)

# Print summary
for r in results
    p = round.(Float64.(r.best_params); digits=3)
    println(r.label, "  gr=", p[1], "  N_max=", p[3],
            "  AICc=", round(r.best_aic; digits=1))
end

# Save three CSV files
paths = save_results(results, "output/plate/"; prefix="ecoli")
println("Summary:       ", paths.summary)
println("Fitted curves: ", paths.fitted_curves)
println("All models:    ", paths.all_models)
```

The `_summary.csv` columns are:
`label, cluster, best_model, n_params, param_names, aic, loss, param_1, param_2, …`

---

## Multistart

For difficult landscapes with many local minima, run k-means-style multistart:

```julia
opts = FitOptions(
    smooth      = true,
    multistart  = true,
    n_restart   = 50,   # run optimizer 50 times with TikTak initialisation
)
results = kinbiont_fit(data, spec, opts)
```

Multistart significantly increases runtime. Use when single-start fits converge
to implausible parameter values.

---

## Subsetting and re-fitting

```julia
data = GrowthData("data_examples/ecoli_sample.csv")

# Fit only two wells
sub = data[["Curve11730", "Curve11731"]]
results = kinbiont_fit(sub, spec, opts)

# After clustering, fit only growing wells
clustered = kinbiont_fit(data, FitOptions(cluster=true, n_clusters=4,
                                          cluster_prescreen_constant=true))
d = clustered.data
growing = d[d.labels[d.clusters .!= 4]]
results = kinbiont_fit(growing, spec, opts)
```

---

## Log-linear fitting

`LogLinModel` fits a straight line to the log-transformed data in the exponential
phase. It requires no parameter guess and is fast.

```julia
spec_loglin = ModelSpec(
    [LogLinModel()],
    [Float64[]],    # no parameters needed
)

results = kinbiont_fit(data, spec_loglin, FitOptions())
r = results[1]
println("slope (≈ growth rate): ", r.best_params[2])
```

`LogLinModel` can be combined with other models in one `ModelSpec` for automatic
comparison:

```julia
spec = ModelSpec(
    [LogLinModel(), MODEL_REGISTRY["aHPM"]],
    [Float64[], [1.0, 1.0, 1.0, 1.0]];
    lower = [nothing, [0.0, 0.0, 0.0, 0.0]],
    upper = [nothing, [50., 50., 50., 50.]],
)
```

---

## Legacy API

Kinbiont also has a lower-level API: `fitting_one_well_ODE_constrained`,
`fitting_one_well_Log_Lin`, `segmentation_ODE_file`, etc. These functions
still work and are documented in the [API Reference](@ref). New code should
prefer `kinbiont_fit`.
