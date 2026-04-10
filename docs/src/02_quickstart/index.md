# [Quick Start](@id quickstart)

```@contents
Pages = ["index.md"]
Depth = 3
```

Kinbiont's new API has three core objects and one entry point:

| Object | Role |
|---|---|
| `GrowthData` | Holds your curves, time points, and labels |
| `FitOptions` | All preprocessing + fitting settings (40+ options, all have defaults) |
| `ModelSpec` | Which model(s) to fit and their starting parameters |
| `kinbiont_fit` | Run everything: preprocess → fit → return results |

Choose your entry point below.

---

## Track 1 — I have a CSV file

The fastest path. Your CSV must have time in the first column and one well per remaining column:

```
Time_h,Curve11728,Curve11729,...
0.25,0.017,0.006,...
0.5,0.017,0.006,...
```

The example below uses a subset of the *E. coli* Keio knockout growth curve dataset (Reding-Roman et al., *Nature Scientific Data* 2026, doi:10.1038/s41597-026-07075-9).

```julia
using Kinbiont

# 1. Load ── first column = time, remaining columns = wells
data = GrowthData("data_examples/ecoli_sample.csv")
# data.curves  → 10 × 400 Matrix{Float64}  (n_curves × n_timepoints)
# data.times   → Vector{Float64} of length 400
# data.labels  → ["Curve11728", "Curve11729", …]

# 2. Specify model and starting parameters
#    MODEL_REGISTRY holds all built-in models; inspect with keys(MODEL_REGISTRY)
spec = ModelSpec(
    [MODEL_REGISTRY["aHPM"]],         # aHPM: asymmetric Huang-Pirt model
    [[1.0, 1.0, 1.0, 1.0]];           # initial guess: [gr, exit_lag_rate, N_max, shape]
    lower = [[0.0, 0.0, 0.0, 0.0]],
    upper = [[50.0, 50.0, 50.0, 50.0]],
)

# 3. Configure: smooth with rolling average before fitting
opts = FitOptions(smooth=true, smooth_method=:rolling_avg, smooth_pt_avg=7)

# 4. Fit all curves
results = kinbiont_fit(data, spec, opts)

# 5. Inspect
for r in results
    p = round.(Float64.(r.best_params); digits=3)
    println(r.label, " → ", r.best_model.name,
            "  gr=", p[1], "  N_max=", p[3],
            "  AICc=", round(r.best_aic; digits=1))
end

# 6. Save results to three CSV files
paths = save_results(results, "output/ecoli/"; prefix="knockout")
# → output/ecoli/knockout_summary.csv        (one row per curve)
# → output/ecoli/knockout_fitted_curves.csv  (long-format observed vs fitted)
# → output/ecoli/knockout_all_models.csv     (all candidate models compared)
```

!!! tip "Finding available models"
    ```julia
    collect(keys(MODEL_REGISTRY))          # list all model names
    MODEL_REGISTRY["aHPM"].param_names     # ["gr", "exit_lag_rate", "N_max", "shape"]
    ```

!!! note "Auto-bounds warning"
    If you omit `lower`/`upper`, Kinbiont warns and sets bounds as `guess/10` to `guess*10`. Always provide explicit bounds for better results.

!!! note "What `save_results` writes"
    - `_summary.csv`: one row per curve — label, cluster, best model name, AICc, loss, and one `param_k` column per parameter.
    - `_fitted_curves.csv`: long-format `label, time, observed, fitted`.
    - `_all_models.csv`: one row per (curve, candidate model) — useful for model comparison.

---

## Track 2 — I'm building from a DataFrame

Use this when your data is already in Julia memory (e.g. after filtering a DataFrame, merging plates, or computing derived quantities).

```julia
using Kinbiont, CSV, DataFrames

df = CSV.read("data_examples/ecoli_sample.csv", DataFrame)

times  = Float64.(df.Time_h)           # first column
labels = names(df)[2:end]              # well names from column headers
curves = Matrix{Float64}(df[:, 2:end])'  # transpose: n_curves × n_timepoints

data = GrowthData(curves, times, labels)
# → identical to GrowthData("data_examples/ecoli_sample.csv")
```

The transpose `'` is critical — `GrowthData` expects rows to be curves and columns to be time points. After this, use the same `kinbiont_fit` call as Track 1.

!!! tip "Filtering wells before fitting"
    ```julia
    keep = filter(l -> startswith(l, "Curve117"), labels)
    data = data[keep]   # GrowthData subsetting — returns a new GrowthData
    ```

---

## Track 3 — Full round-trip (simulate → save → reload → fit)

Self-contained example — no external files needed. Good for exploring the API or testing a new model.

```julia
using Kinbiont, CSV, DataFrames, Random
Random.seed!(42)

# ── 1. Simulate a noisy aHPM growth curve ────────────────────────────────────
#    params: gr=0.4, exit_lag_rate=0.1, N_max=1.2, shape=2.0
sim   = ODE_sim("aHPM", [0.05], 0.0, 48.0, 1.0, [0.4, 0.1, 1.2, 2.0])
times = Float64.(sim.t)
od    = Float64.(reduce(hcat, sim.u)[1, :]) .+ 0.005 .* randn(length(times))

# ── 2. Save to CSV ───────────────────────────────────────────────────────────
CSV.write("sim_growth.csv", DataFrame(Time_h=times, well_A1=od))

# ── 3. Reload via GrowthData CSV constructor ─────────────────────────────────
data = GrowthData("sim_growth.csv")

# ── 4. Fit ───────────────────────────────────────────────────────────────────
spec = ModelSpec(
    [MODEL_REGISTRY["aHPM"]],
    [[1.0, 1.0, 1.0, 1.0]];
    lower = [[0.0, 0.0, 0.0, 0.0]],
    upper = [[50.0, 50.0, 50.0, 50.0]],
)
results = kinbiont_fit(data, spec, FitOptions(smooth=true))

# ── 5. Inspect ───────────────────────────────────────────────────────────────
r = results[1]
println("Model:      ", r.best_model.name)
for (name, val) in zip(r.best_model.param_names, r.best_params)
    println("  ", rpad(name, 16), round(Float64(val); digits=4))
end
println("AICc:       ", round(r.best_aic; digits=2))
println("Loss:       ", round(r.loss; sigdigits=3))
```

---

## Next steps

- **Preprocessing in depth** — smoothing methods, blank subtraction, stationary phase → [Preprocessing](@ref preprocessing)
- **Clustering** — group curves by shape before fitting → [Clustering](@ref clustering)
- **All model types** — ODE, NL, log-linear, model selection → [Fitting](@ref fitting)
- **ML downstream** — decision tree and symbolic regression on fitted parameters → [ML Downstream](@ref ml)
