# [Data & Types](@id data)

```@contents
Pages = ["index.md"]
Depth = 3
```

This page describes the data types used throughout Kinbiont's new API.
For a runnable introduction, see [Quick Start](@ref quickstart).

---

## `GrowthData`

Immutable container for a set of growth curves.

```julia
# Fields
data.curves    # Matrix{Float64}: n_curves × n_timepoints
data.times     # Vector{Float64}: shared time points
data.labels    # Vector{String}: one label per curve
data.clusters  # Union{Nothing, Vector{Int}}: populated by preprocess()
data.centroids # Union{Nothing, Matrix{Float64}}: n_clusters × n_timepoints (z-scored)
data.wcss      # Union{Nothing, Float64}: within-cluster sum of squares
```

**Three constructors:**

```julia
# 1. From a CSV file (first column = time, remaining = wells)
data = GrowthData("experiment.csv")

# 2. From a matrix (n_curves × n_timepoints orientation)
curves = Matrix{Float64}(df[:, 2:end])'
data   = GrowthData(curves, times, labels)

# 3. After loading with DataFrames
df     = CSV.read("experiment.csv", DataFrame)
times  = Float64.(df[:, 1])
labels = names(df)[2:end]
curves = Matrix{Float64}(df[:, 2:end])'
data   = GrowthData(curves, times, labels)
```

**Subsetting:**

```julia
# Return a new GrowthData with only the listed wells
sub = data[["Curve11728", "Curve11729"]]
```

---

## CSV format

First column: time. Remaining columns: one curve per column, header = well name.

```
Time_h,A1,A2,A3
0.0,0.09,0.09,0.087
1.0,0.08,0.011,0.012
2.0,0.011,0.18,0.1
```

Separator must be a comma (`,`).

---

## Annotation format

Optional two-column CSV (no header). Maps well names to labels:
- `"b"` → blank well (used for `blank_from_labels=true`)
- `"X"` → discard (excluded from all analysis)
- Any other string → biological label; identical labels are treated as replicates

```
A1,b
A2,X
A3,condition_1
A4,condition_1
A5,condition_2
```

---

## `FitOptions`

All configuration in one `@kwdef` struct. Every field has a default — only set
what you need.

```julia
# Minimal: all defaults
opts = FitOptions()

# Common pattern
opts = FitOptions(
    smooth               = true,
    smooth_method        = :rolling_avg,
    blank_subtraction    = true,
    blank_from_labels    = true,
    cut_stationary_phase = true,
    cluster              = true,
    n_clusters           = 4,
    loss                 = "RE",
)
```

Field groups and their defaults:

| Group | Key fields | Defaults |
|---|---|---|
| Smoothing | `smooth`, `smooth_method`, `smooth_pt_avg`, `lowess_frac`, `boxcar_window`, `gaussian_h_mult` | `false`, `:lowess`, `7`, `0.05`, `5`, `2.0` |
| Blank subtraction | `blank_subtraction`, `blank_value`, `blank_from_labels` | `false`, `0.0`, `false` |
| Replicates | `average_replicates` | `false` |
| Negatives | `correct_negatives`, `negative_method`, `negative_threshold` | `false`, `:remove`, `0.01` |
| Scattering | `scattering_correction`, `calibration_file`, `scattering_method` | `false`, `""`, `:interpolation` |
| Stationary | `cut_stationary_phase`, `stationary_percentile_thr`, `stationary_win_size` | `false`, `0.05`, `5` |
| Clustering | `cluster`, `n_clusters`, `cluster_trend_test`, `cluster_prescreen_constant`, `cluster_exp_prototype`, `kmeans_seed` | `false`, `3`, `true`, `false`, `false`, `0` |
| Fitting | `loss`, `multistart`, `n_restart`, `optimizer`, `integrator` | `"RE"`, `false`, `50`, BBO, Tsit5 |

See the [API Reference](@ref) for the full field list.

---

## `ModelSpec`

Specifies which models to fit and their starting parameters.

```julia
spec = ModelSpec(
    models  = [MODEL_REGISTRY["aHPM"], MODEL_REGISTRY["logistic"]],
    params  = [[1.0, 1.0, 1.0, 1.0], [0.5, 1.0]];
    lower   = [[0.0,0.0,0.0,0.0], [0.0,0.0]],
    upper   = [[50.,50.,50.,50.], [10.,10.]],
)
```

`lower` and `upper` are optional. Omitting them triggers an auto-bounds warning
and uses `guess/10` to `guess*10`.

---

## `MODEL_REGISTRY`

Global `Dict{String, AbstractGrowthModel}` of all built-in models.

```julia
sort(collect(keys(MODEL_REGISTRY)))   # list all names
MODEL_REGISTRY["aHPM"].param_names   # inspect parameters
```

---

## Result types

**`GrowthFitResults`** — returned by `kinbiont_fit`:
```julia
results.data     # GrowthData: the preprocessed data that was fitted
results.results  # Vector{CurveFitResult}: one per curve
```

**`CurveFitResult`** — one per fitted curve:
```julia
r.label          # String: curve identifier
r.best_model     # AbstractGrowthModel: selected by AICc
r.best_params    # Vector{Any}: fitted parameters
r.best_aic       # Float64: AICc of best model
r.fitted_curve   # Vector{Float64}: model values at each time point
r.times          # Vector{Float64}: time points
r.loss           # Float64: final objective value
r.all_results    # Vector{NamedTuple}: raw result for every candidate model
```

`GrowthFitResults` supports indexing and iteration:
```julia
results[1]           # CurveFitResult for the first curve
for r in results     # iterate over all curves
    println(r.label)
end
length(results)      # number of curves fitted
```

---

## CSV format for multi-dimensional data

For `fit_ODEs_System`, `RN_fit`, and `fit_Cybernetic_models`, the input is an
`(m+1) × n_timepoints` matrix where row 1 is time and rows 2..m+1 are the
quantities being fitted:

```
0.0   2.0   4.0   6.0  …
0.09  0.11  0.14  0.18 …   ← species 1
0.10  0.20  0.30  0.40 …   ← species 2
```

---

## ML input format

`downstream_decision_tree_regression` and `downstream_symbolic_regression` take:

1. **Kinbiont_results matrix** — columns = curves, rows = `[label_row, well, model, param_1, …, loss]`
2. **feature_matrix** — rows = curves, columns = `[well_id, feature_1, feature_2, …]`, first row = header

See [ML Downstream](@ref ml) for a full working example built from `save_results` output.
