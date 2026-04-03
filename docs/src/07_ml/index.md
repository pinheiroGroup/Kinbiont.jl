# [ML Downstream](@id ml)

```@contents
Pages = ["index.md"]
Depth = 3
```

Kinbiont's downstream ML connects fitted kinetic parameters to experimental
conditions. The workflow is:

```
kinbiont_fit(...)
    → save_results(...) → _summary.csv
        → load CSV → build Kinbiont_results matrix
            → downstream_decision_tree_regression / downstream_symbolic_regression
```

`save_results` is the bridge: it writes `_summary.csv` with one row per curve,
one column per parameter, ready to feed into the ML functions.

To run all examples:

```julia
using Kinbiont, CSV, DataFrames, DecisionTree, SymbolicRegression
using Plots, Random, AbstractTrees, MLJDecisionTreeInterface
```

---

## Decision tree regression — real data

This example links *E. coli* BW25113 growth parameters to medium composition,
using batch-fit results from 200 chemical media conditions (Reding-Roman et al.,
*Nature Microbiology* 2024, doi:10.1038/s41564-024-01626-9).

```julia
using Kinbiont, CSV, DataFrames, DecisionTree, Random

# ── 1. Load pre-fitted aHPM parameters ───────────────────────────────────────
fit_df = CSV.read("data_examples/batch_fit_sample.csv", DataFrame)
filter!(r -> r.converged, fit_df)

# ── 2. Load medium composition (5 compounds matched to fit_df) ────────────────
feat_df = CSV.read("data_examples/medium_composition_sample.csv", DataFrame)

# ── 3. Join on curve label ────────────────────────────────────────────────────
joined = innerjoin(fit_df, feat_df; on=:label => :Label)

# ── 4. Build Kinbiont_results matrix (rows = [label, well, model, params, loss])
n           = nrow(joined)
param_cols  = [:gr, :exit_lag_rate, :N_max, :shape]
Kinbiont_results = Matrix{Any}(vcat(
    permutedims(fill("batch_fit", n)),     # row 1: experiment label (unused)
    permutedims(string.(joined.label)),    # row 2: well ID
    permutedims(fill("aHPM", n)),          # row 3: model name
    Matrix{Any}(joined[:, param_cols])',   # rows 4–7: gr, exit_lag_rate, N_max, shape
    permutedims(joined.loss),              # row 8: loss
))

# ── 5. Build feature_matrix (first row = header, remaining = data) ─────────────
feat_cols    = [:Glucose_mM, :K2HPO4_mM, :KH2PO4_mM, :Na2HPO4_mM, :NH4Cl_mM]
feat_header  = permutedims(vcat(["label"], string.(feat_cols)))
feat_data    = hcat(string.(joined.label), Matrix{Any}(joined[:, feat_cols]))
feature_matrix = Matrix{Any}(vcat(feat_header, feat_data))

# ── 6. Decision tree: predict growth rate (row 4) from medium compounds ────────
seed = Random.seed!(1234)
dt   = downstream_decision_tree_regression(
    Kinbiont_results,
    feature_matrix,
    4;                       # row 4 = gr
    do_cross_validation = true,
    n_folds_cv          = 5,
    verbose             = true,
)

# dt[1] = fitted tree
# dt[2] = cross-validation R²
# dt[3] = feature importances
println("Cross-validation R² = ", round(dt[2]; digits=3))
println("Feature importances:")
for (feat, imp) in zip(string.(feat_cols), dt[3])
    println("  ", rpad(feat, 14), round(imp; digits=3))
end
```

---

## Decision tree regression — synthetic data

Simulate 100 growth curves under 8 antibiotic combination conditions. Fit each
curve, then use a decision tree to recover the concentration table.

```julia
using Kinbiont, Plots, Random, DecisionTree
Random.seed!(42)

# The "unknown" mapping from antibiotic combination to growth rate scaling
function abx_effect(abx::Vector)
    map = Dict(
        (1,0,0)=>1.0, (0,1,0)=>0.5, (0,0,1)=>0.3,
        (1,1,0)=>0.0, (1,0,1)=>0.3, (0,1,1)=>0.0,
        (1,1,1)=>0.0, (0,0,0)=>1.0,
    )
    return map[Tuple(abx)]
end

base_gr = 0.05
n_exp   = 100
abx_mat = rand(0:1, n_exp, 3)   # random antibiotic combinations

spec = ModelSpec(
    [MODEL_REGISTRY["aHPM"]],
    [[1.0, 1.0, 1.0, 1.0]];
    lower=[[0.0,0.0,0.0,0.0]], upper=[[10.,10.,10.,10.]])
opts = FitOptions(smooth=true)

all_results = []
for i in 1:n_exp
    gr_scaled = base_gr * abx_effect(abx_mat[i, :])
    sim  = ODE_sim("aHPM", [0.05], 0.0, 800.0, 10.0,
                   [gr_scaled, 1.0, 0.5, 2.0])
    times = Float64.(sim.t)
    od    = Float64.(reduce(hcat, sim.u)[1, :]) .+ 0.005 .* randn(length(times))
    data  = GrowthData(reshape(od, 1, length(od)), times, [string(i)])
    res   = kinbiont_fit(data, spec, opts)
    push!(all_results, res[1])
end

# Build Kinbiont_results matrix
n = length(all_results)
Kinbiont_results = Matrix{Any}(vcat(
    permutedims(fill("synth", n)),
    permutedims([r.label for r in all_results]),
    permutedims(fill("aHPM", n)),
    hcat([r.best_params for r in all_results]...),
    permutedims([r.loss for r in all_results]),
))

# Feature matrix
labels_vec  = string.(1:n_exp)
feat_header = ["label" "abx_1" "abx_2" "abx_3"]
feat_data   = hcat(labels_vec, abx_mat)
feature_matrix = Matrix{Any}(vcat(feat_header, feat_data))

seed = Random.seed!(1234)
dt = downstream_decision_tree_regression(
    Kinbiont_results, feature_matrix, 4;   # row 4 = gr
    do_cross_validation=true, n_folds_cv=10, verbose=true)

# Visualise the tree
wt = DecisionTree.wrap(dt[1], (featurenames=["abx_1","abx_2","abx_3"],))
Plots.plot(wt, 0.9, 0.2; size=(1400,700), connect_labels=["yes","no"])
```

---

## Symbolic regression — synthetic data

Discover the relationship between a continuous experimental feature and the
fitted growth rate.

```julia
using Kinbiont, SymbolicRegression, Plots, Random
Random.seed!(42)

unknown_response(x) = 1 / (1 + x)

feature_range = 0.0:0.4:4.0
base_params   = [0.1, 1.0, 50.0, 1.0]

spec = ModelSpec(
    [MODEL_REGISTRY["aHPM"]],
    [[1.0, 1.0, 1.0, 1.0]];
    lower=[[0.0,0.0,0.0,0.0]], upper=[[10.,10.,100.,10.]])
opts = FitOptions(smooth=true)

all_results = []
for f in feature_range
    p = copy(base_params); p[1] *= unknown_response(f)
    sim   = ODE_sim("aHPM", [0.05], 0.0, 800.0, 5.0, p)
    times = Float64.(sim.t)
    od    = Float64.(reduce(hcat, sim.u)[1, :]) .+ 0.01 .* randn(length(times))
    data  = GrowthData(reshape(od,1,length(od)), times, [string(f)])
    push!(all_results, kinbiont_fit(data, spec, opts)[1])
end

n = length(all_results)
Kinbiont_results = Matrix{Any}(vcat(
    permutedims(fill("synth", n)),
    permutedims([r.label for r in all_results]),
    permutedims(fill("aHPM", n)),
    hcat([r.best_params for r in all_results]...),
    permutedims([r.loss for r in all_results]),
))

feature_matrix = permutedims(hcat([string(f), f] for f in feature_range)...)

options = SymbolicRegression.Options(
    binary_operators=[+, /, *, -],
    unary_operators=[],
    parsimony=0.05,
    maxsize=10,
)

gr_sy_reg = downstream_symbolic_regression(Kinbiont_results, feature_matrix, 4; options=options)

scatter(collect(feature_range), [Float64(r.best_params[1]) for r in all_results];
        xlabel="Feature value", ylabel="Growth rate", label="Fitted gr")
plot!(collect(feature_range), unique(gr_sy_reg[3][:, 2]); label="Best equation", linewidth=2)
```
