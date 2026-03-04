# =============================================================================
# Example 3 — Batch Fitting with Preprocessing and Clustering
# =============================================================================
# Demonstrates the full pipeline:
#   1. Load growth curves from a CSV file (plate-reader format)
#   2. Preprocess: blank subtraction → negative correction → smoothing → clustering
#   3. Fit all curves with model selection
#   4. Inspect cluster assignments and per-cluster parameter distributions
#
# The preprocessing step uses:
#   - StatsBase.zscore (not a custom implementation)
#   - Clustering.kmeans (not a custom k-means)
#   - HypothesisTests.MannKendallTest for flat-curve detection
# =============================================================================

using Kinbiont
using OrdinaryDiffEq
using Distributions
using Statistics
using Random

Random.seed!(99)

# ---------------------------------------------------------------------------
# Simulate a small "plate" of 12 curves: 3 biological conditions × 4 replicates
# Each condition has distinct growth kinetics
# ---------------------------------------------------------------------------

tstart, tmax, delta_t = 0.0, 60.0, 1.5
times = collect(tstart:delta_t:tmax)
n_tp  = length(times)

# Condition parameters: [gr, N_max] for logistic
conditions = [
    ("fast_grower",  [0.30, 1.5]),
    ("slow_grower",  [0.10, 1.2]),
    ("no_growth",    [0.005, 0.1]),  # effectively flat — will be flagged by Mann-Kendall
]

n_replicates = 4
n_curves     = length(conditions) * n_replicates
curves       = Matrix{Float64}(undef, n_curves, n_tp)
labels       = String[]
blank_val    = 0.05   # simulated blank OD

i = 1
for (cond_name, params) in conditions
    for rep in 1:n_replicates
        sim   = ODE_sim("logistic", [0.05], tstart, tmax, delta_t, params)
        noise = rand(Uniform(-0.02, 0.02), n_tp)
        # Add blank offset to simulate raw plate reader data
        curves[i, :] = reduce(hcat, sim.u)[1, :] .+ noise .+ blank_val
        push!(labels, "$(cond_name)_rep$(rep)")
        i += 1
    end
end

raw_data = GrowthData(curves, times, labels)

# ---------------------------------------------------------------------------
# Preprocessing options
# ---------------------------------------------------------------------------

opts = FitOptions(
    # Blank subtraction
    blank_subtraction   = true,
    blank_value         = blank_val,

    # Negative value handling after blank subtraction
    correct_negatives   = true,
    negative_method     = :thr_correction,
    negative_threshold  = 1e-4,

    # Smoothing
    smooth              = true,
    smooth_method       = :lowess,
    lowess_frac         = 0.1,

    # Clustering: k=3, use Mann-Kendall to isolate flat curves
    cluster             = true,
    n_clusters          = 3,
    cluster_trend_test  = true,

    # Fitting
    loss                = "RE",
)

# ---------------------------------------------------------------------------
# Preprocess only (without fitting) — useful for QC inspection
# ---------------------------------------------------------------------------

preprocessed = preprocess(raw_data, opts)

println("=== Preprocessing results ===")
println("  Cluster assignments:")
for (label, cluster_id) in zip(preprocessed.labels, preprocessed.clusters)
    println("    $(rpad(label, 22)) → cluster $(cluster_id)")
end

# Per-cluster curve count
cluster_ids = sort(unique(preprocessed.clusters))
println("\n  Cluster sizes:")
for cid in cluster_ids
    n = count(==(cid), preprocessed.clusters)
    tag = cid > opts.n_clusters ? "  (flat curves — Mann-Kendall p≥0.05)" : ""
    println("    cluster $cid: $n curves$tag")
end

# ---------------------------------------------------------------------------
# Fit all curves (preprocessing is applied internally by fit())
# For efficiency we can also pass the already-preprocessed data
# ---------------------------------------------------------------------------

spec = ModelSpec(
    [MODEL_REGISTRY["NL_logistic"], MODEL_REGISTRY["NL_Gompertz"], MODEL_REGISTRY["NL_exponential"]],
    [
        [1.2, 0.2, 0.05],
        [1.2, 0.2, 20.0],
        [0.05, 0.2],
    ],
)

# Pass preprocessed data directly — fit() skips preprocessing when
# smooth/cluster/etc. are false (default), so we set them all off here
# and use the already-preprocessed data we computed above.
fit_only_opts = FitOptions(loss = "RE")   # no preprocessing — data is already clean

results = fit(preprocessed, spec, fit_only_opts)

# ---------------------------------------------------------------------------
# Results summary
# ---------------------------------------------------------------------------

println("\n=== Fitting results ===")
println(rpad("Label", 26), rpad("Cluster", 10), rpad("Model", 22), rpad("AICc", 10), "Params")

for (r, cluster_id) in zip(results, preprocessed.clusters)
    params_str = join(round.(Float64.(r.best_params); digits=3), ", ")
    println(
        rpad(r.label, 26),
        rpad(string(cluster_id), 10),
        rpad(r.best_model.name, 22),
        rpad(round(r.best_aic; digits=1), 10),
        params_str,
    )
end

# ---------------------------------------------------------------------------
# Per-cluster parameter statistics for the winning model
# ---------------------------------------------------------------------------

println("\n=== Parameter statistics per cluster ===")
for cid in cluster_ids
    cid > opts.n_clusters && continue   # skip the flat-curve cluster
    cluster_results = [r for (r, c) in zip(results, preprocessed.clusters) if c == cid]
    isempty(cluster_results) && continue

    println("\n  Cluster $cid (n=$(length(cluster_results)) curves):")
    # Show mean ± std of first parameter (gr proxy) across replicates
    first_params = [Float64(r.best_params[1]) for r in cluster_results]
    println(
        "    param[1] mean ± std: ",
        round(mean(first_params); digits=4), " ± ", round(std(first_params); digits=4)
    )
    println("    best models: ", join(unique([r.best_model.name for r in cluster_results]), ", "))
end
