# [Clustering](@id clustering)

```@contents
Pages = ["index.md"]
Depth = 3
```

Kinbiont clusters growth curves by **shape** using k-means on z-scored curves.
Z-scoring makes clustering scale-independent — a fast-growing well that reaches
OD 2.0 and a slow-growing well that reaches OD 0.4 can still be grouped by
whether their dynamics are sigmoidal, exponential, or flat.

Clustering is performed inside `preprocess()` **before** blank subtraction, so
that non-growing (blank) wells are still distinguishable by their raw signal and
correctly identified as constant.

---

## Basic usage

Clustering-only mode — no `ModelSpec` needed:

```julia
using Kinbiont

data = GrowthData("data_examples/ecoli_sample.csv")

opts    = FitOptions(cluster=true, n_clusters=4)
results = kinbiont_fit(data, opts)

# Assignments and shape prototypes are stored in results.data
println(results.data.clusters)   # Vector{Int}: one label per curve, in 1..4
println(results.data.wcss)       # Float64: within-cluster sum of squares
# results.data.centroids          # 4 × n_timepoints Matrix{Float64} (z-scored)
```

Cluster labels are always integers in `1..n_clusters` (or up to `n_clusters+1`
when `cluster_exp_prototype=true` allocates an extra slot).

---

## Choosing k — elbow plot

Run clustering for a range of k values and plot the within-cluster sum of squares
(WCSS). The "elbow" — where WCSS stops falling steeply — suggests the optimal k.

Data: E. coli Keio knockout strains, *Nature Scientific Data* 2026,
doi:10.1038/s41597-026-07075-9.

```julia
using Kinbiont, Plots

data = GrowthData("data_examples/ecoli_sample.csv")

ks   = 2:8
wcss = [preprocess(data,
            FitOptions(cluster=true, n_clusters=k,
                       cluster_prescreen_constant=true)).wcss
        for k in ks]

scatter(ks, wcss;
        xlabel="Number of clusters (k)", ylabel="WCSS",
        title="Elbow plot — E. coli knockout growth curves",
        legend=false, linewidth=2, marker=:circle)
```

!!! tip "Reading the elbow"
    Pick the k where adding another cluster gives diminishing WCSS reduction.
    If the curve is monotone with no clear bend, try `cluster_prescreen_constant=true`
    to handle non-growing wells separately.

---

## Handling non-growing wells

Two strategies for identifying flat / non-growing curves:

### Strategy 1: `cluster_trend_test` (default)

After k-means, an OLS slope t-test flags curves with no significant linear trend
(p ≥ 0.05) and re-labels them to `n_clusters`. K-means runs on all `n_clusters - 1`
dynamic groups first.

```julia
opts = FitOptions(
    cluster            = true,
    n_clusters         = 4,
    cluster_trend_test = true,   # default — no need to set explicitly
)
```

### Strategy 2: `cluster_prescreen_constant` (recommended for noisy plates)

Pre-screens before k-means using a quantile-ratio criterion. Curves where
`q_high / q_low ≤ cluster_tol_const` are pinned to label `n_clusters`; k-means
then runs only on dynamic curves. More biologically meaningful because k-means
never sees flat wells.

```julia
opts = FitOptions(
    cluster                    = true,
    n_clusters                 = 4,
    cluster_prescreen_constant = true,
    cluster_tol_const          = 1.5,   # increase → fewer "constant" calls
    cluster_q_low              = 0.05,
    cluster_q_high             = 0.95,
)
```

!!! note "Which to use?"
    `cluster_prescreen_constant=true` is generally better for microplate data
    where blank wells and non-growing knockouts are common. Use `cluster_trend_test`
    for datasets where all wells are expected to grow.

---

## Exponential prototype cluster

Add a dedicated cluster for curves that look exponential rather than sigmoidal.
Kinbiont builds z-scored exponential prototypes (bases 2⁶..2¹⁶) and reassigns
each curve to this cluster if it is closer to an exponential prototype than to
its k-means centroid.

```julia
opts = FitOptions(
    cluster                    = true,
    n_clusters                 = 4,
    cluster_prescreen_constant = true,
    cluster_exp_prototype      = true,
    # Label allocation:
    # n_clusters     → constant wells (from prescreen)
    # n_clusters - 1 → exponential prototype cluster
    # 1..n_clusters-2 → k-means dynamic groups
)
```

When the preferred exponential slot is already occupied by k-means,
the exponential prototype is assigned label `max(labels) + 1`.

---

## Reproducibility

K-means uses a fixed random seed (`42`) by default. To set your own:

```julia
opts = FitOptions(
    cluster        = true,
    n_clusters     = 4,
    kmeans_seed    = 1234,   # any non-zero Int → deterministic
    kmeans_n_init  = 20,     # run k-means 20 times, keep lowest WCSS
)
```

`kmeans_seed = 0` (the default) also uses seed 42 (not the global RNG).

---

## Inspecting and subsetting clusters

```julia
data = GrowthData("data_examples/ecoli_sample.csv")
opts = FitOptions(cluster=true, n_clusters=4, cluster_prescreen_constant=true)
res  = kinbiont_fit(data, opts)
d    = res.data   # preprocessed GrowthData with cluster assignments

# How many curves per cluster?
for k in 1:4
    n = sum(d.clusters .== k)
    println("Cluster $k: $n curves")
end

# Plot centroids (z-scored shape prototypes)
using Plots
p = plot(xlabel="Time (h)", ylabel="OD (z-scored)",
         title="Cluster centroids")
for k in 1:size(d.centroids, 1)
    n = sum(d.clusters .== k)
    plot!(p, d.times, d.centroids[k, :]; label="Cluster $k (n=$n)", linewidth=2)
end
display(p)

# Subset to curves in cluster 2 only
idx      = d.clusters .== 2
cluster2 = d[d.labels[idx]]   # returns a new GrowthData
```

---

## Two-step workflow: cluster then fit

Cluster first to identify growing wells; fit only those.

```julia
using Kinbiont

data = GrowthData("data_examples/ecoli_sample.csv")

# Step 1: Cluster to identify non-growing wells
cluster_opts = FitOptions(
    cluster                    = true,
    n_clusters                 = 4,
    cluster_prescreen_constant = true,
)
clustered = kinbiont_fit(data, cluster_opts)
d = clustered.data

# Step 2: Keep only growing wells (all clusters except n_clusters = constant)
growing_labels = d.labels[d.clusters .!= 4]
growing        = d[growing_labels]

# Step 3: Fit growing wells
spec = ModelSpec(
    [MODEL_REGISTRY["aHPM"]],
    [[1.0, 1.0, 1.0, 1.0]];
    lower = [[0.0, 0.0, 0.0, 0.0]],
    upper = [[50.0, 50.0, 50.0, 50.0]],
)
fit_opts = FitOptions(smooth=true, smooth_method=:rolling_avg)
results  = kinbiont_fit(growing, spec, fit_opts)

save_results(results, "output/growing_only/")
```
