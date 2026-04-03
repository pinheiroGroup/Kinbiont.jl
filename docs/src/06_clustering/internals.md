# [Clustering in depth](@id clustering_internals)

```@contents
Pages = ["internals.md"]
Depth = 3
```

This page explains how Kinbiont's clustering pipeline works step by step.
For a practical how-to, see [Clustering](@ref clustering).

---

## Algorithm overview

The clustering pipeline runs inside `preprocess()` in this order:

1. **Z-score each curve** row-wise → scale-independent shape representation
2. **Identify constant wells** (optional) via quantile-ratio pre-screening or OLS trend test
3. **Run k-means** on the dynamic subset (or all curves if no pre-screening)
4. **Exponential prototype reassignment** (optional) — post-hoc step
5. **Compute centroids** by averaging the z-scored curves in each final cluster

Steps 2 and 4 are optional and controlled by `FitOptions` fields.

---

## Step 1 — Row-wise z-scoring

Each curve is z-scored independently across its own time points before k-means.
This removes OD scale differences so that two wells with different carrying
capacities but the same sigmoidal shape end up at the same point in shape space.

If a curve's standard deviation is below `1e-12` (effectively constant), the
entire row is replaced by zeros. This means perfectly flat wells collapse to a
single point at the origin, consistent with shape-based comparison.

```julia
using Kinbiont, Plots, Random
Random.seed!(1)

times = collect(0.0:2.0:48.0)
n_tp  = length(times)

# Two logistic curves at different OD scales — same shape
high_od = 0.05 .+ 1.5 ./ (1 .+ exp.(-0.35 .* (times .- 20))) .+ 0.005 .* randn(n_tp)
low_od  = 0.05 .+ 0.3 ./ (1 .+ exp.(-0.35 .* (times .- 20))) .+ 0.001 .* randn(n_tp)
flat    = fill(0.06, n_tp) .+ 0.001 .* randn(n_tp)

data = GrowthData(
    vcat(high_od', low_od', flat'),
    times,
    ["high_OD", "low_OD", "flat"],
)

proc = preprocess(data, FitOptions(cluster=true, n_clusters=2, kmeans_seed=1))

# high_OD and low_OD should cluster together despite different scales
println("Cluster assignments: ", proc.clusters)
# → high_OD and low_OD get the same cluster label

p = plot(xlabel="Time (h)", ylabel="OD (z-scored)", title="Z-scored curves")
for i in 1:3
    # Reconstruct z-scored version for display
    row = data.curves[i, :]
    s   = std(row)
    z   = s < 1e-12 ? zeros(length(row)) : (row .- mean(row)) ./ s
    plot!(p, times, z; label=data.labels[i], linewidth=2)
end
display(p)
```

---

## Step 2 — Identifying non-growing wells

Two mutually exclusive strategies. **Pre-screening takes priority** if both
`cluster_prescreen_constant` and `cluster_trend_test` are `true`.

### Quantile-ratio pre-screening

For each curve the code computes a low quantile `q_low` (default 5th percentile)
and a high quantile `q_high` (default 95th percentile).

- If `q_low > 0`: curve is constant when `q_high ≤ cluster_tol_const × q_low`
- If `q_low ≤ 0`: a fallback rule based on the absolute quantile range is used

Pre-screened constant curves are pinned to label `n_clusters` **before** k-means
runs. K-means then sees only the dynamic subset.

```julia
times = collect(0.0:2.0:48.0)
n_tp  = length(times)
Random.seed!(2)

# Mix: 6 flat + 4 logistic
curves = vcat(
    [fill(0.07, n_tp)' .+ 0.001 .* randn(n_tp)' for _ in 1:6]...,
    [(0.04 .+ 0.8 ./ (1 .+ exp.(-0.35 .* (times .- 18))))' .+ 0.01 .* randn(n_tp)' for _ in 1:4]...,
)
data = GrowthData(curves, times, [string(i) for i in 1:10])

proc = preprocess(data, FitOptions(
    cluster                    = true,
    n_clusters                 = 3,          # constant wells → label 3
    cluster_prescreen_constant = true,
    cluster_tol_const          = 1.9,
    cluster_q_low              = 0.10,
    cluster_q_high             = 0.90,
    kmeans_seed                = 42,
))

println("Assignments (curves 1–6 are flat): ", proc.clusters)
# → curves 1–6 should all have label 3
```

### OLS slope trend test

After plain k-means, each curve's OLS slope is computed analytically and tested
with a two-tailed Student t-test (n − 2 degrees of freedom). If `p ≥ 0.05`,
the curve is reassigned to the flat cluster label `n_clusters`.

This mode is better suited to datasets where all wells are expected to grow but
some have very weak or noisy signal that escapes a quantile filter.

```julia
proc = preprocess(data, FitOptions(
    cluster            = true,
    n_clusters         = 3,
    cluster_trend_test = true,   # default; prescreen takes priority if both set
    kmeans_seed        = 42,
))
println("Trend-test assignments: ", proc.clusters)
```

---

## Step 3 — K-means

K-means is run `kmeans_n_init` times (default 1) with a deterministic
`MersenneTwister`. When `kmeans_seed = 0` (default), seed 42 is used. The same
RNG object is advanced across all repetitions, so the full sequence of starts is
reproducible for a given seed.

The run with the lowest total within-cluster cost (WCSS) is kept.

```julia
# More restarts → more stable partition, higher cost
proc = preprocess(data, FitOptions(
    cluster       = true,
    n_clusters    = 3,
    kmeans_seed   = 123,
    kmeans_n_init = 20,     # default is 1
))
println("WCSS: ", proc.wcss)
```

!!! note "WCSS is from k-means, not from final labels"
    `proc.wcss` always reflects the k-means cost, computed **before** any
    exponential prototype reassignment. It is the right quantity to use for
    elbow plots.

---

## Step 4 — Exponential prototype reassignment (optional)

When `cluster_exp_prototype=true`, Kinbiont builds a bank of z-scored exponential
curves with bases 2⁶ through 2¹⁶, normalised to unit length. Each curve is then
compared to:

1. Its currently assigned k-means centroid (squared Euclidean distance in z-space)
2. The nearest exponential prototype (same metric)

If the exponential prototype is closer, the curve is reassigned to the
exponential cluster.

**How the exponential label is chosen:**

- Preferred label = `n_clusters - 1` when `cluster_prescreen_constant=true`
- Preferred label = `n_clusters` otherwise
- If that label is already occupied by a k-means cluster, the code uses
  `max(current_labels) + 1` instead

This means final labels can exceed `n_clusters` in edge cases.

```julia
times = collect(0.0:2.0:48.0)
n_tp  = length(times)
Random.seed!(3)

function exp_curve(t, rate, amp, base)
    tnorm = (t .- minimum(t)) ./ (maximum(t) - minimum(t))
    base .+ amp .* (exp.(rate .* tnorm) .- 1.0) ./ (exp(rate) - 1.0)
end

logistic_curve(t) = 0.04 .+ 0.8 ./ (1 .+ exp.(-0.35 .* (t .- 20)))

curves = vcat(
    [logistic_curve(times)' .+ 0.01 .* randn(n_tp)' for _ in 1:5]...,
    [exp_curve(times, 5.0, 0.9, 0.04)' .+ 0.01 .* randn(n_tp)' for _ in 1:5]...,
)
data = GrowthData(curves, times, [string(i) for i in 1:10])

proc = preprocess(data, FitOptions(
    cluster               = true,
    n_clusters            = 2,
    cluster_exp_prototype = true,
    kmeans_seed           = 10,
))

println("Assignments (1–5 logistic, 6–10 exponential): ", proc.clusters)
# → curves 6–10 should be reassigned to the exponential prototype cluster
```

---

## Step 5 — Centroid computation

Centroids are computed by averaging the **z-scored** curves in each cluster.
They represent shape prototypes, not average OD trajectories. An empty cluster
produces a zero centroid row.

```julia
# Access centroids after clustering
proc = preprocess(data, FitOptions(cluster=true, n_clusters=3, kmeans_seed=1))

p = plot(xlabel="Time index", ylabel="Z-score",
         title="Cluster centroids (z-normalised space)")
for k in 1:size(proc.centroids, 1)
    n = sum(proc.clusters .== k)
    plot!(p, proc.times, proc.centroids[k, :];
          label="Cluster $k (n=$n)", linewidth=2)
end
display(p)

# If you need original-scale centroids, average the raw curves yourself:
for k in 1:maximum(proc.clusters)
    idx = findall(proc.clusters .== k)
    if !isempty(idx)
        mean_od = vec(mean(data.curves[idx, :], dims=1))
        println("Original-scale centroid for cluster $k: ", round.(mean_od[1:3]; digits=3), " …")
    end
end
```

---

## Edge case: empty curve matrix

If `GrowthData` is constructed with zero curves, `preprocess` returns a
`GrowthData` with empty `clusters` and `centroids`, and `wcss = 0.0`. Downstream
code can safely check `length(proc.clusters) == 0` to detect this.
