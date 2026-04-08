# Clustering analysis of growth curves

This section provides examples of the clustering module used in the unified preprocessing pipeline.
The purpose of this module is to group growth curves according to their **shape**, rather than their absolute OD magnitude.
This is useful for distinguishing different classes of temporal behaviour among the growing curves.

The clustering procedure is based on **k-means applied to row-wise z-scored curves**, with optional extensions for:

- pre-screening constant curves before k-means,
- post-hoc reassignment of flat curves based on a linear trend test,
- forced shift-detection with post-hoc reassignment of curves that are closer to exponential prototypes than to their assigned k-means centroid.

```@contents
Pages = ["index.md"]
Depth = 3
```

## Clustering of growth curves

In this section, we present different examples of how to use the clustering part of the preprocessing pipeline.
To run these examples, you will typically need the following Julia packages:

```julia
using Clustering
using StatsBase
using Distributions
using Random
using Statistics
using Plots
```

The clustering logic is implemented internally through the following functions:

```julia
_kmeans_best(X, k, opts)
_cluster(curves, times, opts)
_prescreen_constant(curves, opts)
_cluster_with_prescreen(curves, zscored_all, opts)
_compute_centroids(curves, labels, n_clusters)
_build_exp_prototypes(times)
_apply_exp_prototype_labels(zscored, labels, exp_protos, centroids_norm, exp_label)
_apply_trend_labels(curves, times, labels, flat_id)
_zscore_rows(curves)
_exp_prototype_label(opts, labels)
```

The main internal entry point is:

```julia
_cluster(curves, times, opts)
```

This function returns:

```julia
(labels, centroids, wcss)
```

where:

- `labels` is the vector of cluster assignments,
- `centroids` is the matrix of cluster centroids in **z-normalised space**,
- `wcss` is the within-cluster sum of squares returned by the **k-means step**.

When exponential prototype relabeling is enabled, `wcss` is **not recomputed** after relabeling: 
it still refers to the k-means solution obtained before the exponential reassignment step.

## Basic clustering on z-scored curves

In this first example, we simulate a small set of synthetic growth curves with different temporal shapes 
and cluster them using plain k-means on z-scored trajectories.

The purpose of this example is to show the simplest use of the clustering module, 
without constant pre-screening, trend testing, or exponential prototype relabeling.

First, we define the time axis:

```julia
times = collect(0.0:2.0:48.0)
n_tp = length(times)
```

Next, we define some helper functions to generate synthetic curves with different behaviours:

```julia
function logistic_like_curve(t, lag, rate, amp, baseline)
    return baseline .+ amp ./ (1 .+ exp.(-rate .* (t .- lag)))
end

function flat_curve(t, level)
    return fill(level, length(t))
end

function delayed_growth_curve(t, lag, rate, amp, baseline)
    y = baseline .+ amp ./ (1 .+ exp.(-rate .* (t .- lag)))
    y[t .< lag] .= baseline
    return y
end
```

We now generate a synthetic dataset composed of flat curves, standard logistic-like curves, and delayed-growth curves:

```julia
Random.seed!(1234)

curves = Matrix{Float64}(undef, 12, n_tp)

for i in 1:4
    curves[i, :] = flat_curve(times, 0.08) .+ randn(n_tp) .* 0.002
end

for i in 5:8
    curves[i, :] = logistic_like_curve(times, 20.0, 0.35, 0.8, 0.05) .+ randn(n_tp) .* 0.01
end

for i in 9:12
    curves[i, :] = delayed_growth_curve(times, 28.0, 0.45, 0.9, 0.05) .+ randn(n_tp) .* 0.01
end
```

We can visualise the simulated curves:

```julia
p = plot(xlabel="Time", ylabel="OD", size=(500, 350), label=nothing)
for i in 1:size(curves, 1)
    plot!(p, times, curves[i, :], label=nothing)
end
display(p)
```

We now define the clustering options and run the clustering:

```julia
opts = FitOptions(
    cluster=true,
    n_clusters=3,
    cluster_prescreen_constant=false,
    cluster_trend_test=false,
    cluster_exp_prototype=false,
    kmeans_n_init=10,
    kmeans_seed=1234,
    kmeans_max_iters=1000,
    kmeans_tol=1e-6
)

labels, centroids, wcss = _cluster(curves, times, opts)
```

The clustering result can be inspected through:

```julia
labels
centroids
wcss
```

We can also plot the curves coloured by their assigned cluster:

```julia
for k in 1:maximum(labels)
    idx = findall(labels .== k)
    p = plot(title="Cluster $k", xlabel="Time", ylabel="OD", size=(400, 300))
    for i in idx
        plot!(p, times, curves[i, :], label=nothing)
    end
    display(p)
end
```

In this example, the flat curves and the two kinds of growing curves are expected 
to be separated mainly according to their **shape**, because the clustering is 
performed after row-wise z-score normalisation.

## Behaviour of row-wise z-score normalisation

The helper used internally is:

```julia
_zscore_rows(curves)
```

Each curve is z-scored independently across its own time points.
If a curve is effectively constant, the whole row is replaced by zeros:

```julia
if s < 1e-12
    out[i, :] .= 0.0
else
    out[i, :] = zscore(row)
end
```

This means that perfectly flat curves become all-zero trajectories in z-normalised space, 
which is consistent with the intended shape-based comparison.

## Constant pre-screening of flat curves

In many datasets, some wells are clearly non-growing or nearly constant.
In this case, it can be useful to identify these curves before running 
k-means, so that the dynamic clusters are learned only from the genuinely varying trajectories.

This behaviour is activated with:

```julia
cluster_prescreen_constant=true
```

The constant detection is based on a quantile-ratio criterion implemented by:

```julia
_prescreen_constant(curves, opts)
```

For each curve, the code computes:

- a low quantile `q_low`,
- a high quantile `q_high`.

If the low quantile is positive, a curve is classified as constant when:

```julia
q_high <= cluster_tol_const * q_low
```

If the low quantile is non-positive, the code switches to a fallback rule based on whether the overall quantile range is negligible.

### Constant pre-screening: synthetic example

In this example, we simulate a dataset with many flat curves and a smaller set of dynamic curves.

```julia
Random.seed!(2222)

curves = Matrix{Float64}(undef, 15, n_tp)

for i in 1:7
    curves[i, :] = flat_curve(times, 0.06) .+ randn(n_tp) .* 0.001
end

for i in 8:11
    curves[i, :] = logistic_like_curve(times, 18.0, 0.30, 0.75, 0.04) .+ randn(n_tp) .* 0.01
end

for i in 12:15
    curves[i, :] = delayed_growth_curve(times, 30.0, 0.50, 0.95, 0.04) .+ randn(n_tp) .* 0.01
end
```

We define clustering options with constant pre-screening enabled:

```julia
opts = FitOptions(
    cluster=true,
    n_clusters=3,
    cluster_prescreen_constant=true,
    cluster_q_low=0.10,
    cluster_q_high=0.90,
    cluster_tol_const=1.9,
    cluster_trend_test=false,
    cluster_exp_prototype=false,
    kmeans_n_init=10,
    kmeans_seed=42,
    kmeans_max_iters=1000,
    kmeans_tol=1e-6
)

labels, centroids, wcss = _cluster(curves, times, opts)
```

In this mode:

- constant curves are assigned directly to cluster `n_clusters`,
- k-means is run only on the dynamic subset,
- the dynamic curves receive labels `1:(n_clusters-1)`.

The corresponding labels can be inspected through:

```julia
labels
```

To visualise the result:

```julia
for k in 1:maximum(labels)
    idx = findall(labels .== k)
    p = plot(title="Cluster $k", xlabel="Time", ylabel="OD", size=(400, 300))
    for i in idx
        plot!(p, times, curves[i, :], label=nothing)
    end
    display(p)
end
```

A special edge case is also handled explicitly in the code: 
if all curves are classified as constant, no k-means step is run and `wcss` remains `0.0`.

## Trend-test reassignment of flat curves

A different strategy for handling flat curves is to first run 
k-means and then reassign curves whose overall slope is not statistically significant.

This behaviour is activated with:

```julia
cluster_trend_test=true
```

In this case, the code:

1. runs k-means on all z-scored curves using `max(1, n_clusters - 1)` clusters,
2. then applies a post-hoc trend test on the **original curves**,
3. reassigns curves with non-significant slope to a dedicated flat cluster with label `n_clusters`.

The helper used for this step is:

```julia
_apply_trend_labels(curves, times, labels, flat_id)
```

This function computes the OLS slope analytically and evaluates a two-tailed 
t-test on that slope using a Student t distribution with `n - 2` degrees of freedom.

The reassignment rule is:

- if `p >= 0.05`, the curve is reassigned to `flat_id`,
- otherwise, it keeps its k-means label.

### Trend-test reassignment: synthetic example

In this example, we simulate curves that are not perfectly constant, but still have very weak overall trend.

```julia
Random.seed!(3333)

curves = Matrix{Float64}(undef, 15, n_tp)

for i in 1:5
    curves[i, :] = 0.08 .+ 0.002 .* sin.(times ./ 3) .+ randn(n_tp) .* 0.002
end

for i in 6:10
    curves[i, :] = logistic_like_curve(times, 16.0, 0.35, 0.7, 0.05) .+ randn(n_tp) .* 0.01
end

for i in 11:15
    curves[i, :] = delayed_growth_curve(times, 30.0, 0.55, 0.9, 0.05) .+ randn(n_tp) .* 0.01
end
```

We now activate the trend-test mode:

```julia
opts = FitOptions(
    cluster=true,
    n_clusters=3,
    cluster_prescreen_constant=false,
    cluster_trend_test=true,
    cluster_exp_prototype=false,
    kmeans_n_init=10,
    kmeans_seed=7,
    kmeans_max_iters=1000,
    kmeans_tol=1e-6
)

labels, centroids, wcss = _cluster(curves, times, opts)
```

The cluster labels can be inspected through:

```julia
labels
```

And visualised with:

```julia
for k in 1:maximum(labels)
    idx = findall(labels .== k)
    p = plot(title="Cluster $k", xlabel="Time", ylabel="OD", size=(400, 300))
    for i in idx
        plot!(p, times, curves[i, :], label=nothing)
    end
    display(p)
end
```

This strategy is useful when flat curves are noisy enough to escape the 
quantile-ratio filter, but still do not show a statistically significant linear trend.

## Priority between pre-screening and trend-test modes

The code does **not** combine constant pre-screening and trend-test reassignment as two sequential steps.

The internal priority is:

1. if `cluster_prescreen_constant=true`, use pre-screening mode;
2. else if `cluster_trend_test=true`, use trend-test mode;
3. otherwise, use plain k-means.

After this choice, exponential prototype relabeling may still be applied as an additional final step.

So, if both `cluster_prescreen_constant=true` and `cluster_trend_test=true`, the trend-test branch is **ignored**.

## Exponential prototype relabeling

The clustering module also supports a post-processing step that identifies curves 
closer to a family of idealised exponential shapes than to their assigned k-means centroid.

This behaviour is activated with:

```julia
cluster_exp_prototype=true
```

In this case, the code:

1. builds a set of z-scored exponential prototypes,
2. computes the squared distance between each curve and the nearest exponential prototype,
3. computes the squared distance to the curve's currently assigned centroid,
4. reassigns the curve if the exponential prototype is closer.

The helper functions used are:

```julia
_build_exp_prototypes(times)
_apply_exp_prototype_labels(zscored, labels, exp_protos, centroids_norm, exp_label)
_exp_prototype_label(opts, labels)
```

### How the exponential label is chosen

The exponential cluster label is not arbitrary.

The code first computes a preferred label:

- `n_clusters - 1` when `cluster_prescreen_constant=true`,
- `n_clusters` otherwise.

If that preferred label is already present in the current labels, the code allocates a new label equal to `maximum(labels) + 1`.

So when exponential prototype relabeling is enabled, the final labels may extend beyond `1:opts.n_clusters`.

### Exponential prototype relabeling: synthetic example

We simulate a dataset containing flat curves, logistic-like curves, and strongly exponential trajectories.

```julia
function exponential_like_curve(t, rate, amp, baseline)
    tnorm = (t .- minimum(t)) ./ (maximum(t) - minimum(t))
    return baseline .+ amp .* ((exp.(rate .* tnorm) .- 1.0) ./ (exp(rate) - 1.0))
end
```

```julia
Random.seed!(4444)

curves = Matrix{Float64}(undef, 15, n_tp)

for i in 1:5
    curves[i, :] = flat_curve(times, 0.05) .+ randn(n_tp) .* 0.002
end

for i in 6:10
    curves[i, :] = logistic_like_curve(times, 20.0, 0.30, 0.8, 0.04) .+ randn(n_tp) .* 0.01
end

for i in 11:15
    curves[i, :] = exponential_like_curve(times, 5.0, 0.9, 0.04) .+ randn(n_tp) .* 0.01
end
```

We now activate exponential prototype relabeling:

```julia
opts = FitOptions(
    cluster=true,
    n_clusters=3,
    cluster_prescreen_constant=false,
    cluster_trend_test=false,
    cluster_exp_prototype=true,
    kmeans_n_init=10,
    kmeans_seed=100,
    kmeans_max_iters=1000,
    kmeans_tol=1e-6
)

labels, centroids, wcss = _cluster(curves, times, opts)
```

The output can be inspected through:

```julia
labels
centroids
wcss
```

To visualise the result:

```julia
for k in 1:maximum(labels)
    idx = findall(labels .== k)
    p = plot(title="Cluster $k", xlabel="Time", ylabel="OD", size=(400, 300))
    for i in idx
        plot!(p, times, curves[i, :], label=nothing)
    end
    display(p)
end
```

## Centroid computation

After the final cluster labels are determined, the centroids are computed by averaging the z-scored curves assigned to each cluster.

This is performed in the actual workflow by:

```julia
centroids = _compute_centroids(zscored_all, labels, n_effective)
```

The helper itself is generic:

```julia
_compute_centroids(curves, labels, n_clusters)
```

It simply averages the rows of the matrix that it receives.
Therefore, `_compute_centroids` does **not** perform z-scoring on its own.
The centroids are in z-normalised space only because, in the real pipeline, the function is called on `zscored_all`.

The returned centroid matrix is therefore expressed in **z-normalised space**, not in the original OD scale.
This means that the centroids represent **shape prototypes**, not average absolute OD trajectories.

If a cluster is empty, the corresponding centroid row remains a vector of zeros.

### Centroid computation: simple example

```julia
zscored = _zscore_rows(curves)
centroids = _compute_centroids(zscored, labels, maximum(labels))
```

The centroids can be plotted directly:

```julia
p = plot(title="Cluster centroids (z-normalised space)", xlabel="Time index", ylabel="Z-score", size=(500, 350))
for k in 1:size(centroids, 1)
    plot!(p, 1:size(centroids, 2), centroids[k, :], label="Cluster $k")
end
display(p)
```

If needed, original-scale prototypes can be reconstructed manually by averaging the original curves within each cluster:

```julia
for k in 1:maximum(labels)
    idx = findall(labels .== k)
    if !isempty(idx)
        original_centroid = vec(mean(curves[idx, :], dims=1))
        println("Original-scale centroid for cluster $k:")
        println(original_centroid)
    end
end
```

## Repeated k-means initialisations

The k-means routine is wrapped by:

```julia
_kmeans_best(X, k, opts)
```

This function runs k-means `opts.kmeans_n_init` times and returns the result with the lowest total cost.

The implementation uses a deterministic `MersenneTwister`, with:

- seed `42` when `opts.kmeans_seed == 0`,
- the user-provided seed otherwise.

The same RNG object is then reused across repeated initialisations, so the sequence of starts is deterministic for a given seed.

### Example of repeated initialisation

```julia
opts = FitOptions(
    cluster=true,
    n_clusters=4,
    kmeans_n_init=20,
    kmeans_seed=123,
    kmeans_max_iters=5000,
    kmeans_tol=1e-8
)

labels, centroids, wcss = _cluster(curves, times, opts)
```

This is useful when the user wants a more stable partition and wants to reduce the chance of keeping a poor local optimum.

## Edge case: empty input

The clustering code explicitly handles the case where the input contains zero curves.

If `size(curves, 1) == 0`, the function returns:

- an empty label vector,
- a zero centroid matrix of size `n_clusters × n_timepoints`,
- `wcss = 0.0`.

This behaviour ensures that downstream code can still rely on a consistent return type.

## Full clustering workflow example

In this final example, we simulate a mixed dataset and compare the behaviour of different clustering modes.

We first generate the dataset:

```julia
Random.seed!(5555)

curves = Matrix{Float64}(undef, 18, n_tp)

for i in 1:6
    curves[i, :] = flat_curve(times, 0.06) .+ randn(n_tp) .* 0.002
end

for i in 7:12
    curves[i, :] = logistic_like_curve(times, 18.0, 0.35, 0.8, 0.04) .+ randn(n_tp) .* 0.01
end

for i in 13:18
    curves[i, :] = exponential_like_curve(times, 4.5, 0.9, 0.04) .+ randn(n_tp) .* 0.01
end
```

We first run plain clustering:

```julia
opts_plain = FitOptions(
    cluster=true,
    n_clusters=3,
    cluster_prescreen_constant=false,
    cluster_trend_test=false,
    cluster_exp_prototype=false,
    kmeans_n_init=10,
    kmeans_seed=11
)

labels_plain, centroids_plain, wcss_plain = _cluster(curves, times, opts_plain)
```

Then clustering with constant pre-screening:

```julia
opts_const = FitOptions(
    cluster=true,
    n_clusters=3,
    cluster_prescreen_constant=true,
    cluster_q_low=0.10,
    cluster_q_high=0.90,
    cluster_tol_const=1.9,
    cluster_trend_test=false,
    cluster_exp_prototype=false,
    kmeans_n_init=10,
    kmeans_seed=11
)

labels_const, centroids_const, wcss_const = _cluster(curves, times, opts_const)
```

Finally clustering with exponential prototype relabeling:

```julia
opts_exp = FitOptions(
    cluster=true,
    n_clusters=3,
    cluster_prescreen_constant=false,
    cluster_trend_test=false,
    cluster_exp_prototype=true,
    kmeans_n_init=10,
    kmeans_seed=11
)

labels_exp, centroids_exp, wcss_exp = _cluster(curves, times, opts_exp)
```

The results can be compared through:

```julia
labels_plain
wcss_plain

labels_const
wcss_const

labels_exp
wcss_exp
```

This workflow illustrates how the choice of clustering mode affects the final partition of the growth curves.

## Summary

The clustering module groups curves according to **shape**, not raw OD scale, by applying k-means to row-wise z-scored trajectories.

Depending on the chosen options, the user can:

- run plain k-means on all curves,
- reserve one cluster for clearly constant curves using quantile-ratio pre-screening,
- reserve one cluster for flat curves using a slope significance test,
- add an exponential-like cluster by comparing curves against analytical prototypes.

These options make the clustering stage flexible and suitable 
for both exploratory analyses and more structured downstream workflows in growth-curve analysis.
