# IrregularGrowthData — Design Spec
_Date: 2026-04-08_

## Background

`GrowthData` assumes all curves share one common time vector. Two real scenarios break this:

1. Curves are sampled at genuinely different time points (different plate readers, different start times, missing measurements).
2. `negative_method = :remove` drops leading/trailing points per curve, leaving each curve with a different valid range.

A collaborator proposed normalizing each curve's time axis to [0, 1] (removing absolute scale), building a union grid in [0, 1], resampling all curves onto it, then clustering. This spec defines how to integrate that idea into KinBiont as a first-class type.

---

## Scope

- New type `IrregularGrowthData` in a new file `src/api/irregular.jl`.
- New `preprocess` method dispatching on `IrregularGrowthData` (clustering only).
- New test file `test/api/test_irregular.jl`.
- New example + validation script `examples/04_irregular_time_clustering.jl` that runs both the standalone approach (from `newClusteringIdea.jl`) and the KinBiont approach on the same data and compares cluster assignments.
- No changes to `GrowthData`, `FitOptions`, `preprocess(::GrowthData, ...)`, or any fitting code.

---

## Architecture

### New file: `src/api/irregular.jl`

#### Private helpers

```julia
_normalize_times_01(times::Vector{Float64}) -> Vector{Float64}
```
Maps `times` to [0, 1] via `(t .- tmin) ./ (tmax - tmin)`. Returns `zeros` when span ≤ 0.

```julia
_build_union_grid(times01_list::Vector{Vector{Float64}}; step::Float64=0.01) -> Vector{Float64}
```
Snaps every normalized time point to the nearest `step` multiple, collects unique values, ensures 0.0 and 1.0 are included. Returns a sorted vector.

```julia
_interp_linear(x::Vector{Float64}, y::Vector{Float64}, x_new::Vector{Float64}) -> Vector{Float64}
```
Piecewise linear interpolation. `x` must be sorted. Points outside `[x[1], x[end]]` are clamped to the boundary values. Uses `searchsortedfirst` for correctness (fixes the walking-pointer limitation in the collaborator's version).

```julia
_resample_to_union_grid(
    times01_list::Vector{Vector{Float64}},
    values_list::Vector{Vector{Float64}},
    union_grid::Vector{Float64},
) -> Matrix{Float64}
```
For each curve: snap its normalized times to `step`, deduplicate collisions by keeping the first, sort, then call `_interp_linear` onto `union_grid`. Returns an `n_curves × length(union_grid)` matrix.

#### Public struct

```julia
struct IrregularGrowthData
    raw_curves  :: Vector{Vector{Float64}}   # original values per curve
    raw_times   :: Vector{Vector{Float64}}   # original (unnormalized) times per curve
    labels      :: Vector{String}
    curves      :: Matrix{Float64}           # n_curves × n_grid, resampled on [0,1]
    times       :: Vector{Float64}           # normalized [0,1] union grid
    clusters    :: Union{Nothing, Vector{Int}}
    centroids   :: Union{Nothing, Matrix{Float64}}
    wcss        :: Union{Nothing, Float64}
end
```

Inner constructor signature:
```julia
IrregularGrowthData(
    raw_curves  :: Vector{Vector{Float64}},
    raw_times   :: Vector{Vector{Float64}},
    labels      :: Vector{String};
    step        :: Float64 = 0.01,
    clusters    = nothing,
    centroids   = nothing,
    wcss        = nothing,
)
```

Validation checks (all throw `ArgumentError`):
- `length(raw_curves) == length(raw_times) == length(labels)`
- Every `raw_times[i]` has the same length as `raw_curves[i]`
- Every `raw_times[i]` has length ≥ 2

Constructor steps:
1. Normalize each `raw_times[i]` → `times01_list`
2. Build union grid from `times01_list` with `step`
3. Resample all curves onto union grid → `curves` matrix
4. Store all fields

#### Subsetting

```julia
Base.getindex(data::IrregularGrowthData, labels::Vector{String}) -> IrregularGrowthData
```
Returns a new `IrregularGrowthData` with only the named curves. Clustering fields are dropped (same as `GrowthData`). Throws `ArgumentError` on unknown labels. Re-uses the already-computed `curves` and `times` (no re-resampling needed since the union grid is shared).

---

### Edit: `src/api/preprocessing.jl`

New method:
```julia
function preprocess(data::IrregularGrowthData, opts::FitOptions)::IrregularGrowthData
```

Behavior:
- Issues `@warn` if any of `blank_subtraction`, `correct_negatives`, or `smooth` are `true`, telling the user to apply those on raw data before constructing `IrregularGrowthData`.
- Runs `_cluster(data.curves, data.times, opts)` when `opts.cluster = true`.
- Returns a new `IrregularGrowthData` with the same `raw_curves`, `raw_times`, `labels`, `curves`, `times` but updated `clusters`, `centroids`, `wcss`. The inner constructor is bypassed on the return path (use the keyword-arg form directly, no re-resampling).
- `average_replicates` and `scattering_correction` are also warned and skipped.

---

### Edit: `src/Kinbiont.jl`

Add `include("api/irregular.jl")` after `include("api/types.jl")`.

`irregular.jl` ends with `export IrregularGrowthData`.

---

### New file: `test/api/test_irregular.jl`

Test groups:

1. **Constructor validation** — length mismatches, single-point curves throw `ArgumentError`; valid input stores raw fields correctly.
2. **`_normalize_times_01`** — output in [0,1], handles zero-span.
3. **`_build_union_grid`** — always contains 0.0 and 1.0; step respected; sorted.
4. **`_interp_linear`** — exact at known points; clamped outside range; linear between points.
5. **Resampling round-trip** — construct `IrregularGrowthData` from 3 curves with known shapes; verify `data.curves[i, :]` matches expected values at grid points.
6. **Clustering via `preprocess`** — clusters populated, length matches n_curves, labels in 1..n_clusters, WCSS ≥ 0.
7. **`getindex` subsetting** — correct curves selected; clusters dropped; unknown label throws.
8. **Unsupported opts warning** — `@test_warn` that blank_subtraction, correct_negatives, smooth each emit a warning.
9. **Pure function** — `preprocess` does not mutate input `raw_curves`.

---

### Edit: `test/runtests.jl`

Add `include("api/test_irregular.jl")`.

---

### New file: `examples/04_irregular_time_clustering.jl`

Standalone script with two sections:

**Section A — Standalone approach** (inline from `newClusteringIdea.jl`, self-contained):
- Generates `n_curves = 240` synthetic curves with `seed = 2026` using the collaborator's shape generators.
- Runs normalization → union grid → resample → zscore → kmeans (`k=6`, 15 inits, seed=2026).
- Stores `assignments_standalone`.

**Section B — KinBiont approach**:
- Uses the same generated `times_list`, `values_list`, `true_labels`.
- Constructs `IrregularGrowthData(values_list, times_list, true_labels; step=0.01)`.
- Calls `preprocess(data, FitOptions(cluster=true, n_clusters=6, kmeans_seed=2026, kmeans_n_init=15, cluster_trend_test=false))`.
- Stores `assignments_kinbiont`.

**Section C — Comparison**:
- Prints a contingency table (`n_clusters × n_clusters` matrix where `M[i,j]` = number of curves assigned to standalone cluster `i` and KinBiont cluster `j`).
- Prints the fraction of curves where assignments agree after resolving the label permutation via the Hungarian algorithm (`HungarianAlgorithm.jl`) or by finding the best-matching permutation for small `k`.
- Prints per-shape-class assignment counts to show that both methods recover the true labels equally well.

Note: perfect agreement is expected when both use identical seeds and the same `Clustering.kmeans` backend.

---

## What changes from `newClusteringIdea.jl`

| Original | Fate in KinBiont |
|---|---|
| `logistic_curve`, `exponential_curve`, `make_shape_curve` | Moved to example file only |
| `generate_irregular_times`, `generate_dataset` | Moved to example file only |
| `normalize_times_01` | → `_normalize_times_01` in `irregular.jl` |
| `build_union_grid` | → `_build_union_grid` in `irregular.jl` |
| `interpolate_on_grid` | → `_interp_linear` in `irregular.jl`; uses `searchsortedfirst` instead of walking pointer |
| `resample_all_curves` | Inlined into `IrregularGrowthData` constructor via `_resample_to_union_grid` |
| `zscore_rows` | Deleted — `_zscore_rows` already in `preprocessing.jl` |
| `run_kmeans_best` | Deleted — `_kmeans_best` already in `preprocessing.jl` |
| `plot_all_clusters_z`, `plot_all_clusters_original` | Moved to example file only |
| Main script block | Becomes Section A of `examples/04_irregular_time_clustering.jl` |

---

## Files changed / created

| File | Change |
|---|---|
| `src/api/irregular.jl` | **New** |
| `src/Kinbiont.jl` | Add `include` + export |
| `src/api/preprocessing.jl` | Add `preprocess(::IrregularGrowthData, ...)` |
| `test/api/test_irregular.jl` | **New** |
| `test/runtests.jl` | Add `include` |
| `examples/04_irregular_time_clustering.jl` | **New** |
