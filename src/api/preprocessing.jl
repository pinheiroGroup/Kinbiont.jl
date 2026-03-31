# =============================================================================
# Unified Preprocessing Pipeline
# =============================================================================
# preprocess(data::GrowthData, opts::FitOptions) -> GrowthData
#
# Pure function — always returns a new GrowthData, never mutates the input.
# Delegates all actual computation to the existing functions already in
# pre_processing_functions.jl, and uses:
#   - StatsBase.zscore   for z-score normalisation (replaces custom zscore_vector)
#   - Clustering.kmeans  for k-means clustering (replaces custom implementation)
#   - linear slope t-test (Statistics stdlib) for flat-curve trend detection
# =============================================================================

using Clustering: kmeans, assignments
using StatsBase: zscore
using Distributions: TDist, cdf
using Random: MersenneTwister

# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

"""
    preprocess(data::GrowthData, opts::FitOptions) -> GrowthData

Apply the preprocessing pipeline defined by `opts` to `data`.
Returns a **new** `GrowthData`; the input is never modified.

Steps applied in order (each is a no-op when the corresponding option is `false`):
1. Multiple-scattering OD correction (`opts.scattering_correction`)
2. K-means clustering on z-scored curves (`opts.cluster`) — intentionally
   performed **before** blank subtraction so that blank/non-growing wells are
   still distinguishable by their raw signal.
3. Blank subtraction (`opts.blank_subtraction`)
4. Negative-value correction (`opts.correct_negatives`)
5. Smoothing (`opts.smooth`)

# Example
```julia
opts = FitOptions(smooth=true, cluster=true, n_clusters=4)
processed = preprocess(raw_data, opts)
# processed.clusters now holds a cluster id for every curve
```
"""
function preprocess(data::GrowthData, opts::FitOptions)::GrowthData
    curves = data.curves   # n_curves × n_tp, never mutated in place
    times  = data.times
    labels = data.labels

    # Replicate averaging collapses curves with identical labels before any other step
    curves, labels = _apply_replicate_averaging(curves, labels, opts)

    curves = _apply_scattering_correction(curves, times, opts)

    # Clustering is performed on scattering-corrected but otherwise raw data so
    # that blank subtraction does not hide the distinction between growing and
    # non-growing wells (which is precisely what clustering is meant to find).
    clusters, centroids, wcss = opts.cluster ? _cluster(curves, times, opts) : (nothing, nothing, nothing)

    curves = _apply_blank_subtraction(curves, data.labels, opts)
    curves = _apply_negative_correction(curves, times, opts)
    curves, times = _apply_smoothing(curves, times, opts)   # Gaussian may change times

    return GrowthData(curves, times, labels, clusters, centroids, wcss)
end

# ---------------------------------------------------------------------------
# Step 0 — Replicate averaging
# ---------------------------------------------------------------------------

"""
    _apply_replicate_averaging(curves, labels, opts) -> (Matrix{Float64}, Vector{String})

When `opts.average_replicates=true`, average all curves that share the same label
into a single row. Wells labelled `"b"` (blank) or `"X"` (discard) are excluded
from averaging and dropped from the output — they are not biologically meaningful
replicates. Returns the merged curves matrix and the deduplicated label vector.

When `opts.average_replicates=false`, returns the inputs unchanged.
"""
function _apply_replicate_averaging(
    curves::Matrix{Float64},
    labels::Vector{String},
    opts::FitOptions,
)::Tuple{Matrix{Float64}, Vector{String}}
    opts.average_replicates || return curves, labels

    skip = ("b", "X")
    unique_labels = unique(l for l in labels if l ∉ skip)
    isempty(unique_labels) && return curves[Int[], :], String[]

    n_tp  = size(curves, 2)
    out   = Matrix{Float64}(undef, length(unique_labels), n_tp)
    for (i, lbl) in enumerate(unique_labels)
        idx = findall(==(lbl), labels)
        out[i, :] = vec(mean(curves[idx, :], dims=1))
    end
    return out, unique_labels
end

# ---------------------------------------------------------------------------
# Step 1 — Multiple-scattering OD correction
# ---------------------------------------------------------------------------

function _apply_scattering_correction(
    curves::Matrix{Float64},
    times::Vector{Float64},
    opts::FitOptions,
)::Matrix{Float64}
    opts.scattering_correction || return curves
    isempty(opts.calibration_file) && error(
        "scattering_correction=true requires a calibration_file path in FitOptions"
    )

    method = string(opts.scattering_method)
    corrected = similar(curves)
    for i in axes(curves, 1)
        # correction_OD_multiple_scattering expects a 2×n_tp matrix [times; values]
        curve_mat = Matrix(transpose(hcat(times, curves[i, :])))
        result = correction_OD_multiple_scattering(curve_mat, opts.calibration_file; method)
        corrected[i, :] = result[2, :]
    end
    return corrected
end

# ---------------------------------------------------------------------------
# Step 2 — Blank subtraction
# ---------------------------------------------------------------------------

function _apply_blank_subtraction(
    curves::Matrix{Float64},
    labels::Vector{String},
    opts::FitOptions,
)::Matrix{Float64}
    opts.blank_subtraction || return curves

    if opts.blank_from_labels
        blank_val = _blank_from_labels(curves, labels)
    else
        blank_val = opts.blank_value
    end

    return curves .- blank_val
end

"""
    _blank_from_labels(curves, labels) -> Float64

Compute the mean OD across all timepoints and all wells whose label is `"b"`.
Warns and returns 0.0 when no blank wells are found.
"""
function _blank_from_labels(curves::Matrix{Float64}, labels::Vector{String})::Float64
    blank_idx = findall(==("b"), labels)
    if isempty(blank_idx)
        @warn "blank_from_labels=true but no wells labelled \"b\" found; blank value set to 0.0"
        return 0.0
    end
    return mean(curves[blank_idx, :])
end

# ---------------------------------------------------------------------------
# Step 3 — Negative-value correction
# ---------------------------------------------------------------------------

function _apply_negative_correction(
    curves::Matrix{Float64},
    times::Vector{Float64},
    opts::FitOptions,
)::Matrix{Float64}
    opts.correct_negatives || return curves

    method = string(opts.negative_method)
    corrected = similar(curves)
    for i in axes(curves, 1)
        curve_mat = Matrix(transpose(hcat(times, curves[i, :])))
        result = negative_value_correction(
            curve_mat,
            Float64[];
            method,
            thr_negative = opts.negative_threshold,
        )
        # negative_value_correction may shorten the curve when method=:remove
        if size(result, 2) == length(times)
            corrected[i, :] = result[2, :]
        else
            @warn "Curve $(i): negative removal changed length; curve kept as-is"
            corrected[i, :] = curves[i, :]
        end
    end
    return corrected
end

# ---------------------------------------------------------------------------
# Step 4 — Smoothing
# Returns (smoothed_curves, times); times may change when :gaussian + time_grid
# ---------------------------------------------------------------------------

function _apply_smoothing(
    curves::Matrix{Float64},
    times::Vector{Float64},
    opts::FitOptions,
)::Tuple{Matrix{Float64}, Vector{Float64}}
    opts.smooth || return curves, times
    opts.smooth_method == :none && return curves, times

    n_curves = size(curves, 1)
    n_curves == 0 && return curves, times

    opts.smooth_method == :boxcar && return _apply_boxcar_smoothing(curves, times, opts)

    if opts.smooth_method == :boxcar
        return _apply_boxcar_smoothing(curves, times, opts)
    end

    smoothing_str = _smoothing_symbol_to_string(opts.smooth_method)

    # Process first curve to determine output length (smoothing may shorten data,
    # and :gaussian with a custom time grid may produce a different number of points)
    curve_mat1 = Matrix(transpose(hcat(times, curves[1, :])))
    result1 = smoothing_data(
        curve_mat1;
        method             = smoothing_str,
        pt_avg             = opts.smooth_pt_avg,
        thr_lowess         = opts.lowess_frac,
        gaussian_h_mult    = opts.gaussian_h_mult,
        gaussian_time_grid = opts.gaussian_time_grid,
    )
    n_out     = size(result1, 2)
    new_times = Vector{Float64}(result1[1, :])
    smoothed  = Matrix{Float64}(undef, n_curves, n_out)
    smoothed[1, :] = result1[2, :]

    for i in 2:n_curves
        curve_mat = Matrix(transpose(hcat(times, curves[i, :])))
        result = smoothing_data(
            curve_mat;
            method             = smoothing_str,
            pt_avg             = opts.smooth_pt_avg,
            thr_lowess         = opts.lowess_frac,
            gaussian_h_mult    = opts.gaussian_h_mult,
            gaussian_time_grid = opts.gaussian_time_grid,
        )
        ni = size(result, 2)
        if ni == n_out
            smoothed[i, :] = result[2, :]
        else
            # Unexpected length difference — fall back to first n_out values or padding
            m = min(ni, n_out)
            smoothed[i, 1:m]       = result[2, 1:m]
            smoothed[i, m+1:n_out] .= result[2, end]
        end
    end
    return smoothed, new_times
end

_smoothing_symbol_to_string(s::Symbol) = Dict(
    :lowess      => "lowess",
    :rolling_avg => "rolling_avg",
    :gaussian    => "gaussian",
    :none        => "NO",
)[s]

# ---------------------------------------------------------------------------
# Boxcar (symmetric moving average) smoothing
# Keeps the original time grid; each point is replaced by the mean of the
# symmetric window [j-half, j+half], clamped at the array boundaries.
# ---------------------------------------------------------------------------

"""
    _apply_boxcar_smoothing(curves, times, opts) -> (Matrix{Float64}, Vector{Float64})

Apply a centered boxcar (symmetric moving-average) filter to every curve.
Unlike `:rolling_avg`, the time grid is unchanged and no points are dropped.
Window width is `opts.boxcar_window` (must be ≥ 1).
"""
function _apply_boxcar_smoothing(
    curves::Matrix{Float64},
    times::Vector{Float64},
    opts::FitOptions,
)::Tuple{Matrix{Float64}, Vector{Float64}}
    w = max(1, opts.boxcar_window)
    half = w ÷ 2
    n_curves, n_tp = size(curves)
    smoothed = Matrix{Float64}(undef, n_curves, n_tp)

    for i in 1:n_curves
        curve = @view curves[i, :]
        for j in 1:n_tp
            lo = max(1, j - half)
            hi = min(n_tp, j + half)
            smoothed[i, j] = mean(@view curve[lo:hi])
        end
    end

    return smoothed, times   # time grid is preserved
end

# ---------------------------------------------------------------------------
# Step 5 — K-means clustering on z-scored curves
# ---------------------------------------------------------------------------

# Run k-means `opts.kmeans_n_init` times and return the result with the lowest WCSS.
# Uses a seeded MersenneTwister when `opts.kmeans_seed != 0` for reproducibility.
function _kmeans_best(X::AbstractMatrix{Float64}, k::Int, opts::FitOptions)
    rng    = opts.kmeans_seed == 0 ? MersenneTwister(42) : MersenneTwister(opts.kmeans_seed)
    n_init = max(1, opts.kmeans_n_init)
    best   = kmeans(X, k; maxiter=opts.kmeans_max_iters, tol=opts.kmeans_tol, rng=rng)
    for _ in 2:n_init
        r = kmeans(X, k; maxiter=opts.kmeans_max_iters, tol=opts.kmeans_tol, rng=rng)
        r.totalcost < best.totalcost && (best = r)
    end
    return best
end

"""
    _cluster(curves, times, opts) -> (labels, centroids, wcss)

Cluster growth curves and return per-curve labels, per-cluster shape centroids,
and the within-cluster sum of squares (WCSS / total SSE from k-means).

Labels are always in `1..opts.n_clusters`. `centroids` is an
`n_clusters × n_timepoints` matrix in **z-normalised space**: each row is the
mean of the z-scored curves assigned to that cluster. This captures the *shape*
of each cluster independently of absolute OD magnitude. To recover original-space
prototypes, compute the mean of `data.curves[data.clusters .== k, :]` for each `k`.
`wcss` is the total SSE from the k-means step (0.0 when all curves were classified
as constant).

Modes (evaluated in priority order):
1. `cluster_prescreen_constant=true`: quantile-ratio pre-screening identifies
   non-growing wells before k-means, which then runs on dynamic curves only.
2. `cluster_trend_test=true`: OLS slope t-test post-hoc re-labels flat curves.
3. Neither: plain k-means on all curves.

When `cluster_exp_prototype=true`, an additional exponential shape cluster is
added after k-means: each curve is tested against pre-built z-scored exponential
prototypes and re-labelled if it is closer to them than to any k-means centroid.
"""
function _cluster(
    curves::Matrix{Float64},
    times::Vector{Float64},
    opts::FitOptions,
)::Tuple{Vector{Int}, Matrix{Float64}, Float64}
    size(curves, 1) == 0 && return Int[], zeros(Float64, opts.n_clusters, size(curves, 2)), 0.0

    # Z-score all curves once; used for both k-means and centroid computation.
    zscored_all = _zscore_rows(curves)

    if opts.cluster_prescreen_constant
        labels, wcss = _cluster_with_prescreen(curves, zscored_all, opts)
    elseif opts.cluster_trend_test
        k_dynamic = max(1, opts.n_clusters - 1)
        result = _kmeans_best(zscored_all', k_dynamic, opts)
        labels = assignments(result)
        wcss   = result.totalcost
        labels = _apply_trend_labels(curves, times, labels, opts.n_clusters)
    else
        result = _kmeans_best(zscored_all', opts.n_clusters, opts)
        labels = assignments(result)
        wcss   = result.totalcost
    end

    if opts.cluster_exp_prototype
        exp_label      = _exp_prototype_label(opts, labels)
        exp_protos     = _build_exp_prototypes(times)
        centroids_norm = _compute_centroids(zscored_all, labels, opts.n_clusters)
        labels         = _apply_exp_prototype_labels(zscored_all, labels, exp_protos,
                                                     centroids_norm, exp_label)
    end

    n_effective = isempty(labels) ? opts.n_clusters : maximum(labels)
    # Centroids in z-normalised space (shape prototypes, scale-independent).
    centroids = _compute_centroids(zscored_all, labels, n_effective)
    return labels, centroids, wcss
end

# The exponential prototype occupies a label that is NOT already in use by k-means.
# Preferred slot: n_clusters-1 (constant pre-screening mode) or n_clusters (plain mode).
# If that slot is already occupied by k-means, allocate max(labels)+1 instead.
function _exp_prototype_label(opts::FitOptions, labels::Vector{Int})
    preferred = opts.cluster_prescreen_constant ? opts.n_clusters - 1 : opts.n_clusters
    preferred ∈ labels ? maximum(labels) + 1 : preferred
end

# ---------------------------------------------------------------------------
# Constant pre-screening helpers
# ---------------------------------------------------------------------------

"""
    _prescreen_constant(curves, opts) -> BitVector

Identify non-growing (constant) curves using a quantile-ratio criterion.
Curve `i` is flagged when its `q_high` quantile ≤ `tol_const × q_low` quantile,
meaning its signal never grows meaningfully above its baseline.
"""
function _prescreen_constant(curves::Matrix{Float64}, opts::FitOptions)::BitVector
    W = size(curves, 1)
    const_mask = falses(W)
    for i in 1:W
        row = curves[i, :]
        lm  = quantile(row, opts.cluster_q_low)
        hm  = quantile(row, opts.cluster_q_high)
        const_mask[i] = if lm > 0
            hm <= opts.cluster_tol_const * lm
        else
            # baseline is at or below zero — flat if the total range is negligible
            abs(hm - lm) < 1e-6 * (abs(lm) + 1.0)
        end
    end
    return const_mask
end

"""
    _cluster_with_prescreen(curves, zscored_all, opts) -> (Vector{Int}, Float64)

Pre-screen constant curves using original-space quantile ratio, then run k-means
on the z-normalised dynamic subset. Constant curves receive label `n_clusters`;
dynamic curves get labels `1..n_clusters-1`.
Returns labels and the WCSS from the k-means step (0.0 when all curves are constant).
"""
function _cluster_with_prescreen(
    curves::Matrix{Float64},
    zscored_all::Matrix{Float64},
    opts::FitOptions,
)::Tuple{Vector{Int}, Float64}
    W           = size(curves, 1)
    const_mask  = _prescreen_constant(curves, opts)   # quantile test on raw curves
    dynamic_idx = findall(.!const_mask)
    labels      = fill(opts.n_clusters, W)   # default: constant
    wcss        = 0.0

    if !isempty(dynamic_idx) && opts.n_clusters > 1
        k_dynamic = opts.n_clusters - 1
        # k-means runs on z-normalised dynamic curves
        result    = _kmeans_best(zscored_all[dynamic_idx, :]', k_dynamic, opts)
        km_labels = assignments(result)
        wcss      = result.totalcost
        for (pos, idx) in enumerate(dynamic_idx)
            labels[idx] = km_labels[pos]
        end
    end

    return labels, wcss
end

# ---------------------------------------------------------------------------
# Centroid computation
# ---------------------------------------------------------------------------

"""
    _compute_centroids(curves, labels, n_clusters) -> Matrix{Float64}

Compute the mean curve (in z-normalised space) for each cluster.
Returns an `n_clusters × n_timepoints` matrix; empty clusters produce a zero row.
"""
function _compute_centroids(
    curves::Matrix{Float64},
    labels::Vector{Int},
    n_clusters::Int,
)::Matrix{Float64}
    n_tp      = size(curves, 2)
    centroids = zeros(Float64, n_clusters, n_tp)
    for k in 1:n_clusters
        idxs = findall(==(k), labels)
        isempty(idxs) && continue
        centroids[k, :] = vec(mean(curves[idxs, :], dims=1))
    end
    return centroids
end

# ---------------------------------------------------------------------------
# Exponential prototype helpers
# ---------------------------------------------------------------------------

"""
    _build_exp_prototypes(times) -> Vector{Vector{Float64}}

Build z-scored exponential shape prototypes using bases 2⁶..2¹⁶.
Each prototype is `(b^τ - 1)/(b - 1)` evaluated at normalised time τ ∈ [0,1],
then z-score normalised so distances are shape-only comparisons.
"""
function _build_exp_prototypes(times::Vector{Float64})::Vector{Vector{Float64}}
    tmin, tmax = minimum(times), maximum(times)
    dt = max(tmax - tmin, 1e-12)
    tau = (times .- tmin) ./ dt

    protos = Vector{Vector{Float64}}()
    for exp_b in (2^6, 2^7, 2^8, 2^9, 2^10, 2^11, 2^12, 2^13, 2^14, 2^15, 2^16)
        bf = Float64(exp_b)
        y  = [(bf^t - 1.0) / (bf - 1.0) for t in tau]
        push!(protos, _zscore_vec(y))
    end
    return protos
end

@inline function _zscore_vec(x::Vector{Float64})::Vector{Float64}
    mu = mean(x); s = std(x; corrected=false)
    s < 1e-12 ? zeros(length(x)) : (x .- mu) ./ s
end

"""
    _apply_exp_prototype_labels(zscored, labels, exp_protos, centroids_norm, exp_label)

For each curve, compare its squared distance to the nearest exponential prototype
against the distance to the nearest k-means centroid. If the exponential is closer,
reassign the curve to `exp_label`.
"""
function _apply_exp_prototype_labels(
    zscored::Matrix{Float64},
    labels::Vector{Int},
    exp_protos::Vector{Vector{Float64}},
    centroids_norm::Matrix{Float64},
    exp_label::Int,
)::Vector{Int}
    new_labels = copy(labels)
    for i in axes(zscored, 1)
        xi = @view zscored[i, :]
        # distance to nearest exponential prototype
        d_exp = minimum(sum((xi .- p) .^ 2) for p in exp_protos)
        # distance to assigned k-means centroid
        lbl   = labels[i]
        d_km  = sum((xi .- @view centroids_norm[lbl, :]) .^ 2)
        if d_exp < d_km
            new_labels[i] = exp_label
        end
    end
    return new_labels
end

# Z-score each row (curve) over its time points; handle constant rows gracefully
function _zscore_rows(curves::Matrix{Float64})::Matrix{Float64}
    out = similar(curves)
    for i in axes(curves, 1)
        row = curves[i, :]
        s = std(row)
        if s < 1e-12
            out[i, :] .= 0.0   # constant curve → all zeros after z-score
        else
            out[i, :] = zscore(row)   # StatsBase.zscore
        end
    end
    return out
end

# Assign a dedicated cluster id to curves with no significant linear trend.
# Uses a two-tailed t-test on the OLS slope (p ≥ 0.05 → flat).
# `flat_id` is passed in by the caller; it must already be within 1..n_clusters.
function _apply_trend_labels(
    curves::Matrix{Float64},
    times::Vector{Float64},
    labels::Vector{Int},
    flat_id::Int,
)::Vector{Int}
    new_labels = copy(labels)
    n          = length(times)
    n < 3      && return new_labels
    t_centered = times .- mean(times)
    ss_t       = sum(t_centered .^ 2)
    ss_t <= 0  && return new_labels

    for i in axes(curves, 1)
        y         = curves[i, :]
        slope     = sum(t_centered .* y) / ss_t
        y_hat     = mean(y) .+ slope .* t_centered
        residuals = y .- y_hat
        s2        = sum(residuals .^ 2) / (n - 2)
        se_slope  = sqrt(s2 / ss_t)
        if se_slope <= 0 || !isfinite(se_slope)
            if abs(slope) < 1e-12
                new_labels[i] = flat_id
            end
            continue
        end
        t_stat   = slope / se_slope
        # Two-tailed p-value from the t-distribution with n-2 degrees of freedom
        p_approx = 2 * (1 - cdf(TDist(n - 2), abs(t_stat)))
        if p_approx >= 0.05
            new_labels[i] = flat_id
        end
    end
    return new_labels
end

# ---------------------------------------------------------------------------
# Stationary phase cutoff detection
# ---------------------------------------------------------------------------

"""
    _find_stationary_cutoff(data_mat, opts) -> Int

Given a `2×n` data matrix `[times; od_values]` and a `FitOptions`, return the
column index at which the stationary phase begins. Returns `size(data_mat, 2)`
(full length) when no stationary phase is detected.

The algorithm:
1. Restrict to time points where OD > `opts.stationary_thr_od`.
2. Compute the specific growth rate (SGR) over the restricted data.
3. Find the first window of `opts.stationary_win_size` consecutive points after
   the SGR maximum where all SGR values are below
   `maximum(SGR) × opts.stationary_percentile_thr`.
4. Snap the cutoff forward to the OD peak within the next
   `opts.stationary_win_size` time points (the OD may still rise after the SGR
   threshold is crossed).
"""
function _find_stationary_cutoff(data_mat::Matrix{Float64}, opts::FitOptions)::Int
    n = size(data_mat, 2)
    index_od = findall(data_mat[2, :] .> opts.stationary_thr_od)
    isempty(index_od) && return n

    data_t = data_mat[:, index_od]
    sgr = specific_gr_evaluation(data_t, opts.stationary_pt_smooth_derivative)
    thr = maximum(sgr) * opts.stationary_percentile_thr
    max_idx = argmax(sgr)
    win = opts.stationary_win_size

    # Guard: if fewer than win_size points exceed the threshold, the curve is
    # flat / not growing — skip detection and use the full dataset (matches the
    # behaviour of the legacy find_stationary_phase).
    length(findall(sgr .> thr)) > win || return n

    for i in max_idx:(length(sgr) - win)
        if all(sgr[i:(i + win - 1)] .< thr)
            base_idx   = i + index_od[1] - 1
            search_end = min(base_idx + win, n)
            peak_offset = argmax(data_mat[2, base_idx:search_end])
            return base_idx + peak_offset - 1
        end
    end

    return n
end

export preprocess
