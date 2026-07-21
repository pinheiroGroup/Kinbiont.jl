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

using Clustering: kmeans, kmedoids, hclust, cutree, dbscan, assignments
using StatsBase: zscore
using Distributions: Normal, TDist, cdf
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

    # Resolve the scalar or point-by-point blank before replicate averaging,
    # because averaging removes "b"-labelled wells.
    blank = _resolve_blank_value(curves, labels, opts)

    curves, labels = _apply_replicate_averaging(curves, labels, opts)

    curves = _apply_scattering_correction(curves, times, opts)

    # Clustering is performed on scattering-corrected but otherwise raw data so
    # that blank subtraction does not hide the distinction between growing and
    # non-growing wells (which is precisely what clustering is meant to find).
    clusters, centroids, wcss = opts.cluster ? _cluster(curves, times, opts) : (nothing, nothing, nothing)

    curves = _apply_blank_subtraction(curves, blank, opts)
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

"""
    apply_blank_timeseries(curves, blank_timeseries; method=:pointbypoint, floor=nothing)

Apply one blank trace to every curve. `:pointbypoint` subtracts the trace,
`:shift` subtracts its finite mean and translates each row uniformly, and
`:clip` subtracts its finite mean and floors the result. Passing a numeric
`floor` also makes `:pointbypoint` translate each row above that floor; the
default `nothing` preserves the legacy subtraction-only behaviour.
"""
function apply_blank_timeseries(
    curves::Matrix{Float64},
    blank_timeseries::Vector{Float64};
    method::Symbol=:pointbypoint,
    floor::Union{Nothing, Float64}=nothing,
)::Matrix{Float64}
    n_tp = size(curves, 2)
    length(blank_timeseries) >= n_tp || throw(ArgumentError(
        "blank_timeseries length $(length(blank_timeseries)) < n_timepoints $n_tp"
    ))
    ts = blank_timeseries[1:n_tp]

    if method in (:pointbypoint, :pointwise, :time_avg)
        ts_safe = [isfinite(v) ? v : 0.0 for v in ts]
        corrected = curves .- reshape(ts_safe, 1, :)
        return floor === nothing ? corrected : _shift_rows_to_floor(corrected, floor)
    elseif method == :shift
        finite = filter(isfinite, ts)
        isempty(finite) && return copy(curves)
        return _shift_rows_to_floor(curves .- mean(finite), something(floor, 0.0))
    elseif method == :clip
        finite = filter(isfinite, ts)
        isempty(finite) && return copy(curves)
        return max.(curves .- mean(finite), something(floor, 0.0))
    end
    throw(ArgumentError("Unknown blank correction method: $method"))
end

function _shift_rows_to_floor(curves::Matrix{Float64}, floor::Float64)::Matrix{Float64}
    shifted = copy(curves)
    for i in axes(shifted, 1)
        finite_values = filter(isfinite, shifted[i, :])
        isempty(finite_values) && continue
        delta = max(floor - minimum(finite_values), 0.0)
        shifted[i, :] .+= delta
    end
    return shifted
end

"""
    _resolve_blank_value(curves, labels, opts) -> Union{Float64,Vector{Float64}}

Determine the scalar or point-by-point blank to subtract. Must be called
**before** replicate averaging, because averaging drops `"b"`-labelled wells.

Returns 0.0 when `opts.blank_subtraction` is false.
"""
function _resolve_blank_value(
    curves::Matrix{Float64},
    labels::Vector{String},
    opts::FitOptions,
)::Union{Float64, Vector{Float64}}
    opts.blank_subtraction || return 0.0

    method = opts.blank_method
    if method in (:pointbypoint, :pointwise, :time_avg)
        if opts.blank_from_labels
            return _blank_timeseries_from_labels(curves, labels)
        end
        opts.blank_timeseries === nothing && throw(ArgumentError(
            "blank_method=:pointbypoint requires blank_timeseries or blank_from_labels=true"
        ))
        length(opts.blank_timeseries) == size(curves, 2) || throw(ArgumentError(
            "blank_timeseries length $(length(opts.blank_timeseries)) does not match " *
            "the number of time points $(size(curves, 2))"
        ))
        return [isfinite(v) ? v : 0.0 for v in opts.blank_timeseries]
    elseif method in (:global, :average, :avg_blank, :avg_subtraction, :shift, :clip)
        return opts.blank_from_labels ? _blank_from_labels(curves, labels) : opts.blank_value
    end
    throw(ArgumentError("Unknown blank_method: $(opts.blank_method)"))
end

function _apply_blank_subtraction(
    curves::Matrix{Float64},
    blank::Union{Float64, Vector{Float64}},
    opts::FitOptions,
)::Matrix{Float64}
    opts.blank_subtraction || return curves
    method = opts.blank_method
    if method in (:pointbypoint, :pointwise, :time_avg)
        return apply_blank_timeseries(curves, blank; method=:pointbypoint, floor=opts.blank_floor)
    elseif method in (:shift, :clip)
        blank_trace = fill(blank, size(curves, 2))
        return apply_blank_timeseries(curves, blank_trace; method, floor=opts.blank_floor)
    end
    return curves .- blank
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

"""
    _blank_timeseries_from_labels(curves, labels) -> Vector{Float64}

Compute the per-timepoint mean across wells labelled `"b"`. Non-finite values
are ignored; a time point with no finite blank measurement falls back to the
mean of the other timepoint means, or 0.0 when no finite blank data exist.
"""
function _blank_timeseries_from_labels(
    curves::Matrix{Float64},
    labels::Vector{String},
)::Vector{Float64}
    blank_idx = findall(==("b"), labels)
    if isempty(blank_idx)
        @warn "blank_from_labels=true but no wells labelled \"b\" found; blank timeseries set to 0.0"
        return zeros(Float64, size(curves, 2))
    end

    blank = Vector{Float64}(undef, size(curves, 2))
    for j in axes(curves, 2)
        values = filter(isfinite, curves[blank_idx, j])
        blank[j] = isempty(values) ? NaN : mean(values)
    end
    finite_means = filter(isfinite, blank)
    fallback = isempty(finite_means) ? 0.0 : mean(finite_means)
    for j in eachindex(blank)
        isfinite(blank[j]) || (blank[j] = fallback)
    end
    return blank
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
Window width is `opts.boxcar_window` (an odd integer of at least 3).
"""
function _apply_boxcar_smoothing(
    curves::Matrix{Float64},
    times::Vector{Float64},
    opts::FitOptions,
)::Tuple{Matrix{Float64}, Vector{Float64}}
    w = opts.boxcar_window
    w >= 3 || throw(ArgumentError("boxcar_window must be greater than or equal to 3"))
    isodd(w) || throw(ArgumentError("boxcar_window must be odd so the average is centered"))
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
# Uses a deterministic MersenneTwister. When `opts.kmeans_seed == 0`, a default
# fixed seed (42) is used; otherwise the user-provided seed is used.
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
    _cluster(curves, times, opts) -> (labels, centroids, cost)

Cluster growth curves and return per-curve labels, per-cluster shape centroids,
and the method-specific clustering cost stored in the historical `wcss` field.

Without exponential prototype relabeling, labels lie in `1..opts.n_clusters`.
When `cluster_exp_prototype=true`, an additional label may be allocated (up to
`opts.n_clusters + 1`) if the preferred slot is already occupied by k-means.
`centroids` is an `n_effective × n_timepoints` matrix in **z-normalised space**:
each row is the mean of the z-scored curves assigned to that cluster. This captures
the *shape* of each cluster independently of absolute OD magnitude. To recover
original-space prototypes, compute the mean of `data.curves[data.clusters .== k, :]`
for each `k`.
The cost is WCSS for k-means and hierarchical clustering, total distance to assigned
medoids for k-medoids, and `0.0` for DBSCAN. Sentinel curves separated by a non-growing
criterion are excluded from that cost.

When the quantile pre-screen and slope trend test are enabled together, their union
forms one non-growing group. For methods with a predefined `k`, clustering runs on
the remaining dynamic curves; DBSCAN applies the same reassignment post-hoc.

When `cluster_exp_prototype=true`, an additional exponential shape cluster is
added after k-means: each curve is tested against pre-built z-scored exponential
prototypes and re-labelled if it is closer to them than to any k-means centroid.
"""
function _pairwise_euclidean(X::Matrix{Float64})::Matrix{Float64}
    n = size(X, 1)
    D = Matrix{Float64}(undef, n, n)
    for i in 1:n, j in 1:n
        D[i, j] = sqrt(sum((X[i, k] - X[j, k])^2 for k in axes(X, 2)))
    end
    return D
end

function _cluster(
    curves::Matrix{Float64},
    times::Vector{Float64},
    opts::FitOptions,
)::Tuple{Vector{Int}, Matrix{Float64}, Float64}
    size(curves, 1) == 0 && return Int[], zeros(Float64, opts.n_clusters, size(curves, 2)), 0.0

    zscored_all = _zscore_rows(curves)
    method      = opts.cluster_method

    # When both detectors are enabled, their union forms one non-growing class.
    # DBSCAN applies the same union post-hoc because it has no k.
    if (opts.cluster_prescreen_constant || opts.cluster_trend_test) && method != :dbscan
        non_growing = _non_growing_mask(curves, times, opts)
        labels, wcss = _cluster_with_non_growing_method(
            zscored_all, non_growing, opts
        )
    else
        labels, wcss = _cluster_dispatch(zscored_all, opts)

        # DBSCAN has no k parameter, so the selected non-growing criteria are
        # applied post-hoc using one fresh label above its cluster ids.
        if method == :dbscan &&
           (opts.cluster_prescreen_constant || opts.cluster_trend_test)
            non_growing = _non_growing_mask(curves, times, opts)
            if any(non_growing)
                flat_id = isempty(labels) ? 1 : maximum(labels) + 1
                labels[non_growing] .= flat_id
            end
        end
    end

    # Exponential prototype (kmeans only)
    if opts.cluster_exp_prototype && method == :kmeans
        exp_label      = _exp_prototype_label(opts, labels)
        exp_protos     = _build_exp_prototypes(times)
        centroids_norm = _compute_centroids(zscored_all, labels, opts.n_clusters)
        labels         = _apply_exp_prototype_labels(zscored_all, labels, exp_protos,
                                                     centroids_norm, exp_label)
    end

    # Report WCSS with a single, method-independent definition: the total SSE
    # between every standardized curve and the mean (centroid) of its assigned
    # cluster, computed over ALL curves — including any flat/constant sentinel
    # cluster and DBSCAN noise. This overrides the per-method cost so the value
    # matches the documented "sum of squared Euclidean distances between every
    # standardized trajectory and its assigned cluster centroid": k-medoids
    # otherwise returns a sum of *unsquared* distances to medoids, and the
    # prescreen/trend paths otherwise report the dynamic-subset cost only,
    # excluding the set-aside curves.
    wcss = _wcss(zscored_all, labels)

    pos_labels  = filter(>(0), labels)
    n_effective = (isempty(labels) || isempty(pos_labels)) ? opts.n_clusters : maximum(pos_labels)
    n_effective = max(n_effective, 1)
    centroids   = _compute_centroids(zscored_all, labels, n_effective)
    return labels, centroids, wcss
end

# Within-cluster sum of squares (total SSE): for every standardized curve, the
# squared Euclidean distance to the mean (centroid) of its assigned cluster,
# summed over all curves. Each distinct label defines a group — including a
# reserved flat/constant sentinel cluster and DBSCAN noise (label 0) — and each
# curve is measured against its own group's mean. Method-independent so the
# elbow/sweep plot uses one consistent definition.
function _wcss(zscored::Matrix{Float64}, labels::Vector{Int})::Float64
    (isempty(labels) || size(zscored, 1) == 0) && return 0.0
    total = 0.0
    for lbl in unique(labels)
        idxs = findall(==(lbl), labels)
        sub  = zscored[idxs, :]
        c    = vec(mean(sub, dims=1))
        for i in axes(sub, 1)
            total += sum((sub[i, :] .- c) .^ 2)
        end
    end
    return total
end

function _non_growing_mask(
    curves::Matrix{Float64},
    times::Vector{Float64},
    opts::FitOptions,
)::BitVector
    mask = falses(size(curves, 1))
    opts.cluster_prescreen_constant && (mask .|= _prescreen_constant(curves, opts))
    opts.cluster_trend_test && (mask .|= _flat_curve_mask(
        curves, times; p_threshold=opts.cluster_trend_p_thr
    ))
    return mask
end

"""
    detect_non_growing_indices(curves, times; kwargs...) -> Vector{Int}

Identify curves selected by the clustering non-growing controls. The
quantile-ratio pre-screen and the slope trend test can be enabled separately;
when both are enabled, a curve is selected when either criterion matches.
"""
function detect_non_growing_indices(
    curves::Matrix{Float64},
    times::Vector{Float64};
    prescreen_constant::Bool=false,
    trend_test::Bool=false,
    prescreen_tol::Float64=1.5,
    prescreen_q_low::Float64=0.05,
    prescreen_q_high::Float64=0.95,
    trend_p_threshold::Float64=0.05,
)::Vector{Int}
    opts = FitOptions(
        cluster_prescreen_constant=prescreen_constant,
        cluster_trend_test=trend_test,
        cluster_tol_const=prescreen_tol,
        cluster_q_low=prescreen_q_low,
        cluster_q_high=prescreen_q_high,
        cluster_trend_p_thr=trend_p_threshold,
    )
    return findall(_non_growing_mask(curves, times, opts))
end

# Dispatch to the appropriate clustering algorithm on z-scored curves.
function _cluster_dispatch(
    zscored::Matrix{Float64},
    opts::FitOptions,
)::Tuple{Vector{Int}, Float64}
    method = opts.cluster_method
    k      = opts.n_clusters

    if method == :kmeans
        result = _kmeans_best(zscored', k, opts)
        return assignments(result), result.totalcost

    elseif method == :kmedoids
        D      = _pairwise_euclidean(zscored)
        result = kmedoids(D, k; maxiter = opts.kmeans_max_iters, tol = opts.kmeans_tol)
        ids    = assignments(result)
        return ids, result.totalcost

    elseif method == :hclust
        D      = _pairwise_euclidean(zscored)
        hc     = hclust(D; linkage = opts.cluster_hclust_linkage)
        ids    = cutree(hc; k)
        wcss   = sum(sum((zscored[ids .== c, :] .- mean(zscored[ids .== c, :], dims=1)).^2)
                     for c in unique(ids))
        return ids, wcss

    elseif method == :dbscan
        result = dbscan(zscored', opts.cluster_dbscan_eps;
                        min_neighbors = opts.cluster_dbscan_minpts)
        ids    = assignments(result)   # 0 = noise
        return ids, 0.0

    else
        error("Unknown cluster_method: $method. Choose :kmeans, :kmedoids, :hclust, or :dbscan.")
    end
end

# Trend-test variant: flag flat curves, run the chosen algorithm on the dynamic
# subset only, assign flat curves to the sentinel label n_clusters.
function _cluster_with_trend_method(
    curves::Matrix{Float64},
    times::Vector{Float64},
    zscored_all::Matrix{Float64},
    opts::FitOptions,
)::Tuple{Vector{Int}, Float64}
    opts.n_clusters <= 1 && return _cluster_dispatch(zscored_all, opts)

    W           = size(curves, 1)
    flat_mask   = _flat_curve_mask(curves, times; p_threshold = opts.cluster_trend_p_thr)
    dynamic_idx = findall(.!flat_mask)

    any(flat_mask) || return _cluster_dispatch(zscored_all, opts)

    labels      = fill(opts.n_clusters, W)
    wcss        = 0.0

    if !isempty(dynamic_idx) && opts.n_clusters > 1
        sub_opts = FitOptions(; (f => getfield(opts, f) for f in fieldnames(FitOptions))...,
                                n_clusters                 = opts.n_clusters - 1,
                                cluster_trend_test         = false,
                                cluster_prescreen_constant = false)
        sub_z = zscored_all[dynamic_idx, :]
        sub_labels, wcss = _cluster_dispatch(sub_z, sub_opts)
        for (pos, idx) in enumerate(dynamic_idx)
            labels[idx] = sub_labels[pos]
        end
    end
    return labels, wcss
end

function _cluster_with_non_growing_method(
    zscored_all::Matrix{Float64},
    non_growing::BitVector,
    opts::FitOptions,
)::Tuple{Vector{Int}, Float64}
    opts.n_clusters <= 1 && return _cluster_dispatch(zscored_all, opts)
    any(non_growing) || return _cluster_dispatch(zscored_all, opts)

    dynamic_idx = findall(.!non_growing)
    labels = fill(opts.n_clusters, size(zscored_all, 1))
    wcss = 0.0
    if !isempty(dynamic_idx)
        sub_opts = FitOptions(; (f => getfield(opts, f) for f in fieldnames(FitOptions))...,
            n_clusters=opts.n_clusters - 1,
            cluster_trend_test=false,
            cluster_prescreen_constant=false,
        )
        sub_labels, wcss = _cluster_dispatch(zscored_all[dynamic_idx, :], sub_opts)
        labels[dynamic_idx] = sub_labels
    end
    return labels, wcss
end

# Pre-screening variant: flag constant curves, run the chosen algorithm on the
# dynamic subset only, assign constant curves to the sentinel label n_clusters.
function _cluster_with_prescreen_method(
    curves::Matrix{Float64},
    zscored_all::Matrix{Float64},
    opts::FitOptions,
)::Tuple{Vector{Int}, Float64}
    opts.n_clusters <= 1 && return _cluster_dispatch(zscored_all, opts)

    W           = size(curves, 1)
    const_mask  = _prescreen_constant(curves, opts)
    dynamic_idx = findall(.!const_mask)

    any(const_mask) || return _cluster_dispatch(zscored_all, opts)

    labels      = fill(opts.n_clusters, W)
    wcss        = 0.0

    if !isempty(dynamic_idx) && opts.n_clusters > 1
        sub_opts   = FitOptions(; (f => getfield(opts, f) for f in fieldnames(FitOptions))...,
                                  n_clusters          = opts.n_clusters - 1,
                                  cluster_trend_test  = false,
                                  cluster_prescreen_constant = false)
        sub_z      = zscored_all[dynamic_idx, :]
        sub_labels, wcss = _cluster_dispatch(sub_z, sub_opts)
        for (pos, idx) in enumerate(dynamic_idx)
            labels[idx] = sub_labels[pos]
        end
    end
    return labels, wcss
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
against the distance to its assigned k-means centroid. If the exponential prototype
is closer, reassign the curve to `exp_label`.
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

function _flat_curve_mask(
    curves::Matrix{Float64},
    times::Vector{Float64};
    p_threshold::Float64 = 0.05,
    range_threshold::Union{Nothing,Float64} = nothing,
)::BitVector
    length(times) == size(curves, 2) ||
        throw(ArgumentError("times length $(length(times)) != n_timepoints $(size(curves, 2))"))

    mask = falses(size(curves, 1))
    length(times) < 3 && return mask

    for i in axes(curves, 1)
        finite_mask = isfinite.(curves[i, :]) .& isfinite.(times)
        nf = sum(finite_mask)
        nf < 3 && continue

        y = curves[i, finite_mask]
        t = times[finite_mask]
        if range_threshold !== nothing && maximum(y) - minimum(y) < range_threshold
            mask[i] = true
            continue
        end

        t_centered = t .- mean(t)
        ss_t = sum(t_centered .^ 2)
        if ss_t <= 1e-12
            mask[i] = true
            continue
        end

        slope = sum(t_centered .* y) / ss_t
        y_hat = mean(y) .+ slope .* t_centered
        residuals = y .- y_hat
        s2 = sum(residuals .^ 2) / (nf - 2)
        se_slope = sqrt(max(s2, 0.0) / ss_t)
        if se_slope <= 0 || !isfinite(se_slope)
            mask[i] = abs(slope) < 1e-12
            continue
        end

        t_stat = slope / se_slope
        # Two-tailed p-value from the t-distribution with nf-2 degrees of freedom
        p_value = 2 * (1 - cdf(TDist(nf - 2), abs(t_stat)))
        mask[i] = p_value >= p_threshold
    end
    return mask
end

# Assign a dedicated cluster id to curves with no significant linear trend.
# Uses a two-tailed t-test on the OLS slope (p ≥ 0.05 → flat).
# `flat_id` is passed in by the caller; it must already be within 1..n_clusters.
function _apply_trend_labels(
    curves::Matrix{Float64},
    times::Vector{Float64},
    labels::Vector{Int},
    flat_id::Int,
    p_threshold::Float64 = 0.05,
)::Vector{Int}
    flat_mask = _flat_curve_mask(curves, times; p_threshold=p_threshold)
    new_labels = copy(labels)
    new_labels[flat_mask] .= flat_id
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
function _find_stationary_cutoff(
    data_mat::Matrix{Float64},
    opts::FitOptions;
    reference_mu_max::Union{Nothing, Float64}=nothing,
)::Int
    n = size(data_mat, 2)
    index_od = findall(data_mat[2, :] .> opts.stationary_thr_od)
    isempty(index_od) && return n

    data_t = data_mat[:, index_od]
    sgr = specific_gr_evaluation(data_t, opts.stationary_pt_smooth_derivative)
    mu_max = isnothing(reference_mu_max) ? maximum(sgr) : reference_mu_max
    thr = mu_max * opts.stationary_percentile_thr
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

"""
    find_stationary_cutoff_from_mu(data_mat, mu_max, opts=FitOptions()) -> Int

Find the stationary-phase cutoff using an externally estimated maximum specific
growth rate `mu_max` as the reference for the SGR threshold. The scan and OD
peak snapping are otherwise identical to [`_find_stationary_cutoff`](@ref).

This is intended for workflows that first estimate `mu_max` with
[`fitting_one_well_Log_Lin`](@ref), then derive an empirical carrying capacity
from `data_mat[2, cutoff]`. It does not alter the legacy log-linear result
vector, whose 16th element remains the q95-based `N_max_emp`.
"""
function find_stationary_cutoff_from_mu(
    data_mat::Matrix{Float64},
    mu_max::Real,
    opts::FitOptions=FitOptions(),
)::Int
    mu = Float64(mu_max)
    isfinite(mu) && mu > 0 || throw(ArgumentError("mu_max must be finite and positive"))
    return _find_stationary_cutoff(data_mat, opts; reference_mu_max=mu)
end

# ---------------------------------------------------------------------------
# preprocess for IrregularGrowthData
# ---------------------------------------------------------------------------

"""
    preprocess(data::IrregularGrowthData, opts::FitOptions) -> IrregularGrowthData

Apply the clustering step defined by `opts` to an `IrregularGrowthData`.
Returns a **new** `IrregularGrowthData`; the input is never modified.

Clustering operates on `data.curves` — the resampled matrix on the normalized [0, 1]
union grid built at construction time.

Other preprocessing steps (blank subtraction, negative correction, smoothing, replicate
averaging, scattering correction) are **not supported** for `IrregularGrowthData` because
they require physical time units. Apply them to the raw measurements before constructing
`IrregularGrowthData`. A warning is emitted for each unsupported option that is `true`.
"""
function preprocess(data::IrregularGrowthData, opts::FitOptions)::IrregularGrowthData
    opts.blank_subtraction    && @warn "blank_subtraction is not supported for IrregularGrowthData; apply it on raw data before construction"
    opts.correct_negatives    && @warn "correct_negatives is not supported for IrregularGrowthData; apply it on raw data before construction"
    opts.smooth               && @warn "smooth is not supported for IrregularGrowthData; apply it on raw data before construction"
    opts.average_replicates   && @warn "average_replicates is not supported for IrregularGrowthData"
    opts.scattering_correction && @warn "scattering_correction is not supported for IrregularGrowthData"

    clusters, centroids, wcss = opts.cluster ?
        _cluster(data.curves, data.times, opts) :
        (nothing, nothing, nothing)

    # Re-use the inner constructor directly: raw data and resampled matrix are
    # unchanged; only the clustering fields are updated.
    return IrregularGrowthData(
        data.raw_curves, data.raw_times, data.labels,
        data.curves, data.times,
        clusters, centroids, wcss,
    )
end

export preprocess
export apply_blank_timeseries
export detect_non_growing_indices
export find_stationary_cutoff_from_mu
