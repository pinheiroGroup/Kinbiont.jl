# =============================================================================
# Unified Preprocessing Pipeline
# =============================================================================
# preprocess(data::GrowthData, opts::FitOptions) -> GrowthData
#
# Pure function — always returns a new GrowthData, never mutates the input.
# Delegates all actual computation to the existing functions already in
# pre_processing_functions.jl, and uses:
#   - StatsBase.zscore       for z-score normalisation
#   - Clustering.kmeans      for k-means (default, cluster_engine=:clustering_jl)
#   - _km_parallel           for k-means (opt-in, cluster_engine=:parallel_julia)
#   - linear slope t-test (Statistics stdlib) for flat-curve trend detection
# =============================================================================

using Clustering: kmeans, assignments
using StatsBase: zscore
using Distributions: Normal, cdf

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

    curves = _apply_scattering_correction(curves, times, opts)

    # Clustering is performed on scattering-corrected but otherwise raw data so
    # that blank subtraction does not hide the distinction between growing and
    # non-growing wells (which is precisely what clustering is meant to find).
    clusters, centroids, wcss = opts.cluster ? _cluster(curves, times, opts) : (nothing, nothing, nothing)

    curves = _apply_blank_subtraction(curves, opts)
    curves = _apply_negative_correction(curves, times, opts)
    curves, times = _apply_smoothing(curves, times, opts)   # Gaussian may change times

    return GrowthData(curves, times, data.labels, clusters, centroids, wcss)
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
    opts::FitOptions,
)::Matrix{Float64}
    opts.blank_subtraction || return curves
    return curves .- opts.blank_value
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

    if opts.smooth_method == :gaussian
        return _apply_gaussian_smoothing(curves, times, opts)
    end

    smoothing_str = _smoothing_symbol_to_string(opts.smooth_method)
    n_curves = size(curves, 1)

    # Process first curve to determine output length (smoothing may shorten data)
    curve_mat1 = Matrix(transpose(hcat(times, curves[1, :])))
    result1 = smoothing_data(
        curve_mat1;
        method     = smoothing_str,
        pt_avg     = opts.smooth_pt_avg,
        thr_lowess = opts.lowess_frac,
    )
    n_out     = size(result1, 2)
    new_times = Vector{Float64}(result1[1, :])
    smoothed  = Matrix{Float64}(undef, n_curves, n_out)
    smoothed[1, :] = result1[2, :]

    for i in 2:n_curves
        curve_mat = Matrix(transpose(hcat(times, curves[i, :])))
        result = smoothing_data(
            curve_mat;
            method     = smoothing_str,
            pt_avg     = opts.smooth_pt_avg,
            thr_lowess = opts.lowess_frac,
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
    :none        => "NO",
)[s]

# ---------------------------------------------------------------------------
# Gaussian kernel smoothing (no external dependencies)
# ---------------------------------------------------------------------------

"""
    _gaussian_kernel(u) -> Float64

Un-normalised Gaussian kernel evaluated at standardised distance `u = (t - tᵢ)/h`.
"""
_gaussian_kernel(u::Float64) = exp(-0.5 * u * u)

"""
    _gaussian_bandwidth(t; h_mult) -> Float64

Estimate bandwidth as `h_mult × median(Δt)` over the finite, sorted time points.
Falls back to 1.0 when fewer than three points are available or the median step
is non-positive.
"""
function _gaussian_bandwidth(t::Vector{Float64}; h_mult::Float64 = 2.0)
    t_finite = filter(isfinite, t)
    length(t_finite) < 3 && return 1.0
    dt = median(diff(sort(t_finite)))
    (isfinite(dt) && dt > 0.0) || return 1.0
    return h_mult * dt
end

"""
    _gaussian_smooth_curve(t, y, tq; h_mult) -> Vector{Float64}

Nadaraya–Watson kernel smoother with a Gaussian kernel. Evaluates the smoothed
curve at the query points `tq`. Non-finite values in `(t, y)` are excluded.
Falls back to the nearest observed value when the total kernel weight is < 1e-12.
"""
function _gaussian_smooth_curve(
    t::Vector{Float64},
    y::Vector{Float64},
    tq::Vector{Float64};
    h_mult::Float64 = 2.0,
)::Vector{Float64}
    mask = isfinite.(t) .& isfinite.(y)
    t2 = t[mask]
    y2 = max.(y[mask], 1e-9)     # clamp to small positive (matches New-api)

    length(t2) == 0 && return fill(0.0, length(tq))
    length(t2) == 1 && return fill(y2[1], length(tq))

    h    = _gaussian_bandwidth(t2; h_mult)
    invh = 1.0 / h
    yhat = Vector{Float64}(undef, length(tq))

    for (j, x) in enumerate(tq)
        ww = _gaussian_kernel.((x .- t2) .* invh)
        s  = sum(ww)
        if s <= 1e-12
            # nearest-point fallback
            yhat[j] = y2[argmin(abs.(t2 .- x))]
        else
            yhat[j] = sum(ww .* y2) / s
        end
    end
    return yhat
end

function _apply_gaussian_smoothing(
    curves::Matrix{Float64},
    times::Vector{Float64},
    opts::FitOptions,
)::Tuple{Matrix{Float64}, Vector{Float64}}
    tq      = something(opts.gaussian_time_grid, times)   # query grid
    n_curves = size(curves, 1)
    smoothed = Matrix{Float64}(undef, n_curves, length(tq))

    for i in axes(curves, 1)
        smoothed[i, :] = _gaussian_smooth_curve(
            times, curves[i, :], tq; h_mult = opts.gaussian_h_mult
        )
    end
    return smoothed, tq
end

# ---------------------------------------------------------------------------
# Step 5 — K-means clustering on z-scored curves
# ---------------------------------------------------------------------------

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

Three modes (evaluated in priority order):
1. `cluster_prescreen_constant=true`: quantile-ratio pre-screening identifies
   non-growing wells before k-means, which then runs on dynamic curves only.
2. `cluster_trend_test=true`: OLS slope t-test post-hoc re-labels flat curves.
3. Neither: plain k-means on all curves.
"""
function _cluster(
    curves::Matrix{Float64},
    times::Vector{Float64},
    opts::FitOptions,
)::Tuple{Vector{Int}, Matrix{Float64}, Float64}
    # Z-score all curves once; used for both k-means and centroid computation.
    zscored_all = _zscore_rows(curves)

    if opts.cluster_prescreen_constant
        labels, wcss = _cluster_with_prescreen(curves, zscored_all, opts)
    elseif opts.cluster_trend_test
        k_dynamic      = max(1, opts.n_clusters - 1)
        labels, wcss   = _run_kmeans(zscored_all, k_dynamic, opts)
        labels         = _apply_trend_labels(curves, times, labels, opts.n_clusters)
    else
        labels, wcss   = _run_kmeans(zscored_all, opts.n_clusters, opts)
    end

    # Centroids in z-normalised space (shape prototypes, scale-independent).
    centroids = _compute_centroids(zscored_all, labels, opts.n_clusters)
    return labels, centroids, wcss
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
        km_labels, wcss = _run_kmeans(zscored_all[dynamic_idx, :], k_dynamic, opts)
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

Compute the mean curve (in original space) for each cluster.
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
# Uses a t-test on the OLS slope (p ≥ 0.05 → flat), implemented with Statistics
# stdlib only — no extra package required.
# `flat_id` is passed in by the caller; it must already be within 1..n_clusters.
function _apply_trend_labels(
    curves::Matrix{Float64},
    times::Vector{Float64},
    labels::Vector{Int},
    flat_id::Int,
)::Vector{Int}
    new_labels = copy(labels)
    n          = length(times)
    t_centered = times .- mean(times)
    ss_t       = sum(t_centered .^ 2)

    for i in axes(curves, 1)
        y         = curves[i, :]
        slope     = sum(t_centered .* y) / ss_t
        y_hat     = mean(y) .+ slope .* t_centered
        residuals = y .- y_hat
        s2        = sum(residuals .^ 2) / (n - 2)
        se_slope  = sqrt(s2 / ss_t)
        t_stat    = slope / se_slope
        # Two-tailed p-value via normal approximation (good for n > 10)
        p_approx  = 2 * (1 - cdf(Normal(), abs(t_stat)))
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

# ---------------------------------------------------------------------------
# Optional parallel k-means engine  (cluster_engine = :parallel_julia)
# ---------------------------------------------------------------------------
# A pure-Julia multi-threaded k-means implementation.  Assignment and centroid-
# update steps use Base.Threads.@threads; the outer restart loop is sequential.
# The interface mirrors Clustering.kmeans: returns (labels, totalcost).

@inline function _km_sqdist(u, v)
    s = 0.0
    @inbounds @simd for j in eachindex(u, v)
        d = u[j] - v[j]; s += d * d
    end
    s
end

function _km_init_random(X::Matrix{Float64}, k::Int, rng::AbstractRNG)
    [copy(X[rand(rng, 1:size(X, 1)), :]) for _ in 1:k]
end

function _km_assign!(labels::Vector{Int}, X::Matrix{Float64},
                     centroids::Vector{Vector{Float64}})::Float64
    n = size(X, 1); k = length(centroids)
    inertia = zeros(Float64, Threads.nthreads())
    Threads.@threads for i in 1:n
        xi = @view X[i, :]
        best_k = 1
        best_d = _km_sqdist(xi, centroids[1])
        @inbounds for j in 2:k
            d = _km_sqdist(xi, centroids[j])
            if d < best_d; best_d = d; best_k = j; end
        end
        labels[i] = best_k
        inertia[Threads.threadid()] += best_d
    end
    sum(inertia)
end

function _km_update!(centroids::Vector{Vector{Float64}},
                     X::Matrix{Float64}, labels::Vector{Int})
    n, m = size(X); k = length(centroids); T = Threads.nthreads()
    sums = [zeros(Float64, m, k) for _ in 1:T]
    cnts = [zeros(Int, k)        for _ in 1:T]
    Threads.@threads for i in 1:n
        t = Threads.threadid(); lbl = labels[i]
        @views sums[t][:, lbl] .+= X[i, :]
        cnts[t][lbl] += 1
    end
    S = sum(sums); C = sum(cnts)
    for j in 1:k
        if C[j] == 0
            centroids[j] .= X[rand(1:n), :]
        else
            centroids[j] .= S[:, j] ./ C[j]
        end
    end
end

"""
    _km_parallel(X, k; n_init, rng) -> (labels::Vector{Int}, totalcost::Float64)

Pure-Julia multi-threaded k-means on the rows of `X`.  Runs `n_init` random
restarts and returns the assignment with the lowest total SSE.  Convergence is
declared when the SSE change between iterations falls below `1e-8`.
"""
function _km_parallel(
    X::Matrix{Float64},
    k::Int;
    n_init::Int        = 3,
    max_iter::Int      = 300,
    tol::Float64       = 1e-8,
    rng::AbstractRNG   = MersenneTwister(42),
)::Tuple{Vector{Int}, Float64}
    n = size(X, 1)
    best_sse    = Inf
    best_labels = zeros(Int, n)
    labels      = zeros(Int, n)

    for _ in 1:n_init
        centroids = _km_init_random(X, k, rng)
        prev_sse  = Inf
        for _ in 1:max_iter
            sse = _km_assign!(labels, X, centroids)
            _km_update!(centroids, X, labels)
            abs(prev_sse - sse) < tol && break
            prev_sse = sse
        end
        sse = _km_assign!(labels, X, centroids)
        if sse < best_sse
            best_sse = sse
            best_labels .= labels
        end
    end
    return best_labels, best_sse
end

# Thin dispatch wrapper: calls the right engine and returns (labels, totalcost).
function _run_kmeans(
    X_zscored::Matrix{Float64},
    k::Int,
    opts::FitOptions,
)::Tuple{Vector{Int}, Float64}
    if opts.cluster_engine === :parallel_julia
        return _km_parallel(X_zscored, k; n_init = opts.cluster_n_init)
    else
        result = kmeans(X_zscored', k)
        return assignments(result), result.totalcost
    end
end

export preprocess
