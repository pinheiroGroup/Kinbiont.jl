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
    clusters = opts.cluster ? _cluster(curves, times, opts) : nothing

    curves = _apply_blank_subtraction(curves, opts)
    curves = _apply_negative_correction(curves, times, opts)
    curves, times = _apply_smoothing(curves, times, opts)   # Gaussian may change times

    return GrowthData(curves, times, data.labels, clusters)
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
    smoothed = similar(curves)
    for i in axes(curves, 1)
        curve_mat = Matrix(transpose(hcat(times, curves[i, :])))
        result = smoothing_data(
            curve_mat;
            method     = smoothing_str,
            pt_avg     = opts.smooth_pt_avg,
            thr_lowess = opts.lowess_frac,
        )
        # rolling_avg reduces length; pad with original tail if needed
        n_out = size(result, 2)
        if n_out == length(times)
            smoothed[i, :] = result[2, :]
        else
            smoothed[i, 1:n_out] = result[2, :]
            smoothed[i, n_out+1:end] = curves[i, n_out+1:end]
        end
    end
    return smoothed, times
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
    _cluster(curves, times, opts) -> Vector{Int}

Cluster growth curves using k-means on z-score-normalised data.
Flat / non-growing curves are identified via a Mann-Kendall trend test and
may be assigned to their own cluster regardless of geometric distance.
"""
function _cluster(
    curves::Matrix{Float64},
    times::Vector{Float64},
    opts::FitOptions,
)::Vector{Int}
    n_curves = size(curves, 1)

    # Z-score each curve independently (StatsBase.zscore over columns = time axis)
    # zscore(x, 2) normalises each row (curve) across time points
    zscored = _zscore_rows(curves)

    # k-means needs features × samples layout
    result = kmeans(zscored', opts.n_clusters)
    labels = assignments(result)   # Vector{Int}, length n_curves

    # Optionally re-label structurally flat curves using Mann-Kendall
    if opts.cluster_trend_test
        labels = _apply_trend_labels(curves, times, labels, opts.n_clusters)
    end

    return labels
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
function _apply_trend_labels(
    curves::Matrix{Float64},
    times::Vector{Float64},
    labels::Vector{Int},
    k::Int,
)::Vector{Int}
    flat_id   = k + 1
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

export preprocess
