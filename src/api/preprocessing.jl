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
2. Blank subtraction (`opts.blank_subtraction`)
3. Negative-value correction (`opts.correct_negatives`)
4. Smoothing (`opts.smooth`)
5. K-means clustering on z-scored curves (`opts.cluster`)

# Example
```julia
opts = FitOptions(smooth=true, cluster=true, n_clusters=4)
processed = preprocess(raw_data, opts)
# processed.clusters now holds a cluster id for every curve
```
"""
function preprocess(data::GrowthData, opts::FitOptions)::GrowthData
    curves = data.curves   # n_curves × n_tp, never mutated in place

    curves = _apply_scattering_correction(curves, data.times, opts)
    curves = _apply_blank_subtraction(curves, opts)
    curves = _apply_negative_correction(curves, data.times, opts)
    curves = _apply_smoothing(curves, data.times, opts)

    clusters = opts.cluster ? _cluster(curves, data.times, opts) : nothing

    return GrowthData(curves, data.times, data.labels, clusters)
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
# ---------------------------------------------------------------------------

function _apply_smoothing(
    curves::Matrix{Float64},
    times::Vector{Float64},
    opts::FitOptions,
)::Matrix{Float64}
    opts.smooth || return curves
    opts.smooth_method == :none && return curves

    smoothing_str = _smoothing_symbol_to_string(opts.smooth_method)
    smoothed = similar(curves)
    for i in axes(curves, 1)
        curve_mat = Matrix(transpose(hcat(times, curves[i, :])))
        result = smoothing_data(
            curve_mat;
            method = smoothing_str,
            pt_avg = opts.smooth_pt_avg,
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
    return smoothed
end

_smoothing_symbol_to_string(s::Symbol) = Dict(
    :lowess     => "lowess",
    :rolling_avg => "rolling_avg",
    :none       => "NO",
)[s]

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
