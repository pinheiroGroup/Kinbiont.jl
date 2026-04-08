# =============================================================================
# IrregularGrowthData — one time vector per curve, normalized to [0, 1]
# =============================================================================

export IrregularGrowthData

function _normalize_times_01(times::Vector{Float64})::Vector{Float64}
    tmin = minimum(times)
    tmax = maximum(times)
    span = tmax - tmin
    span <= 0.0 && return zeros(Float64, length(times))
    return (times .- tmin) ./ span
end
