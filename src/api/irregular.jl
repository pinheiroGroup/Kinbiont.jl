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

function _build_union_grid(
    times01_list::Vector{Vector{Float64}};
    step::Float64 = 0.01,
)::Vector{Float64}
    points = Set{Float64}()
    for t in times01_list
        for v in t
            snapped = round(v / step) * step
            push!(points, clamp(snapped, 0.0, 1.0))
        end
    end
    grid = sort!(collect(points))
    isempty(grid) || grid[1] > 0.0  && pushfirst!(grid, 0.0)
    isempty(grid) || grid[end] < 1.0 && push!(grid, 1.0)
    return grid
end
