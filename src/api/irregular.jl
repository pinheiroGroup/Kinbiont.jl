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

function _interp_linear(
    x::Vector{Float64},
    y::Vector{Float64},
    x_new::Vector{Float64},
)::Vector{Float64}
    n   = length(x)
    out = Vector{Float64}(undef, length(x_new))
    for (k, xi) in enumerate(x_new)
        if xi <= x[1]
            out[k] = y[1]
        elseif xi >= x[end]
            out[k] = y[end]
        else
            hi = searchsortedfirst(x, xi)   # first index where x[hi] >= xi
            lo = hi - 1
            α  = (xi - x[lo]) / (x[hi] - x[lo])
            out[k] = (1.0 - α) * y[lo] + α * y[hi]
        end
    end
    return out
end

function _resample_to_union_grid(
    times01_list::Vector{Vector{Float64}},
    values_list::Vector{Vector{Float64}},
    union_grid::Vector{Float64};
    step::Float64 = 0.01,
)::Matrix{Float64}
    n = length(times01_list)
    X = Matrix{Float64}(undef, n, length(union_grid))

    for i in 1:n
        t_raw = times01_list[i]
        y_raw = values_list[i]

        # Snap normalized times to the same grid step used in build_union_grid
        t_snapped = [clamp(round(v / step) * step, 0.0, 1.0) for v in t_raw]

        # Deduplicate: keep the first occurrence of each snapped time
        seen   = Dict{Float64, Bool}()
        keep   = Vector{Bool}(undef, length(t_snapped))
        for (j, tv) in enumerate(t_snapped)
            keep[j] = !haskey(seen, tv)
            seen[tv] = true
        end
        t_u = t_snapped[keep]
        y_u = y_raw[keep]

        # Sort by time (should already be sorted for regular inputs)
        perm = sortperm(t_u)
        t_u  = t_u[perm]
        y_u  = y_u[perm]

        X[i, :] = _interp_linear(t_u, y_u, union_grid)
    end
    return X
end
