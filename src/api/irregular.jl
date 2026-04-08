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

"""
    IrregularGrowthData

Container for growth curves measured at per-curve irregular time points.

# Fields
- `raw_curves::Vector{Vector{Float64}}`: original OD values, one vector per curve.
- `raw_times::Vector{Vector{Float64}}`: original (unnormalized) time points, one vector per curve.
- `labels::Vector{String}`: identifier for each curve.
- `curves::Matrix{Float64}`: `n_curves × n_grid` resampled matrix on the normalized [0, 1] union grid. Ready for clustering.
- `times::Vector{Float64}`: the [0, 1] union grid shared by all resampled curves.
- `clusters`, `centroids`, `wcss`: populated by [`preprocess`](@ref), `nothing` until then.

Construct with:
```julia
data = IrregularGrowthData(raw_curves, raw_times, labels; step=0.01)
```
where `raw_curves` and `raw_times` are `Vector{Vector{Float64}}` of the same length,
and `step` controls the resolution of the union grid in normalized time.
"""
struct IrregularGrowthData
    raw_curves :: Vector{Vector{Float64}}
    raw_times  :: Vector{Vector{Float64}}
    labels     :: Vector{String}
    curves     :: Matrix{Float64}
    times      :: Vector{Float64}
    clusters   :: Union{Nothing, Vector{Int}}
    centroids  :: Union{Nothing, Matrix{Float64}}
    wcss       :: Union{Nothing, Float64}

    # Inner constructor: validates all fields and stores directly.
    # Called by the outer constructor (after resampling) and by preprocess
    # (when updating only clustering fields).
    function IrregularGrowthData(
        raw_curves, raw_times, labels, curves, times,
        clusters = nothing, centroids = nothing, wcss = nothing,
    )
        n = length(raw_curves)
        length(raw_times)  == n ||
            throw(ArgumentError("raw_times length $(length(raw_times)) ≠ n_curves $n"))
        length(labels)     == n ||
            throw(ArgumentError("labels length $(length(labels)) ≠ n_curves $n"))
        for i in 1:n
            length(raw_curves[i]) == length(raw_times[i]) ||
                throw(ArgumentError("curve $i: raw_curves length $(length(raw_curves[i])) ≠ raw_times length $(length(raw_times[i]))"))
            length(raw_times[i]) >= 2 ||
                throw(ArgumentError("curve $i: time vector must have at least 2 points"))
        end
        size(curves, 1) == n ||
            throw(ArgumentError("curves rows $(size(curves,1)) ≠ n_curves $n"))
        size(curves, 2) == length(times) ||
            throw(ArgumentError("curves cols $(size(curves,2)) ≠ length(times) $(length(times))"))
        clusters === nothing || length(clusters) == n ||
            throw(ArgumentError("clusters length $(length(clusters)) ≠ n_curves $n"))
        new(raw_curves, raw_times, labels, curves, times, clusters, centroids, wcss)
    end
end

# Outer constructor: normalizes times → builds union grid → resamples → stores.
function IrregularGrowthData(
    raw_curves :: Vector{Vector{Float64}},
    raw_times  :: Vector{Vector{Float64}},
    labels     :: Vector{String};
    step       :: Float64 = 0.01,
)
    n = length(raw_curves)
    # Early validation so errors point here, not at the inner constructor
    length(raw_times) == n ||
        throw(ArgumentError("raw_times length $(length(raw_times)) ≠ n_curves $n"))
    length(labels)    == n ||
        throw(ArgumentError("labels length $(length(labels)) ≠ n_curves $n"))
    for i in 1:n
        length(raw_curves[i]) == length(raw_times[i]) ||
            throw(ArgumentError("curve $i: raw_curves and raw_times must have the same length"))
        length(raw_times[i]) >= 2 ||
            throw(ArgumentError("curve $i: time vector must have at least 2 points"))
    end

    times01_list = [_normalize_times_01(raw_times[i]) for i in 1:n]
    union_grid   = _build_union_grid(times01_list; step = step)
    curves_mat   = _resample_to_union_grid(times01_list, raw_curves, union_grid; step = step)

    return IrregularGrowthData(raw_curves, raw_times, labels, curves_mat, union_grid)
end

"""
    data[labels::Vector{String}] -> IrregularGrowthData

Return a new `IrregularGrowthData` containing only the curves whose labels appear in
`labels`, in the order given. Clustering fields (`clusters`, `centroids`,
`wcss`) are not carried over — subsetting invalidates any prior clustering.

Throws `ArgumentError` if any requested label is not present in `data.labels`.

# Example
```julia
data = IrregularGrowthData(raw_curves, raw_times, labels)  # load data
sub = data[["A"]]                                          # single curve
sub_multi = data[["A", "C"]]                               # two curves
```
"""
function Base.getindex(data::IrregularGrowthData, labels::Vector{String})
    missing_labels = filter(l -> !(l in data.labels), labels)
    isempty(missing_labels) || throw(ArgumentError(
        "label(s) not found in IrregularGrowthData: $(join(missing_labels, ", "))"
    ))
    idx = [findfirst(==(l), data.labels) for l in labels]
    return IrregularGrowthData(
        data.raw_curves[idx],
        data.raw_times[idx],
        data.labels[idx],
        data.curves[idx, :],
        data.times,       # shared union grid — no re-resampling needed
    )
end
