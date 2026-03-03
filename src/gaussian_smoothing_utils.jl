# gaussian_smoothing_utils.jl
# Utilities for Gaussian kernel smoothing

function gaussian_kernel(u)
  #= guassian kernel, missing its classical denominator
     because of normalization =#
  return exp(-0.5 * u * u)
end

function compute_bandwidth(t::Vector{Float64}; h_mult::Float64 = 2.0)
    # keep only finite time points
    mask = isfinite.(t)
    t = t[mask]

    # not enough points → use a default bandwidth
    if length(t) < 3
        return 1.0
    end

    # order time points (safety) and estimate a temporal step
    t_sorted = sort(t)
    # most typical temporal step
    dt = median(diff(t_sorted))

    if !isfinite(dt) || dt <= 0.0
        return 1.0
    end

    return h_mult * dt
end

function gaussian_smoothing(t::Vector{Float64}, y::Vector{Float64}, tq::Vector{Float64};
                              h_mult::Float64=2.0)

    mask = isfinite.(t) .& isfinite.(y)
    t2 = t[mask]
    y2 = max.(y[mask], 1e-9)

    if length(t2) == 0
        return fill(0.0, length(tq))
    elseif length(t2) == 1
        return fill(y2[1], length(tq))
    end

    bandwidth = compute_bandwidth(t2; h_mult = h_mult)
    invh = 1.0 / bandwidth

    yhat = Vector{Float64}(undef, length(tq))

    for (j, x) in enumerate(tq)
        ww = gaussian_kernel.((x .- t2) * invh)
        s  = sum(ww)
        if s <= 1e-12 # nearest point fallback
          idx = findmin(abs.(t2 .- x))[2]
          yhat[j] = y2[idx]
        else
          yhat[j] = (ww ⋅ y2) / s
        end
    end

    return yhat
end