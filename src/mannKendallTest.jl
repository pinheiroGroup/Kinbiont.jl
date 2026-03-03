# Stolen Mann Kendall Test, remember to delete this once gone public
using Statistics
using StatsFuns
using Distributions
using StatsBase

function mk_original_test(x; α = 0.05)
    n = length(x)
    s = mk_score(x)
    var_s = variance_s(x,length(x))
    τ = s/(.5*n*(n-1))
    z = z_score(s, var_s)
    
    p, h, trend = p_value(z, α)
    
    return (trend=trend,  h=h, p=p, z=z, τ=τ, s=s, var_s=var_s)
end

function mk_score(x)
    s = 0
    n = length(x)
    
    for i = 1:n-1
        for j = i+1:n
            s += sign(x[j] - x[i])
        end
    end
    
    return s
end

function variance_s(x, n)
    # calculate the unique data
    unique_x = unique(x)
    g = length(unique_x)

    # calculate the var(s)
    if n == g # there is no tie
        var_s = (n*(n-1)*(2*n+5))/18
        
    else # there are some ties in data
        tp = collect(values(countmap(x)))
        var_s = (n * (n - 1) * (2 * n + 5) - sum(tp .* (tp .- 1) .* (2 .* tp .+ 5)))/18
    end
        
    return var_s
end

function z_score(s, var_s)
    z = s == 0 ? 0 : (s - sign(s)) / sqrt(var_s)
    return z
end

function p_value(z, α)
    # two tail test
    p = 2*(1-normcdf(abs(z)))  
    norminv = quantile(Normal(), 1 - α/2)
    h = abs(z) > norminv

    if (z < 0) && h
        trend = "decreasing"
    elseif (z > 0) && h
        trend = "increasing"
    else
        trend = "no trend"
    end
    
    return p, h, trend
end