
function loss_L2_derivative(data, ODE_prob, integrator, p, tsteps)
    sol = solve(
        ODE_prob,
        integrator,
        p=p,
        saveat=tsteps,
        verbose=false,
        abstol=1e-10,
        reltol=1e-10,
    )
    sol_t = reduce(hcat, sol.u)
    sol_time = reduce(hcat, sol.t)
    sol_t = reduce(hcat, sol.u)
    sol_t = sum(sol_t, dims=1)

    itp = interpolate((sol_time,), sol_t, Gridded(Linear()))
    derivative_theo = only.(Interpolations.gradient.(Ref(itp), sol_t))

    itp = interpolate((data[1, :],), data[2, :], Gridded(Linear()))
    derivative_data = only.(Interpolations.gradient.(Ref(itp), data[2, :]))

    if size(derivative_theo)[1] == size(derivative_data)[1]
        lossa =
            NaNMath.sum(abs2.((derivative_theo[2, :] .- derivative_data[1, 1:end]))) /
            length(derivative_data[2, :])
    else
        lossa = 10.0^9 * length(data[2, :])
    end

    return lossa, sol
end

function loss_blank_weighted_L2(data, ODE_prob, integrator, p, tsteps, blank_array)
    sol = solve(
        ODE_prob,
        integrator,
        p=p,
        saveat=tsteps,
        verbose=false,
        abstol=1e-10,
        reltol=1e-10,
    )
    sol_t = reduce(hcat, sol.u)
    sol_t = sum(sol_t, dims=1)

    empirical_blank_distrib = (blank_array .- mean(blank_array))
    # generation of the empirica distrib respect the mean of the noise
    test = StatsBase.fit(Histogram, empirical_blank_distrib)
    binning_distrib = test.edges
    probabiltity_distrib = test.weights ./ sum(test.weights)

    if size(sol_t)[2] == size(data)[2]
        lossa = 0.0

        for ll = 1:size(sol_t)[2]
            dist = (data[2, ll] - sol_t[1, ll])
            index = findfirst(x -> (x > dist), binning_distrib[1])
            if (typeof(index) == Nothing || index > length(probabiltity_distrib))
                lossa = lossa + abs2.(dist) / length(data[2, :])
            else
                prob = probabiltity_distrib[index]
                lossa = lossa + abs2.((1 - prob) * dist) / length(data[2, :])
            end
        end
    else
        lossa = 10.0^9 * length(data[2, :])
    end

    return lossa, sol
end

function loss_L2(data, ODE_prob, integrator, p, tsteps)
    sol = solve(
        ODE_prob,
        integrator,
        p=p,
        saveat=tsteps,
        verbose=false,
        abstol=1e-10,
        reltol=1e-10,
    )
    sol_t = reduce(hcat, sol.u)
    sol_t = sum(sol_t, dims=1)

    if size(sol_t)[2] == size(data)[2]
        lossa = NaNMath.sum(abs2.((data[2, :] .- sol_t[1, 1:end]))) / length(data[2, :])
    else
        lossa = 10.0^9 * length(data[2, :])
    end

    return lossa, sol
end

function loss_L2_log(data, ODE_prob, integrator, p, tsteps)
    sol = solve(
        ODE_prob,
        integrator,
        p=p,
        saveat=tsteps,
        verbose=false,
        abstol=1e-10,
        reltol=1e-10,
    )
    sol_t = reduce(hcat, sol.u)
    sol_t = sum(sol_t, dims=1)

    if size(sol_t)[2] == size(data)[2]
        lossa = log(NaNMath.sum(abs2.((data[2, :] .- sol_t[1, 1:end]))) / length(data[2, :]))
    else
        lossa = 10.0^9 * length(data[2, :])
    end

    return lossa, sol
end


function loss_RE(data, ODE_prob, integrator, p, tsteps)
    sol = solve(
        ODE_prob,
        integrator,
        p=p,
        saveat=tsteps,
        verbose=false,
        abstol=1e-10,
        reltol=1e-10,
    )
    sol_t = reduce(hcat, sol.u)
    sol_t = sum(sol_t, dims=1)

    if size(sol_t)[2] == size(data)[2]
        lossa =
            0.5 * NaNMath.sum(abs2.(1.0 .- (data[2, :] ./ sol_t[1, 1:end]))) /
            length(data[2, :])
    else
        lossa = 10.0^9 * length(data[2, :])
    end

    return lossa, sol
end



function loss_RE_log(data, ODE_prob, integrator, p, tsteps)
    sol = solve(
        ODE_prob,
        integrator,
        p=p,
        saveat=tsteps,
        verbose=false,
        abstol=1e-10,
        reltol=1e-10,
    )
    sol_t = reduce(hcat, sol.u)
    sol_t = sum(sol_t, dims=1)

    if size(sol_t)[2] == size(data)[2]
        lossa = log(
            0.5 * NaNMath.sum(abs2.(1.0 .- (data[2, :] ./ sol_t[1, 1:end]))) /
            length(data[2, :]))
    else
        lossa = 10.0^9 * length(data[2, :])
    end

    return lossa, sol
end


function loss_L2_std_blank(data, ODE_prob, integrator, p, tsteps, blank_array)

    # evaluation of std of  empirica distrib of the blank

    std_blank = Statistics.std(blank_array)
    if std_blank == 0.0
        std_blank = 1.0
    end
    sol = solve(
        ODE_prob,
        integrator,
        p=p,
        saveat=tsteps,
        verbose=false,
        abstol=1e-10,
        reltol=1e-10,
    )
    sol_t = reduce(hcat, sol.u)
    sol_t = sum(sol_t, dims=1)

    if size(sol_t)[2] == size(data)[2]
        lossa = NaNMath.sum(abs2.((data[2, :] .- sol_t[1, 1:end]) ./ std_blank)) / length(data[2, :])
    else
        lossa = 10.0^9 * length(data[2, :])
    end

    return lossa, sol
end







function select_loss_function(loss_name, data, ODE_prob, integrator, tsteps, blank_array)
    loss_functions = Dict(
        "L2" => loss_L2,
        "RE" => loss_RE,
        "blank_weighted_L2" => loss_blank_weighted_L2,
        "L2_derivative" => loss_L2_derivative,
        "L2_std_blank" => loss_L2_std_blank,
        "L2_log" => loss_L2_log,
        "RE_log" => loss_RE_log,)

    if loss_name == "blank_weighted_L2" || loss_name == "L2_std_blank"
        return (p) ->
            loss_functions[loss_name](data, ODE_prob, integrator, p, tsteps, blank_array)
    else
        return (p) -> loss_functions[loss_name](data, ODE_prob, integrator, p, tsteps)
    end
end

export loss_L2_derivative
export loss_blank_weighted_L2
export loss_L2
export loss_L2_log
export loss_RE
export loss_RE_log
export loss_L2_std_blank
export select_loss_export
