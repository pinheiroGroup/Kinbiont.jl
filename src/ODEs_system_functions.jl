
function vectorize_df_results_ODE_System(label_exp, model, res, loss_value)

    if typeof(model) == String

        param_names = ["label_exp" "model" reduce(hcat, ODEs_system_models[model].params) "loss"]

        param_values = [label_exp model reduce(hcat, res) loss_value]


    else

        model_string = "custom"

        nmax_param = size(res)[1]




        param_names_t = reduce(hcat, [string("param_", i) for i in 1:nmax_param])


        param_names = ["label_exp" "model" param_names_t "loss"]

        param_values = [label_exp model_string reduce(hcat, res) loss_value]

    end


    res_matrix = vcat(param_names, param_values)



    return res_matrix

end



function ODEs_system_model_selector2(model, u0, tspan; param=nothing)

    if typeof(model) == String

        model_string = ODEs_system_models[model].name
        model_function = ODEs_system_models[model].func

        ODE_prob = ODEProblem(model_function, u0, tspan, param)


    else

        model_string = "custom"

        ODE_prob = ODEProblem(model, u0, tspan, param)

    end


    return ODE_prob
end

function loss_RE_ODE_Sys(data, index_of_eqs, index_of_data, ODE_prob, Integration_method, p, tsteps)
    sol = solve(
        ODE_prob,
        Integration_method,
        p=p,
        saveat=tsteps,
        verbose=false,
        abstol=1e-10,
        reltol=1e-10,
    )


    sol_t = reduce(hcat, sol.u)



    if size(sol_t)[2] == size(data)[2]

        lossa =
            0.5 * NaNMath.sum(abs2.(1.0 .- (data[index_of_data, :] ./ sol_t[index_of_eqs, :]))) /
            length(data[2, :])

    else

        lossa = 10.0^9 * length(data[2, :])
    end


    return lossa, sol
end


function define_loss_function_odes_system(data, set_of_equation_to_fit, ODE_prob, Integration_method, tsteps)


    if !isnothing(set_of_equation_to_fit)
        index_of_eqs = set_of_equation_to_fit 
        index_of_data = 1:1:(size(set_of_equation_to_fit)[1])
        index_of_data = index_of_data .+1


    else

        index_of_eqs = 1:1:(size(data)[1]-1)
        index_of_data = index_of_eqs .+ 1

    end

    function loss_relative_error(data, index_of_eqs, index_of_data, ODE_prob, Integration_method, p, tsteps)
        sol = solve(
            ODE_prob,
            Integration_method,
            p=p,
            saveat=tsteps,
            verbose=false,
            abstol=1e-10,
            reltol=1e-10,
        )
        sol_t = reduce(hcat, sol.u)
        #sol_t = sum(sol_t, dims=1)
    
        if size(sol_t)[2] == size(data)[2]
            # Calculate relative error: |data - model| / |data|
            # Add small epsilon to avoid division by zero
            epsilon = 1e-10
            relative_errors = abs.((data[index_of_data, :] .- sol_t[index_of_eqs, :])) ./ (abs.(data[index_of_data, :]) .+ epsilon)
            lossa = NaNMath.sum(relative_errors) / length(data[2, :])
        else
            lossa = 10.0^9 * length(data[2, :])
        end
    
        return lossa, sol
    end


    function loss_L2_ODE_Sys(data, index_of_eqs, index_of_data, ODE_prob, Integration_method, p, tsteps)
        sol = solve(
            ODE_prob,
            Integration_method,
            p=p,
            saveat=tsteps,
            verbose=false,
            abstol=1e-10,
            reltol=1e-10,
        )


        sol_t = reduce(hcat, sol.u)



        if size(sol_t)[2] == size(data)[2]

            lossa =NaNMath.sum(abs2.( (data[index_of_data, :] .- sol_t[index_of_eqs, :]))) / length(data[2, :])

        else

            lossa = 10.0^9 * length(data[2, :])
        end


        return lossa, sol
    end



    return (p) ->
    loss_L2_ODE_Sys(data, index_of_eqs, index_of_data, ODE_prob, Integration_method, p, tsteps)

end


"""
    ODEs_system_sim(
    model::Any,
    start_IC::Vector{Float64},
    tstart::Float64, 
    tmax::Float64,
    delta_t::Float64,
    param_of_ode; 
    Integration_method=Tsit5()
    )

Simulates an ordinary differential equations (ODEs) system  system over a given time range.

# Arguments:
- `model::Any`: The ODEs model to be used for simulation. This can be a string identifier if it is an harcoded Kinbiont model or a function.
- `start_IC::Vector{Float64}`: The initial conditions for the ODE system.
- `tstart::Float64`: The starting time of the simulation.
- `tmax::Float64`: The final time of the simulation.
- `delta_t::Float64`: The time step size for simulation output.
- `param_of_ode`: A set of parameters required for the ODE model.

# Optional Arguments:
- `Integration_method=Tsit5()`: The SciML solver method used for integrating the ODE system. Other solvers, such as `KenCarp4(autodiff=true)`, may be used for specific applications.

# Output:
- Returns the simulated time series solution of the ODE system, it the standard SciML format.
"""
function ODEs_system_sim(
    model::Any, #string of the model
    start_IC::Vector{Float64}, # starting condition
    tstart::Float64, # start time of the sim
    tmax::Float64, # final time of the sim
    delta_t::Float64, # 
    param_of_ode; # parameters of the ODE model
    Integration_method=Tsit5(), # which sciml solver of ode
)

    t_steps = tstart:delta_t:tmax
    tspan = (tstart, tmax)
    ODE_prob = ODEs_system_model_selector2(model, start_IC, tspan; param=param_of_ode)


    sim = solve(ODE_prob, Integration_method, saveat=t_steps)

    return sim
end




"""
    fit_ODEs_System(
        data::Matrix{Float64},
        label_exp::String,
        model::Any,
        param, 
        Start_IC;
        set_of_equations_to_fit=nothing,
        Integration_method=Tsit5(),
        multistart=false, 
        n_restart=50, 
        optimizer=BBO_adaptive_de_rand_1_bin_radiuslimited(), 
        auto_diff_method=nothing,
        cons=nothing, 
        opt_params...
        )

Fits an ordinary differential equations (ODEs)  system to  data.

# Arguments:
- `data::Matrix{Float64}`: The dataset containing time points (first row) and observed values (second, third, fourth, ... rows).
- `label_exp::String`: The label for the experiment.
- `model::Any`: The ODEs system used for fitting. Can be a function or a model identifier.
- `param`: Initial guess for the model parameters, typically provided as a vector of `Float64`.
- `Start_IC`: Initial conditions for the ODEs system.

# Optional Arguments:
- `set_of_equations_to_fit=nothing`: Specifies a subset of equations to fit (if applicable). E.g., in a system of 4 equations you have data only on 1 and 4 then you should use `set_of_equations_to_fit = [1,4]`
- `Integration_method=Tsit5()`: The SciML solver used for numerical integration. For stiff models, consider `KenCarp4(autodiff=true)`.
- `type_of_smoothing="rolling_avg"`: Method for smoothing the data. Options: `"NO"`, `"rolling_avg"`, `"lowess"`.
- `thr_lowess=0.05`: Threshold parameter for lowess smoothing.
- `multistart=false`: Flag to enable or disable multistart optimization.
- `n_restart=50`: Number of optimization restarts, applicable if `multistart=true`.
- `optimizer=BBO_adaptive_de_rand_1_bin_radiuslimited()`: The optimization algorithm used for parameter estimation.
- `auto_diff_method=nothing`: Differentiation method for the optimizer, if required.
- `cons=nothing`: Constraints for the optimization process.
- `opt_params...`: Additional optional parameters for the optimizer (e.g., number of iterations).

# Output:
- Returns  a data struct containg the parameter of the best-fit ODEs system parameters along with the corresponding numerical solution.
"""
function fit_ODEs_System(data::Matrix{Float64}, # dataset first row times second row OD
    label_exp::String, #label of the experiment
    model::Any, # ode model to use
    param,
    Start_IC;
    set_of_equations_to_fit=nothing,
    Integration_method=Tsit5(), # selection of sciml Integration_method
    multistart=false,
    n_restart=50,
    optimizer=BBO_adaptive_de_rand_1_bin_radiuslimited(),
    auto_diff_method=nothing,
    cons=nothing,
    opt_params...)

    max_t = data[1, end]
    min_t = data[1, 1]
    tspan = (min_t, max_t)
    tsteps = data[1, :]



    # Defining ODE problem

    ODE_prob = ODEs_system_model_selector2(model, Start_IC, tspan)


    loss_function = define_loss_function_odes_system(data, set_of_equations_to_fit, ODE_prob, Integration_method, tsteps)

    res = KinbiontSolve(loss_function,
        Start_IC,
        param;
        opt=optimizer,
        auto_diff_method=auto_diff_method,
        multistart=multistart,
        n_restart=n_restart,
        cons=cons,
        opt_params...
    )

    remade_solution = solve(remake(ODE_prob, p=res.u), Integration_method, saveat=tsteps)

    loss_value = res.objective

    res_param = vectorize_df_results_ODE_System(label_exp, model, res.u, loss_value)

    Kinbiont_res_odes_system = ("ODEs_System", res_param, remade_solution)

    return Kinbiont_res_odes_system

end





export fit_ODEs_System
export ODEs_system_sim
export ODEs_system_model_selector2
