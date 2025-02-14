using Catalyst
using DiffEqParamEstim
prob_generator_t = (prob,fixed_param) -> remake(prob, u0 = u0, p = replace_in_array2(fixed_param))


STANDARD_PROB_GENERATOR(prob, p) = remake(prob; u0 = eltype(p).(prob.u0), p = replace_in_array2(p))

     
function replace_in_array2(arr)

    param_counter = 0
          return  [isnothing(val) ? p[(param_counter += 1)] : val for val in arr]
end




struct Kinbiont_Reaction_network_res
    network::Any
    param::Any
    fit_times::Any
    fit::Any
    loss::Any
end

function check_bounds_opt_rn(loss_function,
    opt,
    p_guess,
    lb,
    ub)
 
    if SciMLBase.requiresbounds(opt)





        if  isnothing(lb) || isnothing(ub)

            @warn "The used optimization method requires box bounds, Kinbiont.jl will use upper bounds that are 10 times the guess
             and lower bounds that are 10 times lower the guess.
             This choice can be suboptimal. Note that the Kinbiont.jl default optimizer requires box bounds to guaranteed a Real N(t) and positive parameters.
             To avoid to specify the bounds use you can use an optimizer that do not require it, e.g., `optimizer = NOMADOpt()`. 
             Note that numerical instabilities may occur.
             "
             optprob = OptimizationProblem(loss_function,p_guess;lb=p_guess./10,ub=p_guess.*10)

        else
             optprob = OptimizationProblem(loss_function,p_guess;lb=lb,ub=ub)

        end
    
    else

        optprob = OptimizationProblem(loss_function,p_guess)

    end
    
    return optprob

end

"""
    RN_fit(data,
    rn,
    u0,
    param_g;
    Integration_method=Tsit5(),
    optimizer=BBO_adaptive_de_rand_1_bin_radiuslimited(),
    adtype = SciMLBase.NoAD(),
    lb = nothing,
    ub = nothing,
    set_of_equation_to_fit= nothing,
    opt_params...)

Fits a reaction network model to experimental data by optimizing the parameters of the system.

# Arguments:
- `data`: A matrix containing experimental data, where the first row represents time points and subsequent rows represent observed variables.
- `rn`: A reaction network model,  defined using Catalyst.jl, or a string if is an harcoded Kinbiont model.
- `u0`: A vector of initial conditions for the system, following the variable order in the reaction network.
- `param_g`: A vector of initial guesses for the parameters to be estimated. **The order of parameters must match the order in which they are declared in the reaction network model (`rn`).**

# Optional Arguments:
- `Integration_method=Tsit5()`: The numerical solver used for integrating the reaction network.
- `optimizer=BBO_adaptive_de_rand_1_bin_radiuslimited()`: The optimization algorithm for parameter estimation.
- `adtype = SciMLBase.NoAD()`: The automatic differentiation method for gradient computation.
- `lb = nothing`: Lower bounds for the parameters, if applicable.
- `ub = nothing`: Upper bounds for the parameters, if applicable.
- `set_of_equation_to_fit= nothing`: Specifies a subset of equations to be fitted, if only a part of the system needs optimization. **The order of species must match the order in which they are declared in the reaction network model (`rn`).**
- `opt_params...`: Additional parameters for the optimization routine.

# Output:
- Optimized parameter values that best fit the experimental data.
- The fitted model.
"""

function RN_fit(data, 
    rn,
    u0,
    param_g;
    Integration_method=Tsit5(), 
    optimizer=BBO_adaptive_de_rand_1_bin_radiuslimited(),
    adtype = SciMLBase.NoAD(),
    lb = nothing,
    ub = nothing,
    fixed_parameters = nothing,
    set_of_equation_to_fit= nothing,
    opt_params...)

    rn = Reaction_network_model_selector2(rn)

    max_t = data[1, end]
    min_t = data[1, 1]
    tspan = (min_t, max_t)
    tsteps = data[1, :]
    # generate numeric guess (NOTE THAT ORDER IS FIXED AND THE USER SHOULD PAY ATTENTION)
    param_guess = [param_g[i][2] for i in 1:length(param_g)]
    ODE_prob = ODEProblem(rn, u0, tspan, param_g)
    
    # fixing parameters

    data_to_fit =  data[2:end,:]
    # change input lb ub 

    if isnothing(set_of_equation_to_fit) && isnothing(fixed_parameters)

        loss_function = build_loss_objective(ODE_prob, Integration_method, L2Loss(tsteps, data_to_fit), adtype;
                                     opt_params...)
    elseif !isnothing(set_of_equation_to_fit) && isnothing(fixed_parameters)

        loss_function = build_loss_objective(ODE_prob, Integration_method, L2Loss(tsteps, data_to_fit), adtype;
        save_idxs = set_of_equation_to_fit,opt_params...)

    elseif isnothing(set_of_equation_to_fit) && !isnothing(fixed_parameters)

        loss_function = build_loss_objective(ODE_prob, Integration_method, L2Loss(tsteps, data_to_fit), adtype;
        prob_generator = STANDARD_PROB_GENERATOR, opt_params...)


    elseif !isnothing(set_of_equation_to_fit) && !isnothing(fixed_parameters)

        loss_function = build_loss_objective(ODE_prob, Integration_method, L2Loss(tsteps, data_to_fit), adtype;
        save_idxs = set_of_equation_to_fit,prob_generator = STANDARD_PROB_GENERATOR, opt_params...)

    end


    # check bounds of optimizer
    # Set bounds for optimization if needed
    optprob = check_bounds_opt_rn(loss_function,optimizer, param_guess,lb,ub)





    opt_sol = solve(optprob, optimizer)
    res_parameters = opt_sol.u
    loss = opt_sol.objective

    ODE_prob_fitted = remake(ODE_prob; p = opt_sol.u)
    fitted_sol = solve(ODE_prob_fitted, Integration_method)

    Kinbiont_Reaction_network_res = (rn,res_parameters,fitted_sol.t,fitted_sol,loss)


    return Kinbiont_Reaction_network_res
end


function Reaction_network_model_selector(model, u0, tspan; param=nothing, t_steps =nothing)

    if typeof(model) == String

        model_string =  Kinbiont_Reaction_network_models[model].name
        model_function =  Kinbiont_Reaction_network_models[model].network

        if isnothing(t_steps)
             ODE_prob = ODEProblem(model_function, u0, tspan, param)
        else
            ODE_prob = ODEProblem(model_function, u0, tspan, param;saveat= t_steps)

        end

    else

        model_string = "custom"

        if isnothing(t_steps)
            ODE_prob = ODEProblem(model, u0, tspan, param)
       else
           ODE_prob = ODEProblem(model, u0, tspan, param;saveat= t_steps)

       end
    end


    return ODE_prob
end




function Reaction_network_model_selector2(rn)

    if typeof(rn) == String

        model_function =  Kinbiont_Reaction_network_models[rn].network

    else
        model_function = rn
    end


    return model_function
end

"""
    Kinbiont_Reaction_network_sim(model,
    start_IC,
    tstart::Float64,
    tmax::Float64,
    delta_t::Float64,
    param;
    Integration_method=Tsit5())

Simulates a reaction network using the Catalyst.jl package over a specified time range.

# Arguments:
- `model`: Either a string specifying the reaction network model or a `ReactionSystem` from Catalyst.jl.
- `start_IC`: A vector containing the initial conditions for the system (biomass, substrates, and proteins).
- `tstart::Float64`: The start time of the simulation.
- `tmax::Float64`: The final time of the simulation.
- `delta_t::Float64`: The time step for numerical integration.
- `param`: A set of parameters used in the ODE model.

# Optional Arguments:
- `Integration_method=Tsit5()`: The numerical solver from the SciML ecosystem. The default solver is `Tsit5()`, but other solvers (e.g., `KenCarp4(autodiff=true)`) can be used.

# Output:
- A data structure containing:
  1. The time points of the simulation.
  2. The numerical solution of the ODE system for biomass, substrates, and proteins.
"""
function Kinbiont_Reaction_network_sim(model,
    start_IC,
    tstart::Float64, # start time of the sim
    tmax::Float64, # final time of the sim
    delta_t::Float64, # 
    param; # parameters of the ODE model
    Integration_method=Tsit5(), # which sciml solver of ode
    )



    t_steps = tstart:delta_t:tmax
    tspan = (tstart, tmax)
    ODE_prob = Reaction_network_model_selector(model, start_IC, tspan; param=param,t_steps=t_steps)




    # Simulate ODE and plot results.

    sim = solve(ODE_prob, Integration_method, saveat=t_steps)
   return sim
end


export RN_fit
export Kinbiont_Reaction_network_sim
