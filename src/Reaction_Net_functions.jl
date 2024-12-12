using Catalyst







function Reaction_network_model_selector(model, u0, tspan; param=nothing)

    if typeof(model) == String

        model_string =  Kinbiont_Reaction_network_models[model].name
        model_function =  Kinbiont_Reaction_network_models[model].network

        ODE_prob = ODEProblem(model_function, u0, tspan, param)


    else

        model_string = "custom"

        ODE_prob = ODEProblem(model, u0, tspan, param)

    end


    return ODE_prob
end








function Kinbiont_Reaction_network_sim(model,
    start_IC,
    tstart::Float64, # start time of the sim
    tmax::Float64, # final time of the sim
    delta_t::Float64, # 
    param; # parameters of the ODE model
    integrator=KenCarp4(), # which sciml solver of ode
    )



    t_steps = tstart:delta_t:tmax
    tspan = (tstart, tmax)
    ODE_prob = Reaction_network_model_selector(model, start_IC, tspan; param=param)




    # Simulate ODE and plot results.

    sim = solve(ODE_prob, integrator, saveat=t_steps)
   return sim
end







