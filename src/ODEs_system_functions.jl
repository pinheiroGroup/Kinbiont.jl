
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



function define_loss_function_odes_system(data, set_of_equation_to_fit, ODE_prob, integrator, tsteps)


    if !isnothing(set_of_equation_to_fit)
        index_of_eqs = set_of_equation_to_fit 
        index_of_data = set_of_equation_to_fit .+1 


    else

        index_of_eqs = 1:1:(size(data)[1]-1)
        index_of_data = index_of_eqs .+ 1

    end



    function loss_RE_ODE_Sys(data, index_of_eqs, index_of_data, ODE_prob, integrator, p, tsteps)
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



        if size(sol_t)[2] == size(data)[2]

            lossa =
                0.5 * NaNMath.sum(abs2.(1.0 .- (data[index_of_data, :] ./ sol_t[index_of_eqs, :]))) /
                length(data[2, :])

        else

            lossa = 10.0^9 * length(data[2, :])
        end


        return lossa, sol
    end



    return (p) ->
        loss_RE_ODE_Sys(data, index_of_eqs, index_of_data, ODE_prob, integrator, p, tsteps)

end



function ODEs_system_sim(
    model::Any, #string of the model
    start_IC::Vector{Float64}, # starting condition
    tstart::Float64, # start time of the sim
    tmax::Float64, # final time of the sim
    delta_t::Float64, # 
    param_of_ode; # parameters of the ODE model
    integrator=KenCarp4(), # which sciml solver of ode
)

    t_steps = tstart:delta_t:tmax
    tspan = (tstart, tmax)
    ODE_prob = ODEs_system_model_selector2(model, start_IC, tspan; param=param_of_ode)


    sim = solve(ODE_prob, integrator, saveat=t_steps)

    return sim
end


function smoothing_data_ODE_System(
    data,
    set_of_equation_to_fit;
    method=type_of_smoothing,
    pt_avg=pt_avg,
    thr_lowess=thr_lowess
)

end

function fit_ODEs_System(data::Matrix{Float64}, # dataset first row times second row OD
    label_exp::String, #label of the experiment
    model::Any, # ode model to use
    param,
    Start_IC;
    set_of_equations_to_fit=nothing,
    integrator=Tsit5(), # selection of sciml integrator
    pt_avg=1, # numebr of the point to generate intial condition
    smoothing=false, # the smoothing is done or not?
    type_of_smoothing="rolling_avg",
    thr_lowess=0.05,
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

    # smoothing

    if smoothing == true

        data = smoothing_data_ODE_System(
            data,
            set_of_equations_to_fit;
            method=type_of_smoothing,
            pt_avg=pt_avg,
            thr_lowess=thr_lowess
        )

    end

    # Defining ODE problem

    ODE_prob = ODEs_system_model_selector2(model, Start_IC, tspan)


    loss_function = define_loss_function_odes_system(data, set_of_equations_to_fit, ODE_prob, integrator, tsteps)

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

    remade_solution = solve(remake(ODE_prob, p=res.u), integrator, saveat=tsteps)

    loss_value = res.objective

    res_param = vectorize_df_results_ODE_System(label_exp, model, res.u, loss_value)

    Kinbiont_res_odes_system = ("ODEs_System", res_param, remade_solution)

    return Kinbiont_res_odes_system

end





export fit_ODEs_System
export ODEs_system_sim
export ODEs_system_model_selector