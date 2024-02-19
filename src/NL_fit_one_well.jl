



function fit_NL_model(data::Matrix{Float64}, # dataset first row times second row OD
    name_well::String, # name of the well
    label_exp::String, #label of the experiment
    model_function::Any, # ode model to use
    lb_param::Vector{Float64}, # lower bound param
    ub_param::Vector{Float64}; # upper bound param
    u0=lb_param .+ (ub_param .- lb_param) ./ 2,# initial guess param
    optmizator=BBO_adaptive_de_rand_1_bin_radiuslimited(),
    display_plots=true, # display plots in julia or not
    save_plot=false,
    path_to_plot="NA", # where save plots
    pt_avg=1, # numebr of the point to generate intial condition
    pt_smooth_derivative=7,
    smoothing=false, # the smoothing is done or not?
    type_of_smoothing="rolling_avg",
    type_of_loss="RE", # type of used loss
    multiple_scattering_correction=false, # if true uses the given calibration curve to fix the data
    method_multiple_scattering_correction="interpolation",
    calibration_OD_curve="NA",  #  the path to calibration curve to fix the data
    PopulationSize=300,
    maxiters=2000000,
    abstol=0.00001,
    thr_lowess=0.05,
    penality_CI=3.0,
)

    if multiple_scattering_correction == true

        data = correction_OD_multiple_scattering(data, calibration_OD_curve; method=method_multiple_scattering_correction)

    end

    if smoothing == true
        data = smoothing_data(
            data;
            method=type_of_smoothing,
            pt_avg=pt_avg,
            thr_lowess=thr_lowess
        )
    end
    # setting initial conditions
    # TO DO GIVE THE OPTION TO FIX THEM AT LEAST IN KNOWN MODELS
    # TO DO MODEL SELECTOR
    if typeof(model_function) == String

        model_string = NL_models[model_function].name
        model_function = NL_models[model_string].func


    else

        model_string = "custom"


    end

    loss_function =
        select_loss_function_NL(type_of_loss, data, penality_CI, model_function)

    prob = OptimizationProblem(loss_function, u0, data, lb=lb_param, ub=ub_param)

    # Solve the optimization problem
    sol = solve(prob, optmizator, PopulationSize=PopulationSize, maxiters=maxiters, abstol=abstol)
    # evaluate the fitted  model
    fitted_model = model_function(sol, data[1, :])

    if display_plots
        if_display = display
    else
        if_display = identity
    end

    if save_plot
        mkpath(path_to_plot)
    end


    # plotting if required
    if_display(
        Plots.scatter(
            data[1, :],
            data[2, :],
            xlabel="Time",
            ylabel="Arb. Units",
            label=["Data " nothing],
            markersize=2,
            color=:black,
            title=string(label_exp, " ", name_well),
        ),)

    if_display(
        Plots.plot!(
            data[1, :],
            fitted_model,
            xlabel="Time",
            ylabel="Arb. Units",
            label=[string("Fitting ", model_string) nothing],
            c=:red,
        ),
    )
    if save_plot
        png(string(path_to_plot, label_exp, "_", model, "_", name_well, ".png"))
    end

    sol_fin, index_not_zero = remove_negative_value(fitted_model)

    data_th = transpose(hcat(data[1, index_not_zero], sol_fin))


    max_th_gr = maximum(specific_gr_evaluation(Matrix(data_th), pt_smooth_derivative))

    # max empirical gr
    max_em_gr = maximum(specific_gr_evaluation(data, pt_smooth_derivative))
    loss_value = sol.objective


    res_param = [[name_well, model_string], [sol[1:end]], [max_th_gr, max_em_gr, loss_value]]


    res_param = reduce(vcat, reduce(vcat, res_param))


    return res_param, fitted_model
end






function fit_NL_model_with_sensitivity(data::Matrix{Float64}, # dataset first row times second row OD
    name_well::String, # name of the well
    label_exp::String, #label of the experiment
    model_function::Any, # ode model to use
    lb_param::Vector{Float64}, # lower bound param
    ub_param::Vector{Float64}; # upper bound param
    nrep=100,
    optmizator=BBO_adaptive_de_rand_1_bin_radiuslimited(),
    display_plots=true, # display plots in julia or not
    save_plot=false,
    path_to_plot="NA", # where save plots
    pt_avg=1, # numebr of the point to generate intial condition
    pt_smooth_derivative=7,
    smoothing=false, # the smoothing is done or not?
    type_of_smoothing="rolling_avg",
    type_of_loss="RE", # type of used loss
    multiple_scattering_correction=false, # if true uses the given calibration curve to fix the data
    method_multiple_scattering_correction="interpolation",
    calibration_OD_curve="NA",  #  the path to calibration curve to fix the data
    PopulationSize=300,
    maxiters=2000000,
    abstol=0.00001,
    thr_lowess=0.05,
    write_res=false,
    penality_CI=3.0,
)


    if multiple_scattering_correction == true

        data = correction_OD_multiple_scattering(data, calibration_OD_curve; method=method_multiple_scattering_correction)

    end

    if smoothing == true
        data = smoothing_data(
            data;
            method=type_of_smoothing,
            pt_avg=pt_avg,
            thr_lowess=thr_lowess
        )
    end
    # setting initial conditions
    # TO DO GIVE THE OPTION TO FIX THEM AT LEAST IN KNOWN MODELS
    # TO DO MODEL SELECTOR
    if typeof(model_function) == String
        model_string = NL_models[model_function].name
        model_function = NL_models[model_string].func

    else

        model_string = "custom"


    end
    # Define the optimization problem LOSS


    loss_function =
        select_loss_function_NL(type_of_loss, data, penality_CI, model_function)


    max_em_gr = maximum(specific_gr_evaluation(data, pt_smooth_derivative))

    fin_param = initialize_df_results_ode_custom(lb_param)
    param_combination =
        generation_of_combination_of_IC_morris(lb_param, ub_param, nrep)

    for i = 1:size(param_combination)[2]
        u0 = param_combination[:, i]



        prob = OptimizationProblem(loss_function, u0, data, lb=lb_param, ub=ub_param)

        # Solve the optimization problem
        sol = solve(prob, optmizator, PopulationSize=PopulationSize, maxiters=maxiters, abstol=abstol)
        # evaluate the fitted  model
        fitted_model = model_function(sol, data[1, :])

        sol_fin, index_not_zero = remove_negative_value(fitted_model)

        data_th = transpose(hcat(data[1, index_not_zero], sol_fin))
        if size(index_not_zero, 1) > pt_smooth_derivative + 3
            max_th_gr = maximum(specific_gr_evaluation(Matrix(data_th), pt_smooth_derivative))
        else
            max_th_gr = missing
        end

        loss_value = sol.objective


        res_param = [[name_well, model_string], [sol[1:end]], [max_th_gr, max_em_gr, loss_value]]
        res_param = reduce(vcat, reduce(vcat, res_param))



        fin_param = hcat(fin_param, res_param)


    end


    index_best = findmin(fin_param[end, 2:end])[2]


    best_res_param = fin_param[:, index_best+1]

    best_fitted_model = model_function(best_res_param[3:(end-3)], data[1, :])

    if display_plots
        if_display = display
    else
        if_display = identity
    end

    if save_plot
        mkpath(path_to_plot)
    end


    # plotting if required
    if_display(
        Plots.scatter(
            data[1, :],
            data[2, :],
            xlabel="Time",
            ylabel="Arb. Units",
            label=["Data " nothing],
            markersize=2,
            color=:black,
            title=string(label_exp, " ", name_well),
        ),)

    if_display(
        Plots.plot!(
            data[1, :],
            best_fitted_model,
            xlabel="Time",
            ylabel="Arb. Units",
            label=[string("Fitting ", model_string) nothing],
            c=:red,
        ),
    )
    if save_plot
        png(string(path_to_plot, label_exp, "_", model_string, "_", name_well, ".png"))
    end


    if write_res == true
        mkpath(path_to_results)
        CSV.write(
            string(path_to_results, label_exp, "_", model_string, "_results_sensitivity.csv"),
            Tables.table(Matrix(fin_param)),
        )
        CSV.write(
            string(path_to_results, label_exp, "_", model_string, "_configurations_tested.csv"),
            Tables.table(Matrix(param_combination)),
        )
    end


    return best_res_param, best_fitted_model, fin_param, Matrix(param_combination)
end


function fit_NL_model_MCMC_intialization(data::Matrix{Float64}, # dataset first row times second row OD
    name_well::String, # name of the well
    label_exp::String, #label of the experiment
    model_function::Any, # ode model to use
    lb_param::Vector{Float64}, # lower bound param
    ub_param::Vector{Float64}; # upper bound param
    u0=lb_param .+ (ub_param .- lb_param) ./ 2,# initial guess param
    nrep=100,
    optmizator=BBO_adaptive_de_rand_1_bin_radiuslimited(),
    display_plots=true, # display plots in julia or not
    save_plot=false,
    path_to_plot="NA", # where save plots
    pt_avg=1, # numebr of the point to generate intial condition
    pt_smooth_derivative=7,
    smoothing=false, # the smoothing is done or not?
    type_of_smoothing="rolling_avg",
    type_of_loss="RE", # type of used loss
    multiple_scattering_correction=false, # if true uses the given calibration curve to fix the data
    method_multiple_scattering_correction="interpolation",
    calibration_OD_curve="NA",  #  the path to calibration curve to fix the data
    PopulationSize=300,
    maxiters=2000000,
    abstol=0.00001,
    thr_lowess=0.05,
    penality_CI=3.0,
)

    u0_best = copy(u0)
    loss_best = 10^20
    loss_chain = copy(loss_best)
    loss_chain_best = copy(loss_best)
    best_fitted_model = Any
    best_res_param = Any
    if multiple_scattering_correction == true

        data = correction_OD_multiple_scattering(data, calibration_OD_curve; method=method_multiple_scattering_correction)

    end

    if smoothing == true
        data = smoothing_data(
            data;
            method=type_of_smoothing,
            pt_avg=pt_avg,
            thr_lowess=thr_lowess
        )
    end
    # setting initial conditions
    # TO DO GIVE THE OPTION TO FIX THEM AT LEAST IN KNOWN MODELS
    # TO DO MODEL SELECTOR
    if typeof(model_function) == String
        model_string = NL_models[model_function].name
        model_function = NL_models[model_string].func

    else

        model_string = "custom"


    end
    # Define the optimization problem LOSS


    loss_function =
        select_loss_function_NL(type_of_loss, data, penality_CI, model_function)


    max_em_gr = maximum(specific_gr_evaluation(data, pt_smooth_derivative))


    for i = 1:nrep


        index_to_change = rand(1:length(lb_param), 1)[1]
        new_param = rand(Uniform(lb_param[index_to_change], ub_param[index_to_change]), 1)[1]
        u0 = copy(u0_best)
        u0[index_to_change] = new_param
        prob = OptimizationProblem(loss_function, u0, data, lb=lb_param, ub=ub_param)

        # Solve the optimization problem
        sol = solve(prob, optmizator, PopulationSize=PopulationSize, maxiters=maxiters, abstol=abstol)
        # evaluate the fitted  model
        fitted_model = model_function(sol, data[1, :])

        sol_fin, index_not_zero = remove_negative_value(fitted_model)

        data_th = transpose(hcat(data[1, index_not_zero], sol_fin))
        if size(index_not_zero, 1) > pt_smooth_derivative + 3
            max_th_gr = maximum(specific_gr_evaluation(Matrix(data_th), pt_smooth_derivative))
        else
            max_th_gr = missing
        end

        loss_value = sol.objective

        loss_chain = vcat(loss_chain, loss_value)
        res_param = [[name_well, model_string], [sol[1:end]], [max_th_gr, max_em_gr, loss_value]]
        res_param = reduce(vcat, reduce(vcat, res_param))


        if loss_value < loss_best
            loss_best = copy(loss_value)

            best_res_param = copy(res_param)

            best_fitted_model = model_function(best_res_param[3:(end-3)], data[1, :])
            u0_best = copy(u0)
            loss_chain = vcat(loss_chain, loss_value)


        else


            prob = 1 - abs(loss_value - loss_best) / loss_best

            nrand = rand(Uniform(0.0, 1.0), 1)[1]

            if nrand < prob

                loss_best = copy(loss_value)

                best_res_param = copy(res_param)

                best_fitted_model = model_function(best_res_param[3:(end-3)], data[1, :])

                u0_best = copy(u0)

            end

        end
        loss_chain_best = vcat(loss_chain, loss_best)



    end


    if display_plots
        if_display = display
    else
        if_display = identity
    end

    if save_plot
        mkpath(path_to_plot)
    end


    # plotting if required
    if_display(
        Plots.scatter(
            data[1, :],
            data[2, :],
            xlabel="Time",
            ylabel="Arb. Units",
            label=["Data " nothing],
            markersize=2,
            color=:black,
            title=string(label_exp, " ", name_well),
        ),)

    if_display(
        Plots.plot!(
            data[1, :],
            best_fitted_model,
            xlabel="Time",
            ylabel="Arb. Units",
            label=[string("Fitting ", model_string) nothing],
            c=:red,
        ),
    )
    if save_plot
        png(string(path_to_plot, label_exp, "_", model_string, "_", name_well, ".png"))
    end





    return best_res_param, best_fitted_model, loss_chain[2:end], loss_chain_best[2:end]
end






function fit_NL_model_bootstrap(data::Matrix{Float64}, # dataset first row times second row OD
    name_well::String, # name of the well
    label_exp::String, #label of the experiment
    model_function::Any, # ode model to use
    lb_param::Vector{Float64}, # lower bound param
    ub_param::Vector{Float64}; # upper bound param
    nrep=100,
    u0=lb_param .+ (ub_param .- lb_param) ./ 2,# initial guess param
    optmizator=BBO_adaptive_de_rand_1_bin_radiuslimited(),
    display_plots=true, # display plots in julia or not
    save_plot=false,
    size_bootstrap=0.7,
    path_to_plot="NA", # where save plots
    pt_avg=1, # numebr of the point to generate intial condition
    pt_smooth_derivative=7,
    smoothing=false, # the smoothing is done or not?
    type_of_smoothing="rolling_avg",
    type_of_loss="RE", # type of used loss
    multiple_scattering_correction=false, # if true uses the given calibration curve to fix the data
    method_multiple_scattering_correction="interpolation",
    calibration_OD_curve="NA",  #  the path to calibration curve to fix the data
    PopulationSize=300,
    maxiters=2000000,
    abstol=0.00001,
    thr_lowess=0.05,
    write_res=false,
    penality_CI=3.0
)


    if multiple_scattering_correction == true

        data = correction_OD_multiple_scattering(data, calibration_OD_curve; method=method_multiple_scattering_correction)

    end

    if smoothing == true
        data = smoothing_data(
            data;
            method=type_of_smoothing,
            pt_avg=pt_avg,
            thr_lowess=thr_lowess
        )
    end
    # setting initial conditions
    # TO DO GIVE THE OPTION TO FIX THEM AT LEAST IN KNOWN MODELS
    # TO DO MODEL SELECTOR
    if typeof(model_function) == String
        model_string = NL_models[model_function].name
        model_function = NL_models[model_string].func

    else

        model_string = "custom"


    end
    # Define the optimization problem LOSS





    max_em_gr = maximum(specific_gr_evaluation(data, pt_smooth_derivative))

    fin_param = initialize_df_results_ode_custom(lb_param)


    for i = 1:nrep
        idxs = rand(1:1:size(data[2, :], 1), convert.(Int, floor(length(data[2, :]) * size_bootstrap)))
        times_boot = data[1, idxs]
        data_boot = data[2, idxs]

        idxs_2 = sortperm(times_boot)

        times_boot = times_boot[idxs_2]
        data_boot = data_boot[idxs_2]

        data_to_fit = Matrix(transpose(hcat(times_boot, data_boot)))


        loss_function =
            select_loss_function_NL(type_of_loss, data, penality_CI, model_function)

        prob = OptimizationProblem(loss_function, u0, data_to_fit, lb=lb_param, ub=ub_param)

        # Solve the optimization problem
        sol = solve(prob, optmizator, PopulationSize=PopulationSize, maxiters=maxiters, abstol=abstol)
        # evaluate the fitted  model
        fitted_model = model_function(sol, data_to_fit[1, :])
        sol_fin, index_not_zero = remove_negative_value(fitted_model)

        data_th = transpose(hcat(data_to_fit[1, index_not_zero], sol_fin))


        if size(index_not_zero, 1) > pt_smooth_derivative + 3
            max_th_gr = maximum(specific_gr_evaluation(Matrix(data_th), pt_smooth_derivative))
        else
            max_th_gr = missing
        end

        loss_value = sol.objective


        res_param = [[name_well, model_string], [sol[1:end]], [max_th_gr, max_em_gr, loss_value]]
        res_param = reduce(vcat, reduce(vcat, res_param))



        fin_param = hcat(fin_param, res_param)


    end


    index_best = findmin(fin_param[end, 2:end])[2]

    best_res_param = fin_param[:, index_best+1]
    best_fitted_model = model_function(best_res_param[3:(end-3)], data[1, :])

    if display_plots
        if_display = display
    else
        if_display = identity
    end

    if save_plot
        mkpath(path_to_plot)
    end


    # plotting if required
    if_display(
        Plots.scatter(
            data[1, :],
            data[2, :],
            xlabel="Time",
            ylabel="Arb. Units",
            label=["Data " nothing],
            markersize=2,
            color=:black,
            title=string(label_exp, " ", name_well),
        ),)

    if_display(
        Plots.plot!(
            data[1, :],
            best_fitted_model,
            xlabel="Time",
            ylabel="Arb. Units",
            label=[string("Fitting ", model_string) nothing],
            c=:red,
        ),
    )
    if save_plot
        png(string(path_to_plot, label_exp, "_", model_string, "_", name_well, ".png"))
    end


    if write_res == true
        mkpath(path_to_results)
        CSV.write(
            string(path_to_results, label_exp, "_", model_string, "_results_bootstrap.csv"),
            Tables.table(Matrix(fin_param)),
        )

    end
    mean_param = [mean(fin_param[i, 2:end]) for i in 3:axes(fin_param)[1][end]]
    sd_param = [std(fin_param[i, 2:end]) for i in 3:axes(fin_param)[1][end]]

    return best_res_param, best_fitted_model, fin_param, mean_param, sd_param
end

















function NL_model_selection(data::Matrix{Float64}, # dataset first row times second row OD
    name_well::String, # name of the well
    label_exp::String, #label of the experiment
    list_model_function::Any, # ode model to use
    list_lb_param::Any, # lower bound param
    list_ub_param::Any; # upper bound param
    method_of_fitting="MCMC",
    nrep=100,
    list_u0=list_lb_param .+ (list_ub_param .- list_lb_param) ./ 2,# initial guess param
    optmizator=BBO_adaptive_de_rand_1_bin_radiuslimited(),
    display_plots=true, # display plots in julia or not
    save_plot=false,
    size_bootstrap=0.7,
    path_to_plot="NA", # where save plots
    pt_avg=1, # numebr of the point to generate intial condition
    pt_smooth_derivative=7,
    smoothing=false, # the smoothing is done or not?
    type_of_smoothing="rolling_avg",
    type_of_loss="RE", # type of used loss
    multiple_scattering_correction=false, # if true uses the given calibration curve to fix the data
    method_multiple_scattering_correction="interpolation",
    calibration_OD_curve="NA",  #  the path to calibration curve to fix the data
    PopulationSize=300,
    maxiters=2000000,
    abstol=0.00001,
    thr_lowess=0.05,
    write_res=false,
    beta_param=2.0,
    penality_CI=8.0
)
    score_res = ["AIC"]
    top_score = 10^20
    top_model = Vector{Any}
    top_fitted_sol = Vector{Any}

    for mm in 1:eachindex(list_model_function)[end]

        model_to_test = list_model_function[mm]
        lb_param = list_lb_param[mm]
        ub_param = list_ub_param[mm]
        u0 = list_u0[mm]

        if method_of_fitting == "Bootstrap"



            temp_res = fit_NL_model_bootstrap(data, # dataset first row times second row OD
                name_well, # name of the well
                label_exp, #label of the experiment
                model_to_test, # ode model to use
                lb_param, # lower bound param
                ub_param; # upper bound param
                nrep=nrep,
                u0=u0,# initial guess param
                optmizator=optmizator,
                display_plots=display_plots, # display plots in julia or not
                save_plot=save_plot,
                size_bootstrap=size_bootstrap,
                path_to_plot=path_to_plot, # where save plots
                pt_avg=pt_avg, # numebr of the point to generate intial condition
                pt_smooth_derivative=pt_smooth_derivative,
                smoothing=smoothing, # the smoothing is done or not?
                type_of_smoothing=type_of_smoothing,
                type_of_loss=type_of_loss, # type of used loss
                multiple_scattering_correction=multiple_scattering_correction, # if true uses the given calibration curve to fix the data
                method_multiple_scattering_correction=method_multiple_scattering_correction,
                calibration_OD_curve=calibration_OD_curve,  #  the path to calibration curve to fix the data
                PopulationSize=PopulationSize,
                maxiters=maxiters,
                abstol=abstol,
                thr_lowess=thr_lowess,
                write_res=write_res,
                penality_CI=penality_CI
            )

            n_param = length(lb_param)

            temp_AIC = AICc_evaluation(n_param, beta_param, data[2, :], temp_res[2])
            temp = [model_to_test, temp_res[end], temp_AIC]
            score_res = hcat(score_res, temp_AIC)
            if mm == 1

                top_score = copy(temp_AIC)
                top_model = copy(temp_res[1])
                top_fitted_sol = copy(temp_res[2])


            elseif top_score > temp_AIC
                top_score = copy(temp_AIC)
                top_model = copy(temp_res[1])
                top_fitted_sol = copy(temp_res[2])

            end


        elseif method_of_fitting == "Morris_sensitivity"


            temp_res = fit_NL_model_with_sensitivity(data, # dataset first row times second row OD
                name_well, # name of the well
                label_exp, #label of the experiment
                model_to_test, # ode model to use
                lb_param, # lower bound param
                ub_param; # upper bound param
                nrep=nrep,
                optmizator=optmizator,
                display_plots=display_plots, # display plots in julia or not
                save_plot=save_plot,
                path_to_plot=path_to_plot, # where save plots
                pt_avg=pt_avg, # numebr of the point to generate intial condition
                pt_smooth_derivative=pt_smooth_derivative,
                smoothing=smoothing, # the smoothing is done or not?
                type_of_smoothing=type_of_smoothing,
                type_of_loss=type_of_loss, # type of used loss
                multiple_scattering_correction=multiple_scattering_correction, # if true uses the given calibration curve to fix the data
                method_multiple_scattering_correction=method_multiple_scattering_correction,
                calibration_OD_curve=calibration_OD_curve,  #  the path to calibration curve to fix the data
                PopulationSize=PopulationSize,
                maxiters=maxiters,
                abstol=abstol,
                thr_lowess=thr_lowess,
                write_res=write_res,
                penality_CI=penality_CI
            )

            n_param = length(lb_param)

            temp_AIC = AICc_evaluation(n_param, beta_param, data[2, :], temp_res[2])
            temp = [model_to_test, temp_res[end], temp_AIC]

            score_res = hcat(score_res, temp_AIC)

            if mm == 1

                top_score = copy(temp_AIC)
                top_model = copy(temp_res[1])
                top_fitted_sol = copy(temp_res[2])


            elseif top_score > temp_AIC
                top_score = copy(temp_AIC)
                top_model = copy(temp_res[1])
                top_fitted_sol = copy(temp_res[2])

            end

        elseif method_of_fitting == "MCMC"


            temp_res = fit_NL_model_MCMC_intialization(data, # dataset first row times second row OD
                name_well, # name of the well
                label_exp, #label of the experiment
                model_to_test, # ode model to use
                lb_param, # lower bound param
                ub_param; # upper bound param
                nrep=nrep,
                optmizator=optmizator,
                display_plots=display_plots, # display plots in julia or not
                save_plot=save_plot,
                path_to_plot=path_to_plot, # where save plots
                pt_avg=pt_avg, # numebr of the point to generate intial condition
                pt_smooth_derivative=pt_smooth_derivative,
                smoothing=smoothing, # the smoothing is done or not?
                type_of_smoothing=type_of_smoothing,
                type_of_loss=type_of_loss, # type of used loss
                multiple_scattering_correction=multiple_scattering_correction, # if true uses the given calibration curve to fix the data
                method_multiple_scattering_correction=method_multiple_scattering_correction,
                calibration_OD_curve=calibration_OD_curve,  #  the path to calibration curve to fix the data
                PopulationSize=PopulationSize,
                maxiters=maxiters,
                abstol=abstol,
                thr_lowess=thr_lowess,
                penality_CI=penality_CI)

            n_param = length(lb_param)

            temp_AIC = AICc_evaluation(n_param, beta_param, data[2, :], temp_res[2])
            temp = [model_to_test, temp_res[end], temp_AIC]

            score_res = hcat(score_res, temp_AIC)
            if mm == 1

                top_score = copy(temp_AIC)
                top_model = copy(temp_res[1])
                top_fitted_sol = copy(temp_res[2])


            elseif top_score > temp_AIC
                top_score = copy(temp_AIC)
                top_model = copy(temp_res[1])
                top_fitted_sol = copy(temp_res[2])

            end


        else



            temp_res = fit_NL_model(data, # dataset first row times second row OD
                name_well, # name of the well
                label_exp, #label of the experiment
                model_to_test, # ode model to use
                lb_param, # lower bound param
                ub_param; # upper bound param
                u0=u0,# initial guess param
                optmizator=optmizator,
                display_plots=display_plots, # display plots in julia or not
                save_plot=save_plot,
                path_to_plot=path_to_plot, # where save plots
                pt_avg=pt_avg, # numebr of the point to generate intial condition
                pt_smooth_derivative=pt_smooth_derivative,
                smoothing=smoothing, # the smoothing is done or not?
                type_of_smoothing=type_of_smoothing,
                type_of_loss=type_of_loss, # type of used loss
                multiple_scattering_correction=multiple_scattering_correction, # if true uses the given calibration curve to fix the data
                method_multiple_scattering_correction=method_multiple_scattering_correction,
                calibration_OD_curve=calibration_OD_curve,  #  the path to calibration curve to fix the data
                PopulationSize=PopulationSize,
                maxiters=maxiters,
                abstol=abstol,
                thr_lowess=thr_lowess,
                penality_CI=penality_CI
            )


            n_param = length(lb_param)

            temp_AIC = AICc_evaluation(n_param, beta_param, data[2, :], temp_res[2])
            temp = [model_to_test, temp_res[end], temp_AIC]

            score_res = hcat(score_res, temp_AIC)
            if mm == 1

                top_score = copy(temp_AIC)
                top_model = copy(temp_res[1])
                top_fitted_sol = copy(temp_res[2])


            elseif top_score > temp_AIC
                top_score = copy(temp_AIC)
                top_model = copy(temp_res[1])
                top_fitted_sol = copy(temp_res[2])

            end



        end


    end

    # plotting if required
    if display_plots
        if_display = display
    else
        if_display = identity
    end

    if save_plot
        mkpath(path_to_plot)
    end
    if_display(
        Plots.scatter(
            data[1, :],
            data[2, :],
            xlabel="Time",
            ylabel="Arb. Units",
            label=["Data " nothing],
            markersize=2,
            color=:black,
            title=string(label_exp, " ", name_well),
        ),)

    if_display(
        Plots.plot!(
            data[1, :],
            top_fitted_sol,
            xlabel="Time",
            ylabel="Arb. Units",
            label=[string("Fitting ", top_model[2]) nothing],
            c=:red,
        ),
    )
    if save_plot
        png(string(path_to_plot, label_exp, "_", model, "_", name_well, ".png"))
    end


    return top_score, top_model, top_fitted_sol, score_res
end












"""
NL segementation fitting
"""
function selection_NL_fixed_interval(
    data_testing::Matrix{Float64}, # dataset first row times second row OD
    name_well::String, # name of the well
    label_exp::String, #label of the experiment
    list_of_models::Vector{String}, # ode models to use
    list_lb_param::Any, # lower bound param
    list_ub_param::Any, # upper bound param
    interval_changepoints::Any;
    list_u0=list_lb_param .+ (list_ub_param .- list_lb_param) ./ 2,# initial guess param
    type_of_loss="L2", # type of used loss
    optmizator=BBO_adaptive_de_rand_1_bin_radiuslimited(), # selection of optimization method
    method_of_fitting="MCMC", # selection of sciml integrator
    smoothing=false,
    size_bootstrap=0.7,
    nrep=100,
    type_of_smoothing="lowess",
    thr_lowess=0.05,
    pt_avg=1,
    pt_smooth_derivative=0,
    multiple_scattering_correction=false, # if true uses the given calibration curve to fix the data
    method_multiple_scattering_correction="interpolation",
    calibration_OD_curve="NA", #  the path to calibration curve to fix the data
    beta_smoothing_ms=2.0, #  parameter of the AIC penality
    PopulationSize=300,
    maxiters=2000000,
    abstol=0.000000001,
    penality_CI=8.0)

    interval_changepoints = push!(interval_changepoints, data_testing[1, 1])
    interval_changepoints = push!(interval_changepoints, data_testing[1, end])
    interval_changepoints = sort(interval_changepoints)
    bc = [data_testing[1, end], data_testing[2, end]]
    param_out = Vector{Vector{Any}}()
    composed_sol = Type{Any}
    composed_time = Type{Any}
    loss_to_use = ""

    for i = (length(interval_changepoints)):-1:2

        if i == 2
            tspan_array = findall((data_testing[1, :] .<= interval_changepoints[i]))
            data_temp = Matrix(
                transpose(hcat(data_testing[1, tspan_array], data_testing[2, tspan_array])),
            )

            if type_of_loss == "RE"

                loss_to_use = "RE_fixed_end"

            elseif type_of_loss == "L2"

                loss_to_use = "L2_fixed_end"
            else
                loss_to_use = string(type_of_loss)
            end
            data_temp = hcat(data_temp, bc)

        elseif i == (length(interval_changepoints))


            tspan_array = findall((data_testing[1, :] .> interval_changepoints[i-1]))

            data_temp = Matrix(
                transpose(hcat(data_testing[1, tspan_array], data_testing[2, tspan_array])),
            )
            loss_to_use = string(type_of_loss)

        else
            tspan_array_1 = findall((data_testing[1, :] .> interval_changepoints[i-1]))
            tspan_array_2 = findall((data_testing[1, :] .<= interval_changepoints[i]))
            tspan_array = intersect(tspan_array_1, tspan_array_2)
            data_temp = Matrix(
                transpose(hcat(data_testing[1, tspan_array], data_testing[2, tspan_array])),
            )
            if type_of_loss == "RE"

                loss_to_use = "RE_fixed_end"

            elseif type_of_loss == "L2"

                loss_to_use = "L2_fixed_end"
            else


                loss_to_use = string(type_of_loss)
            end
            # imposing_bounduary condition
            data_temp = hcat(data_temp, bc)
        end

        model_selection_results = NL_model_selection(data_temp, # dataset first row times second row OD
            name_well, # name of the well
            label_exp, #label of the experiment
            list_of_models, # ode model to use
            list_lb_param, # lower bound param
            list_ub_param; # upper bound param
            method_of_fitting=method_of_fitting,
            nrep=nrep,
            list_u0=list_u0,# initial guess param
            optmizator=optmizator,
            display_plots=false, # display plots in julia or not
            size_bootstrap=size_bootstrap,
            pt_avg=pt_avg, # numebr of the point to generate intial condition
            pt_smooth_derivative=pt_smooth_derivative,
            smoothing=smoothing, # the smoothing is done or not?
            type_of_smoothing=type_of_smoothing,
            type_of_loss=loss_to_use, # type of used loss
            multiple_scattering_correction=multiple_scattering_correction, # if true uses the given calibration curve to fix the data
            method_multiple_scattering_correction=method_multiple_scattering_correction,
            calibration_OD_curve=calibration_OD_curve,  #  the path to calibration curve to fix the data
            PopulationSize=PopulationSize,
            maxiters=maxiters,
            abstol=abstol,
            thr_lowess=thr_lowess,
            write_res=false,
            beta_param=beta_smoothing_ms,
            penality_CI=penality_CI
        )

        # param of the best model
        temp_res_win = model_selection_results[2]


        time_sol = data_temp[1, :]
        sol_fin = model_selection_results[3]

        value_bonduary = sol_fin[1]
        time_bonduary = time_sol[1]

        bc = [time_bonduary, value_bonduary]

        sol_fin, index_not_zero = remove_negative_value(sol_fin)

        if i == length(interval_changepoints)
            composed_time = copy(time_sol[index_not_zero])
            composed_sol = copy(sol_fin)
            temp_res_win = push!(temp_res_win, i - 1)
            param_out = push!(param_out, temp_res_win)
        else


            composed_time = vcat(time_sol[index_not_zero], composed_time)
            composed_sol = vcat(sol_fin, composed_sol)
            temp_res_win = push!(temp_res_win, i - 1)
            param_out = push!(param_out, temp_res_win)
        end
    end


    return param_out, composed_sol, composed_time

end

function selection_NL_maxiumum_change_points(
    data_testing::Matrix{Float64}, # dataset first row times second row OD
    name_well::String, # name of the well
    label_exp::String, #label of the experiment
    list_of_models::Vector{String}, # ode models to use
    list_lb_param::Any, # lower bound param
    list_ub_param::Any, # upper bound param
    n_change_points::Int;
    list_u0=list_lb_param .+ (list_ub_param .- list_lb_param) ./ 2,# initial guess param
    type_of_loss="L2_fixed_CI", # type of used loss
    optmizator=BBO_adaptive_de_rand_1_bin_radiuslimited(), # selection of optimization method
    method_of_fitting="MCMC", # selection of sciml integrator
    type_of_detection="sliding_win",
    type_of_curve="original",
    smoothing=false,
    nrep=100,
    type_of_smoothing="lowess",
    thr_lowess=0.05,
    pt_avg=1,
    save_plot=false, # do plots or no
    display_plots=false,
    path_to_plot="NA", # where save plots
    win_size=7, # numebr of the point to generate intial condition
    pt_smooth_derivative=0,
    multiple_scattering_correction=false, # if true uses the given calibration curve to fix the data
    method_multiple_scattering_correction="interpolation",
    calibration_OD_curve="NA", #  the path to calibration curve to fix the data
    beta_smoothing_ms=2.0, #  parameter of the AIC penality
    method_peaks_detection="peaks_prominence",
    n_bins=40,
    PopulationSize=300,
    maxiters=2000000,
    abstol=0.000000001,
    dectect_number_cdp=true,
    fixed_cpd=false,
    penality_CI=8.0,
    size_bootstrap=0.7)

    top_aicc = 10^20
    top_param = Vector{Any}
    top_fit = Vector{Any}
    top_time = Vector{Any}
    top_intervals = Vector{Any}

    if multiple_scattering_correction == true
        data_testing = correction_OD_multiple_scattering(data_testing, calibration_OD_curve; method=method_multiple_scattering_correction)
    end

    if smoothing == true
        data_testing = smoothing_data(
            data;
            method=type_of_smoothing,
            pt_avg=pt_avg,
            thr_lowess=thr_lowess
        )
    end




    if dectect_number_cdp == true

        list_change_points_dev = cpd_local_detection(
            data_testing,
            n_change_points;
            type_of_detection=type_of_detection,
            type_of_curve=type_of_curve,
            pt_derivative=pt_smooth_derivative,
            size_win=win_size,
            method=method_peaks_detection,
            number_of_bin=n_bins,
        )

        combination_to_test = generation_of_combination_of_cpds(list_change_points_dev[2],
            n_fix=0)


    elseif fixed_cpd == true

        list_change_points_dev = cpd_local_detection(
            data_testing,
            n_change_points;
            type_of_detection=type_of_detection,
            type_of_curve=type_of_curve,
            pt_derivative=pt_smooth_derivative,
            size_win=win_size,
            method=method_peaks_detection,
            number_of_bin=n_bins,
        )

        combination_to_test = generation_of_combination_of_cpds(list_change_points_dev[2],
            n_fix=n_change_points)



    else
        list_change_points_dev = cpd_local_detection(
            data_testing,
            2 * n_change_points;
            type_of_detection=type_of_detection,
            type_of_curve=type_of_curve,
            pt_derivative=pt_smooth_derivative,
            size_win=win_size,
            method=method_peaks_detection,
            number_of_bin=n_bins,
        )

        combination_to_test = generation_of_combination_of_cpds(list_change_points_dev[2],
            n_fix=n_change_points)


    end

    for i in 1:eachindex(combination_to_test)[end]

        cpd_temp = sort(combination_to_test[i])

        res_this_combination = selection_NL_fixed_interval(
            data_testing, # dataset first row times second row OD
            name_well, # name of the well
            label_exp, #label of the experiment
            list_of_models, # ode models to use
            list_lb_param, # lower bound param
            list_ub_param, # upper bound param
            cpd_temp;
            list_u0=list_u0,
            type_of_loss=type_of_loss, # type of used loss
            optmizator=optmizator, # selection of optimization method
            method_of_fitting=method_of_fitting, # selection of sciml integrator
            smoothing=smoothing,
            type_of_smoothing=type_of_smoothing,
            thr_lowess=thr_lowess,
            pt_avg=pt_avg,
            nrep=nrep,
            size_bootstrap=size_bootstrap,
            pt_smooth_derivative=pt_smooth_derivative,
            multiple_scattering_correction=multiple_scattering_correction, # if true uses the given calibration curve to fix the data
            method_multiple_scattering_correction=method_multiple_scattering_correction,
            calibration_OD_curve=calibration_OD_curve, #  the path to calibration curve to fix the data
            beta_smoothing_ms=beta_smoothing_ms, #  parameter of the AIC penality
            PopulationSize=PopulationSize,
            maxiters=maxiters,
            abstol=abstol,
            penality_CI=penality_CI)


        n_param_full_model = sum([
            length(res_this_combination[1][kk][3:(end-5)]) for
            kk = 1:length(res_this_combination[1])
        ]) + n_change_points





        AICc_full_model = AICc_evaluation(n_param_full_model, beta_smoothing_ms, res_this_combination[3], res_this_combination[2])


        if i == 1

            top_aicc = copy(AICc_full_model)
            top_param = copy(res_this_combination[1])
            top_fit = copy(res_this_combination[2])
            top_time = copy(res_this_combination[3])
            top_intervals = copy(cpd_temp)


        elseif AICc_full_model < top_aicc

            top_aicc = copy(AICc_full_model)
            top_param = copy(res_this_combination[1])
            top_fit = copy(res_this_combination[2])
            top_time = copy(res_this_combination[3])
            top_intervals = copy(cpd_temp)
        end

    end
    if display_plots
        if_display = display
    else
        if_display = identity
    end

    if save_plot
        mkpath(path_to_plot)
    end

    if_display(
        Plots.scatter(
            data_testing[1, :],
            data_testing[2, :],
            xlabel="Time",
            ylabel="Arb. Units",
            label=["Data " nothing],
            markersize=1,
            color=:black,
            title=string(label_exp, " ", name_well),
        ),
    )
    if_display(
        Plots.vline!(
            top_intervals[1:end],
            c=:black,
            label=["change points" nothing],
        ),
    )
    if_display(
        Plots.plot!(
            reduce(vcat, top_time),
            reduce(vcat, top_fit),
            xlabel="Time",
            ylabel="Arb. Units",
            label=[" fitting " nothing],
            color=:red,
            title=string(label_exp, " fitting ", name_well),
        ),
    )
    if save_plot
        png(
            string(
                path_to_plot,
                label_exp,
                "_model_selection_seg_",
                n_change_points,
                "_",
                name_well,
                ".png",
            ),
        )
    end

    return top_param, sort(top_intervals), top_fit, top_time
end