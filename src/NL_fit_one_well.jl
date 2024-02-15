



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
    blank_array=zeros(100), # data of all blanks
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
    blank_array=zeros(100), # data of all blanks
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
    blank_array=zeros(100), # data of all blanks
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
    mean_param = [mean(fin_param[i, 2:end]) for i in 3:size(fin_param, 1)]
    sd_param = [std(fin_param[i, 2:end]) for i in 3:size(fin_param, 1)]

    return best_res_param, best_fitted_model, fin_param, mean_param, sd_param
end

















function NL_model_selection(data::Matrix{Float64}, # dataset first row times second row OD
    name_well::String, # name of the well
    label_exp::String, #label of the experiment
    list_model_function::Any, # ode model to use
    list_lb_param::Vector{Float64}, # lower bound param
    list_ub_param::Vector{Float64}; # upper bound param
    method_of_fitting="Bootstrap",
    nrep=100,
    list_u0=lb_param .+ (ub_param .- lb_param) ./ 2,# initial guess param
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
    blank_array=zeros(100), # data of all blanks
    multiple_scattering_correction=false, # if true uses the given calibration curve to fix the data
    method_multiple_scattering_correction="interpolation",
    calibration_OD_curve="NA",  #  the path to calibration curve to fix the data
    PopulationSize=300,
    maxiters=2000000,
    abstol=0.00001,
    thr_lowess=0.05,
    write_res=false,
    beta_param=2.0,
)
    score_res = ["model", "loss", "AIC"]
    top_score = 10^20
    top_model = Any{}

    for mm in 1:size(list_model_function, 1)

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
                blank_array=blank_array, # data of all blanks
                multiple_scattering_correction=multiple_scattering_correction, # if true uses the given calibration curve to fix the data
                method_multiple_scattering_correction=method_multiple_scattering_correction,
                calibration_OD_curve=calibration_OD_curve,  #  the path to calibration curve to fix the data
                PopulationSize=PopulationSize,
                maxiters=maxiters,
                abstol=abstol,
                thr_lowess=thr_lowess,
                write_res=write_res
            )

            n_param = length(lb_param)

            temp_AIC = AICc_evaluation(n_param, beta_param, data, temp_res[2])
            temp = [model_to_test, temp_res[end], temp_AIC]
            score_res = hcat(score_res, temp_AIC)
            if mm == 1

                top_score = temp_AIC
                top_model = temp_res[1]
                top_fitted_sol = temp_res[2]

            elseif top_score > temp_AIC
                top_score = temp_AIC
                top_model = temp_res[1]
                top_fitted_sol = temp_res[2]

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
                blank_array=blank_array, # data of all blanks
                multiple_scattering_correction=multiple_scattering_correction, # if true uses the given calibration curve to fix the data
                method_multiple_scattering_correction=method_multiple_scattering_correction,
                calibration_OD_curve=calibration_OD_curve,  #  the path to calibration curve to fix the data
                PopulationSize=PopulationSize,
                maxiters=maxiters,
                abstol=abstol,
                thr_lowess=thr_lowess,
                write_res=write_res
            )

            n_param = length(lb_param)

            temp_AIC = AICc_evaluation(n_param, beta_param, data, temp_res[2])
            temp = [model_to_test, temp_res[end], temp_AIC]

            score_res = hcat(score_res, temp_AIC)
            if mm == 1

                top_score = temp_AIC
                top_model = temp_res[1]
                top_fitted_sol = temp_res[2]


            elseif top_score > temp_AIC
                top_score = temp_AIC
                top_model = temp_res[1]
                top_fitted_sol = temp_res[2]

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
                blank_array=blank_array, # data of all blanks
                multiple_scattering_correction=multiple_scattering_correction, # if true uses the given calibration curve to fix the data
                method_multiple_scattering_correction=method_multiple_scattering_correction,
                calibration_OD_curve=calibration_OD_curve,  #  the path to calibration curve to fix the data
                PopulationSize=PopulationSize,
                maxiters=maxiters,
                abstol=abstol,
                thr_lowess=thr_lowess,
            )


            n_param = length(lb_param)

            temp_AIC = AICc_evaluation(n_param, beta_param, data, temp_res[2])
            temp = [model_to_test, temp_res[end], temp_AIC]

            score_res = hcat(score_res, temp_AIC)
            if mm == 1

                top_score = temp_AIC
                top_model = temp_res[1]
                top_fitted_sol = temp_res[2]


            elseif top_score > temp_AIC
                top_score = temp_AIC
                top_model = temp_res[1]
                top_fitted_sol = temp_res[2]

            end


        end


    end

    # plotting if required
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


    return top_score, top_model, top_fitted_sol
end







