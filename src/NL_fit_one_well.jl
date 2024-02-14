



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
        select_loss_function_NL(type_of_loss, data, model_function)

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
    write_res=false
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
        select_loss_function_NL(type_of_loss, data, model_function)



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


        max_th_gr = maximum(specific_gr_evaluation(Matrix(data_th), pt_smooth_derivative))

        # max empirical gr
        loss_value = sol.objective


        res_param = [[name_well, model_string], [sol[1:end]], [max_th_gr, max_em_gr, loss_value]]
        res_param = reduce(vcat, reduce(vcat, res_param))



        fin_param = hcat(fin_param, res_param)


    end


    index_best = findmin(fin_param[end, 2:end])[2]

    best_res_param = fin_param[:, index_best]
    println(best_res_param[3:(end-3)])

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


    return best_res_param, best_fitted_model, fin_param
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
    write_res=false
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
            select_loss_function_NL(type_of_loss, data_to_fit, model_function)
        prob = OptimizationProblem(loss_function, u0, data_to_fit, lb=lb_param, ub=ub_param)

        # Solve the optimization problem
        sol = solve(prob, optmizator, PopulationSize=PopulationSize, maxiters=maxiters, abstol=abstol)
        # evaluate the fitted  model
        fitted_model = model_function(sol, data_to_fit[1, :])
        sol_fin, index_not_zero = remove_negative_value(fitted_model)

        data_th = transpose(hcat(data_to_fit[1, index_not_zero], sol_fin))


        max_th_gr = maximum(specific_gr_evaluation(Matrix(data_th), pt_smooth_derivative))

        # max empirical gr
        loss_value = sol.objective


        res_param = [[name_well, model_string], [sol[1:end]], [max_th_gr, max_em_gr, loss_value]]
        res_param = reduce(vcat, reduce(vcat, res_param))



        fin_param = hcat(fin_param, res_param)


    end


    index_best = findmin(fin_param[end, 2:end])[2]

    best_res_param = fin_param[:, index_best]

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

            if top_score > temp_AIC
                top_score = temp_AIC
                top_model = temp_res[1]

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
            if top_score > temp_AIC
                top_score = temp_AIC
                top_model = temp_res[1]

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
            if top_score > temp_AIC
                top_score = temp_AIC
                top_model = temp_res[1]

            end


        end


    end

    return score_res, top_model
end








"""
fitting dataset function ODE
"""

function fit_NL_model_file(
    label_exp::String, #label of the experiment
    path_to_data::String, # path to the folder to analyze
    path_to_annotation::String,# path to the annotation of the wells
    model::String, # string of the used model
    lb_param::Vector{Float64},# array of the array of the lower bound of the parameters
    ub_param::Vector{Float64}; # array of the array of the upper bound of the parameters
    u0=lb_param .+ (ub_param .- lb_param) ./ 2,# initial guess param
    method_NL_fit="Bootstrap",
    nrep=100,
    optmizator=BBO_adaptive_de_rand_1_bin_radiuslimited(), # selection of optimization method 
    path_to_results="NA", # path where save results
    path_to_plot="NA", # path where to save Plots
    loss_type="RE", # string of the type of the used loss
    smoothing=false, # 1 do smoothing of data with rolling average
    type_of_smoothing="lowess",
    display_plots=true,# display plots in julia or not
    save_plots=false,
    verbose=false, # 1 true verbose
    write_res=false, # write results
    pt_avg=1, # number of points to do smoothing average
    pt_smooth_derivative=7, # number of points to do ssmooth_derivative
    do_blank_subtraction="avg_blank", # string on how to use blank (NO,avg_subtraction,time_avg)
    avg_replicate=false, # if true the average between replicates is fitted. If false all replicate are fitted indipendelitly
    correct_negative="thr_correction", # if "thr_correction" it put a thr on the minimum value of the data with blank subracted, if "blank_correction" uses blank distrib to impute negative values
    thr_negative=0.01,  # used only if correct_negative == "thr_correction"
    multiple_scattering_correction=false, # if true uses the given calibration curve to fix the data
    method_multiple_scattering_correction="interpolation",
    calibration_OD_curve="NA",  #  the path to calibration curve to fix the data
    PopulationSize=300,
    maxiters=2000000,
    abstol=0.00001,
    thr_lowess=0.05,
)


    if write_res == true
        mkpath(path_to_results)
    end

    if save_plots == true
        mkpath(path_to_plot)
    end

    parameter_of_optimization = initialize_df_results_ode_custom(lb_param)




    annotation = CSV.File(string(path_to_annotation), header=false)
    names_of_annotated_df = [annotation[l][1] for l = 1:length(annotation)]
    # selcting blank wells
    properties_of_annotation = [annotation[l][2] for l = 1:length(annotation)]
    list_of_blank = names_of_annotated_df[findall(x -> x == "b", properties_of_annotation)]
    list_of_discarded =
        names_of_annotated_df[findall(x -> x == "X", properties_of_annotation)]
    list_of_blank = Symbol.(list_of_blank)

    # reading files
    dfs_data = CSV.File(path_to_data)

    # shaping df for the inference
    names_of_cols = propertynames(dfs_data)

    # excluding blank data and discarded wells
    names_of_cols = filter!(e -> !(e in list_of_blank), names_of_cols)
    if length(list_of_discarded) > 0
        names_of_cols = filter!(e -> !(e in list_of_discarded), names_of_cols)
    end

    times_data = dfs_data[names_of_cols[1]]
    if length(list_of_blank) > 0
        blank_array = reduce(vcat, [(dfs_data[k]) for k in list_of_blank])
        blank_array = convert(Vector{Float64}, blank_array)

        blank_value = blank_subtraction(
            dfs_data,
            list_of_blank;
            method=do_blank_subtraction
        )
    else
        blank_value = 0.0
    end


    ## considering replicates
    list_replicate = unique(properties_of_annotation)
    list_replicate = filter!(e -> e != "b", list_replicate)

    if avg_replicate == true


        dfs_data, names_of_cols = average_replicate(dfs_data, times_data, properties_of_annotation, names_of_annotated_df)


    end

    # for on the columns to analyze

    for well_name in names_of_cols[2:end]




        if avg_replicate == true

            data_values = copy(dfs_data[!, well_name])

        else
            data_values = copy(dfs_data[well_name])
        end

        # blank subtraction 
        data_values = data_values .- blank_value

        data = Matrix(transpose(hcat(times_data, data_values)))


        # correcting negative values after blank subtraction
        data = negative_value_correction(data,
            blank_array;
            method=correct_negative,
            thr_negative=thr_negative,)

        data = Matrix(data)

        # defining time steps of the inference
        if method_NL_fit == "Bootstrap"



            temp_results_1 = fit_NL_model_bootstrap(data, # dataset first row times second row OD
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
                type_of_loss=loss_type, # type of used loss
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


        elseif method_NL_fit == "Morris_sensitivity"


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
                type_of_loss=loss_type, # type of used loss
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
                type_of_loss=loss_type, # type of used loss
                blank_array=blank_array, # data of all blanks
                multiple_scattering_correction=multiple_scattering_correction, # if true uses the given calibration curve to fix the data
                method_multiple_scattering_correction=method_multiple_scattering_correction,
                calibration_OD_curve=calibration_OD_curve,  #  the path to calibration curve to fix the data
                PopulationSize=PopulationSize,
                maxiters=maxiters,
                abstol=abstol,
                thr_lowess=thr_lowess,
            )





        end


        data = Matrix(data)




        if verbose == true
            println("the results are:")
            println(temp_results_1[1])
        end

        parameter_of_optimization = hcat(parameter_of_optimization, temp_results_1[1])

    end


    if write_res == true

        CSV.write(
            string(path_to_results, label_exp, "_parameters_", model, ".csv"),
            Tables.table(Matrix(parameter_of_optimization)),
        )


    end
    return parameter_of_optimization




end
