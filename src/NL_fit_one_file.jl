
"""
fitting dataset function NL
"""

function fit_NL_model_file(
    label_exp::String, #label of the experiment
    path_to_data::String, # path to the folder to analyze
    path_to_annotation::String,# path to the annotation of the wells
    model::String, # string of the used model
    lb_param::Vector{Float64},# array of the array of the lower bound of the parameters
    ub_param::Vector{Float64}; # array of the array of the upper bound of the parameters
    u0=lb_param .+ (ub_param .- lb_param) ./ 2,# initial guess param
    method_NL_fit="MCMC",
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
    penality_CI=8.0
)


    if write_res == true
        mkpath(path_to_results)
    end

    if save_plots == true
        mkpath(path_to_plot)
    end

    parameter_of_optimization = initialize_df_results_ode_custom(lb_param)




    annotation = CSV.File(string(path_to_annotation), header=false)
    names_of_annotated_df = [annotation[l][1] for l in eachindex(annotation)]
    # selcting blank wells
    properties_of_annotation = [annotation[l][2] for l in eachindex(annotation)]
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
                string(well_name), # name of the well
                label_exp, #label of the experiment
                model, # ode model to use
                lb_param, # lower bound param
                ub_param; # upper bound param
                nrep=nrep,
                u0=u0,# initial guess param
                optmizator=optmizator,
                display_plots=display_plots, # display plots in julia or not
                save_plot=save_plots,
                size_bootstrap=size_bootstrap,
                path_to_plot=path_to_plot, # where save plots
                pt_avg=pt_avg, # numebr of the point to generate intial condition
                pt_smooth_derivative=pt_smooth_derivative,
                smoothing=smoothing, # the smoothing is done or not?
                type_of_smoothing=type_of_smoothing,
                type_of_loss=loss_type, # type of used loss
                multiple_scattering_correction=multiple_scattering_correction, # if true uses the given calibration curve to fix the data
                method_multiple_scattering_correction=method_multiple_scattering_correction,
                calibration_OD_curve=calibration_OD_curve,  #  the path to calibration curve to fix the data
                PopulationSize=PopulationSize,
                maxiters=maxiters,
                abstol=abstol,
                thr_lowess=thr_lowess,
                write_res=write_res,
                penality_CI=penality_CI)


        elseif method_NL_fit == "Morris_sensitivity"


            temp_results_1 = fit_NL_model_with_sensitivity(data, # dataset first row times second row OD
            string(well_name), # name of the well
                label_exp, #label of the experiment
                model, # ode model to use
                lb_param, # lower bound param
                ub_param; # upper bound param
                nrep=nrep,
                optmizator=optmizator,
                display_plots=display_plots, # display plots in julia or not
                save_plot=save_plots,
                path_to_plot=path_to_plot, # where save plots
                pt_avg=pt_avg, # numebr of the point to generate intial condition
                pt_smooth_derivative=pt_smooth_derivative,
                smoothing=smoothing, # the smoothing is done or not?
                type_of_smoothing=type_of_smoothing,
                type_of_loss=loss_type, # type of used loss
                multiple_scattering_correction=multiple_scattering_correction, # if true uses the given calibration curve to fix the data
                method_multiple_scattering_correction=method_multiple_scattering_correction,
                calibration_OD_curve=calibration_OD_curve,  #  the path to calibration curve to fix the data
                PopulationSize=PopulationSize,
                maxiters=maxiters,
                abstol=abstol,
                thr_lowess=thr_lowess,
                write_res=write_res,
                penality_CI=penality_CI)
        elseif method_NL_fit == "MCMC"


            temp_results_1 = fit_NL_model_MCMC_intialization(data, # dataset first row times second row OD
            string(well_name), # name of the well
                label_exp, #label of the experiment
                model, # ode model to use
                lb_param, # lower bound param
                ub_param; # upper bound param
                nrep=nrep,
                optmizator=optmizator,
                display_plots=display_plots, # display plots in julia or not
                save_plot=save_plots,
                path_to_plot=path_to_plot, # where save plots
                pt_avg=pt_avg, # numebr of the point to generate intial condition
                pt_smooth_derivative=pt_smooth_derivative,
                smoothing=smoothing, # the smoothing is done or not?
                type_of_smoothing=type_of_smoothing,
                type_of_loss=loss_type, # type of used loss
                multiple_scattering_correction=multiple_scattering_correction, # if true uses the given calibration curve to fix the data
                method_multiple_scattering_correction=method_multiple_scattering_correction,
                calibration_OD_curve=calibration_OD_curve,  #  the path to calibration curve to fix the data
                PopulationSize=PopulationSize,
                maxiters=maxiters,
                abstol=abstol,
                thr_lowess=thr_lowess,
                penality_CI=penality_CI)


        else



            temp_results_1 = fit_NL_model(data, # dataset first row times second row OD
            string(well_name), # name of the well
                label_exp, #label of the experiment
                model, # ode model to use
                lb_param, # lower bound param
                ub_param; # upper bound param
                u0=u0,# initial guess param
                optmizator=optmizator,
                display_plots=display_plots, # display plots in julia or not
                save_plot=save_plots,
                path_to_plot=path_to_plot, # where save plots
                pt_avg=pt_avg, # numebr of the point to generate intial condition
                pt_smooth_derivative=pt_smooth_derivative,
                smoothing=smoothing, # the smoothing is done or not?
                type_of_smoothing=type_of_smoothing,
                type_of_loss=loss_type, # type of used loss
                multiple_scattering_correction=multiple_scattering_correction, # if true uses the given calibration curve to fix the data
                method_multiple_scattering_correction=method_multiple_scattering_correction,
                calibration_OD_curve=calibration_OD_curve,  #  the path to calibration curve to fix the data
                PopulationSize=PopulationSize,
                maxiters=maxiters,
                abstol=abstol,
                thr_lowess=thr_lowess,
                penality_CI=penality_CI
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





function fit_NL_model_selection_file(
    label_exp::String, #label of the experiment
    path_to_data::String, # path to the folder to analyze
    path_to_annotation::String,# path to the annotation of the wells
    list_model_function::Any, # ode model to use
    list_lb_param::Vector{Float64}, # lower bound param
    list_ub_param::Vector{Float64}; # upper bound param
    method_of_fitting="MCMC",
    nrep=100,
    list_u0=lb_param .+ (ub_param .- lb_param) ./ 2,# initial guess param
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
    beta_param=2.0,
    penality_CI=8.0,
    size_bootstrap=0.7,)


    if write_res == true
        mkpath(path_to_results)
    end

    if save_plots == true
        mkpath(path_to_plot)
    end
    parameter_of_optimization = initialize_res_ms(list_ub_param)




    annotation = CSV.File(string(path_to_annotation), header=false)
    names_of_annotated_df = [annotation[l][1] for l in eachindex(annotation)]
    # selcting blank wells
    properties_of_annotation = [annotation[l][2] for l in eachindex(annotation)]
    list_of_blank = names_of_annotated_df[findall(x -> x == "b", properties_of_annotation)]
    list_of_discarded = names_of_annotated_df[findall(x -> x == "X", properties_of_annotation)]
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


        temp_results_1 = NL_model_selection(data, # dataset first row times second row OD
            string(well_name), # name of the well
            label_exp, #label of the experiment
            list_model_function, # ode model to use
            list_lb_param, # lower bound param
            list_ub_param; # upper bound param
            method_of_fitting=method_of_fitting,
            nrep=nrep,
            list_u0=list_u0,# initial guess param
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
            multiple_scattering_correction=multiple_scattering_correction, # if true uses the given calibration curve to fix the data
            method_multiple_scattering_correction=method_multiple_scattering_correction,
            calibration_OD_curve=calibration_OD_curve,  #  the path to calibration curve to fix the data
            PopulationSize=PopulationSize,
            maxiters=maxiters,
            abstol=abstol,
            thr_lowess=thr_lowess,
            write_res=false,
            beta_param=beta_param,
            penality_CI=penality_CI
        )



        data = Matrix(data)




        if verbose == true
            println("the results are:")
            println(temp_results_1[2])
        end

        parameter_of_optimization = hcat(parameter_of_optimization, temp_results_1[2])

    end


    if write_res == true

        CSV.write(
            string(path_to_results, label_exp, "_parameters_", model, ".csv"),
            Tables.table(Matrix(parameter_of_optimization)),
        )


    end
    return parameter_of_optimization




end



function fit_NL_segmentation_file(
    label_exp::String, #label of the experiment
    path_to_data::String, # path to the folder to analyze
    path_to_annotation::String,# path to the annotation of the wells
    list_model_function::Any, # ode model to use
    list_lb_param::Vector{Vector{Float64}}, # lower bound param
    list_ub_param::Vector{Vector{Float64}}, # upper bound param
    n_change_points::Int;
    method_of_fitting="MCMC",
    nrep=100,
    list_u0=lb_param .+ (ub_param .- lb_param) ./ 2,# initial guess param
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
    size_bootstrap=0.7,
    thr_lowess=0.05,
    dectect_number_cdp=true,
    type_of_detection="sliding_win",
    type_of_curve="original",
    fixed_cpd=false,
    penality_CI=8.0,
    beta_smoothing_ms=2.0,
    win_size=7, # number of the point of cpd sliding win
    n_bins=40,
)


    if write_res == true
        mkpath(path_to_results)
    end

    if save_plots == true
        mkpath(path_to_plot)
    end
    parameter_of_optimization = initialize_res_ms(list_ub_param, number_of_segment=n_change_points)




    annotation = CSV.File(string(path_to_annotation), header=false)
    names_of_annotated_df = [annotation[l][1] for l in eachindex(annotation)]
    # selcting blank wells
    properties_of_annotation = [annotation[l][2] for l in eachindex(annotation)]
    list_of_blank = names_of_annotated_df[findall(x -> x == "b", properties_of_annotation)]
    list_of_discarded = names_of_annotated_df[findall(x -> x == "X", properties_of_annotation)]
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


        temp_results_1 = selection_NL_maxiumum_change_points(
            data, # dataset first row times second row OD
            string(well_name), # name of the well
            label_exp, #label of the experiment
            list_model_function, # ode models to use
            list_lb_param, # lower bound param
            list_ub_param, # upper bound param
            n_change_points;
            list_u0=list_u0,# initial guess param
            type_of_loss=loss_type, # type of used loss
            optmizator=optmizator, # selection of optimization method
            method_of_fitting=method_of_fitting, # selection of sciml integrator
            type_of_detection=type_of_detection,
            type_of_curve=type_of_curve,
            smoothing=smoothing,
            nrep=nrep,
            type_of_smoothing=type_of_smoothing,
            thr_lowess=thr_lowess,
            pt_avg=pt_avg,
            save_plot=save_plots, # do plots or no
            display_plots=display_plots,
            path_to_plot=path_to_plot, # where save plots
            win_size=win_size, # numebr of the point to generate intial condition
            pt_smooth_derivative=pt_smooth_derivative,
            multiple_scattering_correction=multiple_scattering_correction, # if true uses the given calibration curve to fix the data
            method_multiple_scattering_correction=method_multiple_scattering_correction,
            calibration_OD_curve=calibration_OD_curve, #  the path to calibration curve to fix the data
            beta_smoothing_ms=beta_smoothing_ms, #  parameter of the AIC penality
            n_bins=n_bins,
            PopulationSize=PopulationSize,
            maxiters=maxiters,
            abstol=abstol,
            dectect_number_cdp=dectect_number_cdp,
            fixed_cpd=fixed_cpd,
            penality_CI=penality_CI,
            size_bootstrap=size_bootstrap)


        data = Matrix(data)




        if verbose == true
            println("the results are:")
            println(temp_results_1[1])
        end

        results_to_bind = expand_res(
            temp_results_1[1],
            list_lb_param,
            string(well_name),
            label_exp;
            number_of_segment=length(temp_results_1[1]))

        parameter_of_optimization = hcat(parameter_of_optimization, results_to_bind)

    end


    if write_res == true

        CSV.write(
            string(path_to_results, label_exp, "_parameters_ segmentation_NL.csv"),
            Tables.table(Matrix(parameter_of_optimization)),
        )


    end
    return parameter_of_optimization




end