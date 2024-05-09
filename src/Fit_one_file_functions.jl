#######################################################################

"""
plotting functions
"""
function plot_data(
    label_exp::String, #label of the experiment
    path_to_data::String; # path to the folder to analyze
    path_to_annotation::Any = missing,# path to the annotation of the wells
    path_to_plot="NA", # path where to save Plots
    display_plots=true,# display plots in julia or not
    save_plots=false, # save the plot or not
    overlay_plots=true, # true a single plot for all dataset false one plot per well
    do_blank_subtraction="NO", # string on how to use blank (NO,avg_subtraction,time_avg)
    avg_replicate=false, # if true the average between replicates
    correct_negative="thr_correction", # if "thr_correction" it put a thr on the minimum value of the data with blank subracted, if "blank_correction" uses blank distrib to impute negative values
    thr_negative=0.01 ,
    blank_value = 0.0,
    blank_array = [0.0],)
    """
    function that plot a dataset
    """
  
    names_of_annotated_df,properties_of_annotation,list_of_blank, list_of_discarded = reading_annotation(path_to_annotation)
    # reading files
    dfs_data = CSV.File(path_to_data)

    # shaping df for the inference
    names_of_cols = propertynames(dfs_data)

    # excluding blank data and discarded wells
    if length(list_of_blank) > 0
        names_of_cols = filter!(e -> !(e in list_of_blank), names_of_cols)
    end

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

    end


    ## considering replicates
    list_replicate = unique(properties_of_annotation)
    list_replicate = filter!(e -> e != "b", list_replicate)

    if avg_replicate == true


        dfs_data, names_of_cols = average_replicate(dfs_data, times_data, properties_of_annotation, names_of_annotated_df)


    end
    # creating the folder of data if one wants to save
    if save_plots == true
        mkpath(path_to_plot)
    end

    for well_name in names_of_cols[2:end]
        name_well = string(well_name)
        if avg_replicate == true
            data_values = copy(dfs_data[!, well_name])
        else
            data_values = copy(dfs_data[well_name])
        end

        # blank subtraction
        data_values = data_values .- blank_value
        index_missing = findall(ismissing, data_values)
        index_tot =  eachindex(data_values)
        index_tot =  setdiff(index_tot,index_missing)
        data = Matrix(transpose(hcat(times_data[index_tot], data_values[index_tot])))
        # correcting negative values after blank subtraction
        data = negative_value_correction(data,
            blank_array;
            method=correct_negative,
            thr_negative=thr_negative,)

        if display_plots
            if_display = display
        else
            if_display = identity
        end

        if overlay_plots
            if well_name == names_of_cols[2]
                if_display(
                    Plots.plot(
                        data[1, :],
                        data[2, :],
                        xlabel="Time",
                        ylabel="Arb. Units",
                        label=[name_well],
                        title=string(label_exp),
                        legend=:outertopright,
                    ),
                )
            else
                if_display(
                    Plots.plot!(
                        data[1, :],
                        data[2, :],
                        xlabel="Time",
                        ylabel="Arb. Units",
                        label=[name_well],
                        title=string(label_exp),
                        legend=:outertopright,
                    ),
                )
            end

            if save_plots
                png(string(path_to_plot, label_exp, ".png"))
            end
        else
            # save & not plot single plot
            if_display(
                Plots.plot(
                    data[1, :],
                    data[2, :],
                    xlabel="Time",
                    ylabel="Arb. Units",
                    label=["Data " nothing],
                    color=:black,
                    title=string(label_exp, " ", name_well),
                ),
            )
            if save_plots
                png(string(path_to_plot, label_exp, "_", name_well, ".png"))
            end
        end

    end
end

function fit_one_file_Log_Lin(
    label_exp::String, #label of the experiment
    path_to_data::String; # path to the folder to analyze
    path_to_annotation::Any = missing,# path to the annotation of the wells
    path_to_results="NA",# path where save results
    path_to_plot="NA",# path where to save Plots
    display_plots=true,# display plots in julia or not
    save_plots=false, # save the plot or not    verbose=false, # 1 true verbose
    write_res=false, # write results
    type_of_smoothing="rolling_avg", # option, NO, gaussian, rolling avg
    pt_avg=7, # number of points to do smoothing average
    pt_smoothing_derivative=7, # number of poits to smooth the derivative
    pt_min_size_of_win=7, # minimum size of the exp windows in number of smooted points
    type_of_win="maximum", # how the exp. phase win is selected, "maximum" of "global_thr"
    threshold_of_exp=0.9, # threshold of growth rate in quantile to define the exp windows
    do_blank_subtraction="avg_blank", # string on how to use blank (NO,avg_subtraction,time_avg)
    avg_replicate=false, # if true the average between replicates is fitted. If false all replicate are fitted indipendelitly
    correct_negative="thr_correction", # if "thr_correction" it put a thr on the minimum value of the data with blank subracted, if "blank_correction" uses blank distrib to impute negative values
    thr_negative=0.01, # used only if correct_negative == "thr_correction"
    multiple_scattering_correction=false, # if true uses the given calibration curve to fix the data
    method_multiple_scattering_correction="interpolation",
    calibration_OD_curve="NA", #  the path to calibration curve to fix the data
    thr_lowess=0.05, # keyword argument of lowees smoothing
    verbose=false,
    blank_value = 0.0,
    blank_array = [0.0],)








    # TEMPORARY results df
    results_Log_Lin = [
        "label_exp",
        "well_name",
        "t_start",
        "t_end",
        "t_of_max",
        "empirical_max_Growth_rate",
        "Growth_rate",
        "sigma_gr",
        "dt",
        "sigma_confidence_dt_upper",
        "sigma_confidence_dt_lower",
        "intercept",
        "sigma_intercept",
        "Pearson_correlation",
    ]

    if save_plots == true
        mkpath(path_to_plot)
    end
    if write_res == true
        mkpath(path_to_results)
    end
    names_of_annotated_df, properties_of_annotation,list_of_blank, list_of_discarded = reading_annotation(path_to_annotation)


    # reading files
    dfs_data = CSV.File(path_to_data)

    # shaping df for the inference
    names_of_cols = propertynames(dfs_data)

    # excluding blank data and discarded wells
    if length(list_of_blank) > 0
        names_of_cols = filter!(e -> !(e in list_of_blank), names_of_cols)
    end
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

        index_missing = findall(ismissing, data_values)
        index_tot =  eachindex(data_values)
        index_tot =  setdiff(index_tot,index_missing)
        data = Matrix(transpose(hcat(times_data[index_tot], data_values[index_tot])))

        # correcting negative values after blank subtraction
        data = negative_value_correction(data,
            blank_array;
            method=correct_negative,
            thr_negative=thr_negative,)

        data = Matrix(transpose(hcat(data[1, :], data[2, :])))

        # inference


        temp_results_1 = fitting_one_well_Log_Lin(
            data, # dataset first row times second row OD
            string(well_name), # name of the well
            label_exp; #label of the experiment
            display_plots=display_plots,# display plots in julia or not
            save_plot=save_plots, # save the plot or not    verbose=false, # 1 true verbose        
            path_to_plot=path_to_plot, # where save plots
            type_of_smoothing=type_of_smoothing, # option, NO, gaussian, rolling avg
            pt_avg=pt_avg, # number of the point for rolling avg not used in the other cases
            pt_smoothing_derivative=pt_smoothing_derivative, # number of poits to smooth the derivative
            pt_min_size_of_win=pt_min_size_of_win, # minimum size of the exp windows in number of smooted points
            type_of_win=type_of_win, # how the exp. phase win is selected, "maximum" of "global_thr"
            threshold_of_exp=threshold_of_exp, # threshold of growth rate in quantile to define the exp windows
            multiple_scattering_correction=multiple_scattering_correction, # if true uses the given calibration curve to fix the data
            method_multiple_scattering_correction=method_multiple_scattering_correction,
            calibration_OD_curve=calibration_OD_curve, #  the path to calibration curve to fix the data
            thr_lowess=thr_lowess,
        )

        if verbose == true
            println("the results are:")
            println(temp_results_1)
        end

        results_Log_Lin = hcat(results_Log_Lin, temp_results_1)

    end

    if write_res == true

        CSV.write(
            string(path_to_results, label_exp, "_results.csv"),
            Tables.table(Matrix(results_Log_Lin)),
        )

    end

    return results_Log_Lin




end



"""
fitting dataset function ODE
"""

function fit_file_ODE(
    label_exp::String, #label of the experiment
    path_to_data::String, # path to the folder to analyze
    model::String, # string of the used model
    lb_param::Vector{Float64},# array of the array of the lower bound of the parameters
    ub_param::Vector{Float64}; # array of the array of the upper bound of the parameters
    path_to_annotation::Any = missing,# path to the annotation of the wells
    optmizator=BBO_adaptive_de_rand_1_bin_radiuslimited(), # selection of optimization method 
    integrator=Tsit5(), # selection of sciml integrator
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
    blank_value = 0.0,
    blank_array = [0.0],
)


    if write_res == true
        mkpath(path_to_results)
    end

    if save_plots == true
        mkpath(path_to_plot)
    end

    parameter_of_optimization = initialize_df_results(model)



    names_of_annotated_df,properties_of_annotation,list_of_blank, list_of_discarded = reading_annotation(path_to_annotation)


    # reading files
    dfs_data = CSV.File(path_to_data)

    # shaping df for the inference
    names_of_cols = propertynames(dfs_data)

    # excluding blank data and discarded wells
    if length(list_of_blank) > 0
        names_of_cols = filter!(e -> !(e in list_of_blank), names_of_cols)
    end

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

        index_missing = findall(ismissing, data_values)
        index_tot =  eachindex(data_values)
        index_tot =  setdiff(index_tot,index_missing)
        data = Matrix(transpose(hcat(times_data[index_tot], data_values[index_tot])))

        # correcting negative values after blank subtraction
        data = negative_value_correction(data,
            blank_array;
            method=correct_negative,
            thr_negative=thr_negative,)


        # defining time steps of the inference

        max_t = data[1, end]
        min_t = data[1, 1]


        data = Matrix(data)



        # inference


        temp_results_1 = fitting_one_well_ODE_constrained(
            data, # dataset first row times second row OD
            string(well_name), # name of the well
            label_exp, #label of the experiment
            model, # ode model to use 
            lb_param, # lower bound param
            ub_param; # upper bound param
            param=lb_param .+ (ub_param .- lb_param) ./ 2,# initial guess param
            optmizator=optmizator, # selection of optimization method 
            integrator=integrator, # selection of sciml integrator
            path_to_plot=path_to_plot, # where save plots
            pt_avg=pt_avg, # numebr of the point to generate intial condition
            pt_smooth_derivative=pt_smooth_derivative,
            smoothing=smoothing, # the smoothing is done or not?
            type_of_loss=loss_type, # type of used loss 
            blank_array=blank_array, # data of all blanks
            multiple_scattering_correction=multiple_scattering_correction, # if true uses the given calibration curve to fix the data
            method_multiple_scattering_correction=method_multiple_scattering_correction,
            calibration_OD_curve=calibration_OD_curve, #  the path to calibration curve to fix the data
            PopulationSize=PopulationSize,
            maxiters=maxiters,
            abstol=abstol,
            thr_lowess=thr_lowess,
            display_plots=display_plots,# display plots in julia or not
            save_plot=save_plots,
            type_of_smoothing=type_of_smoothing,
        )


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


function fit_file_custom_ODE(
    label_exp::String, #label of the experiment
    path_to_data::String, # path to the folder to analyze
    model::Any, # string of the used model
    lb_param::Vector{Float64},# array of the array of the lower bound of the parameters
    ub_param::Vector{Float64}, # array of the array of the upper bound of the parameters
    n_equation::Int;
    path_to_annotation::Any = missing,# path to the annotation of the wells

    optmizator=BBO_adaptive_de_rand_1_bin_radiuslimited(), # selection of optimization method 
    integrator=Tsit5(), # selection of sciml integrator
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
    blank_value = 0.0,
    blank_array = [0.0],
)


    if write_res == true
        mkpath(path_to_results)
    end

    if save_plots == true
        mkpath(path_to_plot)
    end

    parameter_of_optimization = initialize_df_results_ode_custom(ub_param)

    names_of_annotated_df,properties_of_annotation,list_of_blank, list_of_discarded = reading_annotation(path_to_annotation)


    # reading files
    dfs_data = CSV.File(path_to_data)

    # shaping df for the inference
    names_of_cols = propertynames(dfs_data)

    # excluding blank data and discarded wells
    if length(list_of_blank) > 0
        names_of_cols = filter!(e -> !(e in list_of_blank), names_of_cols)
    end
    
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

        index_missing = findall(ismissing, data_values)
        index_tot =  eachindex(data_values)
        index_tot =  setdiff(index_tot,index_missing)
        data = Matrix(transpose(hcat(times_data[index_tot], data_values[index_tot])))

        # correcting negative values after blank subtraction
        data = negative_value_correction(data,
            blank_array;
            method=correct_negative,
            thr_negative=thr_negative,)




        # inference


        temp_results_1 = fitting_one_well_custom_ODE(
            data, # dataset first row times second row OD
            string(well_name), # name of the well
            label_exp, #label of the experiment
            model, # ode model to use 
            lb_param, # lower bound param
            ub_param, # upper bound param
            n_equation; # number ode in the system
            optmizator=optmizator, # selection of optimization method 
            integrator=integrator, # selection of sciml integrator
            display_plots=display_plots,# display plots in julia or not
            save_plot=save_plots,
            path_to_plot=path_to_plot, # where save plots
            pt_avg=pt_avg, # number of the point to generate intial condition
            pt_smooth_derivative=pt_smooth_derivative,
            smoothing=smoothing, # the smoothing is done or not?
            type_of_loss=loss_type, # type of used loss 
            blank_array=blank_array, # data of all blanks
            multiple_scattering_correction=multiple_scattering_correction, # if true uses the given calibration curve to fix the data
            method_multiple_scattering_correction=method_multiple_scattering_correction,
            calibration_OD_curve=calibration_OD_curve,  #  the path to calibration curve to fix the data
            PopulationSize=PopulationSize,
            maxiters=maxiters,
            abstol=abstol,
            thr_lowess=thr_lowess,
            type_of_smoothing=type_of_smoothing,
        )

        well_results = reduce(vcat, temp_results_1)

        if verbose == true
            println("the results are:")
            println(well_results)
        end

        parameter_of_optimization = hcat(parameter_of_optimization, well_results)

    end


    if write_res == true

        CSV.write(
            string(path_to_results, label_exp, "_parameters_custom_ODE_fitting.csv"),
            Tables.table(Matrix(parameter_of_optimization)),
        )


    end
    return parameter_of_optimization




end



function ODE_model_selection_file(
    label_exp::String, #label of the experiment
    path_to_data::String, # path to the folder to analyze
    models_list::Vector{String}, # ode model to use 
    lb_param_array::Any, # lower bound param
    ub_param_array::Any; # upper bound param
    path_to_annotation::Any = missing,# path to the annotation of the wells
    optmizator=BBO_adaptive_de_rand_1_bin_radiuslimited(), # selection of optimization method 
    integrator=Tsit5(), # selection of sciml integrator
    path_to_results="NA", # path where save results
    path_to_plot="NA", # path where to save Plots
    loss_type="L2", # string of the type of the used loss
    smoothing=false, # 1 do smoothing of data with rolling average
    type_of_smoothing="lowess",
    display_plot_best_model=false, # one wants the results of the best fit to be plotted
    save_plot_best_model=false,
    beta_penality=2.0, # penality for AIC evaluation
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
    correction_AIC=true,
    blank_value = 0.0,
    blank_array = [0.0],
)


    if write_res == true
        mkpath(path_to_results)
    end

    if save_plot_best_model == true
        mkpath(path_to_plot)
    end

    parameter_of_optimization = initialize_res_ms(ub_param_array)



    names_of_annotated_df,properties_of_annotation,list_of_blank, list_of_discarded = reading_annotation(path_to_annotation)

    list_of_discarded = Symbol.(list_of_discarded)

    # reading files
    dfs_data = CSV.File(path_to_data)

    # shaping df for the inference
    names_of_cols = propertynames(dfs_data)

    # excluding blank data and discarded wells


   if length(list_of_blank) > 0
        names_of_cols = filter!(e -> !(e in list_of_blank), names_of_cols)
    end
    
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

        index_missing = findall(ismissing, data_values)
        index_tot =  eachindex(data_values)
        index_tot =  setdiff(index_tot,index_missing)
        data = Matrix(transpose(hcat(times_data[index_tot], data_values[index_tot])))
        # correcting negative values after blank subtraction
        data = negative_value_correction(data,
            blank_array;
            method=correct_negative,
            thr_negative=thr_negative,)



        # inference


        temp_results_1 = ODE_Model_selection(
            data, # dataset first row times second row OD
            string(well_name), # name of the well
            label_exp, #label of the experiment
            models_list, # ode model to use 
            lb_param_array, # lower bound param
            ub_param_array; # upper bound param
            optmizator=optmizator, # selection of optimization method 
            integrator=integrator, # selection of sciml integrator
            pt_avg=pt_avg, # number of the point to generate intial condition
            beta_penality=beta_penality, # penality for AIC evaluation
            smoothing=smoothing, # the smoothing is done or not?
            type_of_smoothing=type_of_smoothing,
            thr_lowess=thr_lowess,
            type_of_loss=loss_type, # type of used loss 
            blank_array=blank_array, # data of all blanks
            display_plot_best_model=display_plot_best_model, # one wants the results of the best fit to be plotted
            save_plot_best_model=save_plot_best_model,
            path_to_plot=path_to_plot,
            pt_smooth_derivative=pt_smooth_derivative,
            multiple_scattering_correction=multiple_scattering_correction, # if true uses the given calibration curve to fix the data
            method_multiple_scattering_correction=method_multiple_scattering_correction,
            calibration_OD_curve=calibration_OD_curve, #  the path to calibration curve to fix the data
            verbose=verbose,
            PopulationSize=PopulationSize,
            maxiters=maxiters,
            abstol=abstol,
            correction_AIC=correction_AIC
        )

        vectorized_temp_results = expand_res(
            temp_results_1[end],
            lb_param_array,
            string(well_name),
            label_exp
        )
        if verbose == true
            println("the results are:")
            println(vectorized_temp_results)
        end

        parameter_of_optimization = hcat(parameter_of_optimization, vectorized_temp_results)

    end


    if write_res == true

        CSV.write(
            string(path_to_results, label_exp, "_parameters_model_selection.csv"),
            Tables.table(Matrix(parameter_of_optimization)),
        )


    end
    return parameter_of_optimization




end





"""
ODE segementation fitting fixed number of cpd for a full file
"""


function selection_ODE_fixed_change_points_file(
    label_exp::String, #label of the experiment
    path_to_data::String, # path to the folder to analyze
    list_of_models::Vector{String}, # ode model to use 
    lb_param_array::Any, # lower bound param
    ub_param_array::Any,# upper bound param
    n_max_change_points::Int;
    path_to_annotation::Any = missing,# path to the annotation of the wells
    optmizator=BBO_adaptive_de_rand_1_bin_radiuslimited(), # selection of optimization method 
    integrator=Tsit5(), # selection of sciml integrator
    type_of_loss="L2", # type of used loss 
    type_of_detection="sliding_win",
    type_of_curve="original",
    do_blank_subtraction="avg_blank",
    correct_negative="thr_correction",
    thr_negative=0.01,
    pt_avg=1, # number of the point to generate intial condition
    smoothing=true, # the smoothing is done or not?
    save_plots=false, # do plots or no
    display_plots=false, # do plots or no
    path_to_plot="NA", # where save plots
    path_to_results="NA",
    win_size=7, # numebr of the point to generate intial condition
    pt_smooth_derivative=0,
    penality_parameter=2.0,
    avg_replicate=false,
    multiple_scattering_correction="false", # if true uses the given calibration curve to fix the data
    method_multiple_scattering_correction="interpolation",
    calibration_OD_curve="NA",  #  the path to calibration curve to fix the data
    write_res=false,
    method_peaks_detection="peaks_prominence",
    bins=40,
    PopulationSize=300,
    maxiters=2000000,
    abstol=0.00001,
    type_of_smoothing="lowess",
    thr_lowess=0.05,
    verbose=false,
    correction_AIC=true,
    blank_value = 0.0,
    blank_array = [0.0],
)


    if write_res == true
        mkpath(path_to_results)
    end

    if save_plots == true
        mkpath(path_to_plot)
    end

    parameter_of_optimization = initialize_res_ms(ub_param_array, number_of_segment=n_max_change_points + 1)



    names_of_annotated_df,properties_of_annotation,list_of_blank, list_of_discarded = reading_annotation(path_to_annotation)


    # reading files
    dfs_data = CSV.File(path_to_data)

    # shaping df for the inference
    names_of_cols = propertynames(dfs_data)

    # excluding blank data and discarded wells
    if length(list_of_blank) > 0
        names_of_cols = filter!(e -> !(e in list_of_blank), names_of_cols)
    end
    
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
        index_missing = findall(ismissing, data_values)
        index_tot =  eachindex(data_values)
        index_tot =  setdiff(index_tot,index_missing)
        data = Matrix(transpose(hcat(times_data[index_tot], data_values[index_tot])))
        # correcting negative values after blank subtraction
        data = negative_value_correction(data,
            blank_array;
            method=correct_negative,
            thr_negative=thr_negative,)



        # inference

        temp_results_1 = selection_ODE_fixed_change_points(
            data, # dataset first row times second row OD
            string(well_name), # name of the well
            label_exp, #label of the experiment
            list_of_models, # ode models to use 
            lb_param_array, # lower bound param
            ub_param_array, # upper bound param
            n_max_change_points;
            type_of_loss=type_of_loss, # type of used loss 
            optmizator=optmizator, # selection of optimization method 
            integrator=integrator, # selection of sciml integrator
            type_of_detection=type_of_detection,
            type_of_curve=type_of_curve,
            smoothing=smoothing,
            pt_avg=pt_avg,
            save_plot=save_plots, # do plots or no
            display_plots=display_plots, # do plots or no
            path_to_plot=path_to_plot, # where save plots
            win_size=win_size, # numebr of the point to generate intial condition
            pt_smooth_derivative=pt_smooth_derivative,
            multiple_scattering_correction=multiple_scattering_correction, # if true uses the given calibration curve to fix the data
            method_multiple_scattering_correction=method_multiple_scattering_correction,
            calibration_OD_curve=calibration_OD_curve, #  the path to calibration curve to fix the data
            beta_smoothing_ms=penality_parameter, #  parameter of the AIC penality
            method_peaks_detection=method_peaks_detection,
            n_bins=bins,
            type_of_smoothing=type_of_smoothing,
            thr_lowess=thr_lowess,
            PopulationSize=PopulationSize,
            maxiters=maxiters,
            abstol=abstol,
            correction_AIC=correction_AIC)

        vectorized_temp_results = expand_res(
            temp_results_1[1],
            lb_param_array,
            string(well_name),
            label_exp;
            number_of_segment=n_max_change_points + 1
        )
        if verbose == true
            println("the results are:")
            println(vectorized_temp_results)
        end

        parameter_of_optimization = hcat(parameter_of_optimization, vectorized_temp_results)

    end


    if write_res == true

        CSV.write(
            string(
                path_to_results,
                label_exp,
                "_parameters_model_cpd_nseg_",
                n_max_change_points + 1,
                ".csv",
            ),
            Tables.table(Matrix(parameter_of_optimization)),
        )


    end
    return parameter_of_optimization




end




function segmentation_ODE_file(
    label_exp::String, #label of the experiment
    path_to_data::String, # path to the folder to analyze
    list_of_models::Vector{String}, # ode model to use 
    lb_param_array::Any, # lower bound param
    ub_param_array::Any,# upper bound param
    n_max_change_points::Int;
    path_to_annotation::Any = missing,# path to the annotation of the wells
    detect_number_cpd=true,
    fixed_cpd=false,
    optmizator=BBO_adaptive_de_rand_1_bin_radiuslimited(), # selection of optimization method 
    integrator=Tsit5(), # selection of sciml integrator
    type_of_loss="L2", # type of used loss 
    type_of_detection="sliding_win",
    type_of_curve="original",
    do_blank_subtraction="avg_blank",
    correct_negative="thr_correction",
    thr_negative=0.01,
    pt_avg=1, # number of the point to generate intial condition
    smoothing=true, # the smoothing is done or not?
    save_plots=false, # do plots or no
    display_plots=false, # do plots or no
    path_to_plot="NA", # where save plots
    path_to_results="NA",
    win_size=7, # numebr of the point to generate intial condition
    pt_smooth_derivative=0,
    penality_parameter=2.0,
    avg_replicate=false,
    multiple_scattering_correction="false", # if true uses the given calibration curve to fix the data
    method_multiple_scattering_correction="interpolation",
    calibration_OD_curve="NA",  #  the path to calibration curve to fix the data
    write_res=false,
    save_all_model=false,
    method_peaks_detection="peaks_prominence",
    n_bins=40,
    PopulationSize=300,
    maxiters=2000000,
    abstol=0.00001,
    type_of_smoothing="lowess",
    thr_lowess=0.05,
    verbose=false,
    correction_AIC=true,
    blank_value = 0.0,
    blank_array = [0.0],)


    vector_AIC = ["AIC"]

    if write_res == true
        mkpath(path_to_results)
    end

    if save_plots == true
        mkpath(path_to_plot)
    end

    parameter_of_optimization = initialize_res_ms(ub_param_array, number_of_segment=n_max_change_points + 1)


    names_of_annotated_df,properties_of_annotation,list_of_blank, list_of_discarded = reading_annotation(path_to_annotation)


    # reading files
    dfs_data = CSV.File(path_to_data)

    # shaping df for the inference
    names_of_cols = propertynames(dfs_data)

    # excluding blank data and discarded wells
    if length(list_of_blank) > 0
        names_of_cols = filter!(e -> !(e in list_of_blank), names_of_cols)
    end
    
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

        index_missing = findall(ismissing, data_values)
        index_tot =  eachindex(data_values)
        index_tot =  setdiff(index_tot,index_missing)
        data = Matrix(transpose(hcat(times_data[index_tot], data_values[index_tot])))

        # correcting negative values after blank subtraction
        data = negative_value_correction(data,
            blank_array;
            method=correct_negative,
            thr_negative=thr_negative,)



        # inference

        temp_results_1 = segmentation_ODE(
            data, # dataset x times y OD/fluorescence
            string(well_name), # name of the well
            label_exp, #label of the experiment
            list_of_models, # ode model to use
            lb_param_array, # lower bound param
            ub_param_array, # upper bound param
            n_max_change_points;
            detect_number_cpd=detect_number_cpd,
            fixed_cpd=fixed_cpd,
            optmizator=optmizator, # selection of optimization method
            integrator=integrator, # selection of sciml integrator
            type_of_loss=type_of_loss, # type of used loss
            type_of_detection=type_of_detection,
            type_of_curve=type_of_curve,
            pt_avg=pt_avg, # number of the point to generate intial condition
            smoothing=smoothing, # the smoothing is done or not?
            save_plot=save_plots, # do plots or no
            display_plot=display_plots, # do plots or no
            path_to_plot=path_to_plot, # where save plots
            path_to_results=path_to_results,
            win_size=win_size, # numebr of the point to generate intial condition
            pt_smooth_derivative=pt_smooth_derivative,
            penality_parameter=penality_parameter,
            multiple_scattering_correction=multiple_scattering_correction, # if true uses the given calibration curve to fix the data
            method_multiple_scattering_correction=method_multiple_scattering_correction,
            calibration_OD_curve=calibration_OD_curve,  #  the path to calibration curve to fix the data
            save_all_model=save_all_model,
            method_peaks_detection=method_peaks_detection,
            n_bins=n_bins,
            PopulationSize=PopulationSize,
            maxiters=maxiters,
            abstol=abstol,
            type_of_smoothing=type_of_smoothing,
            thr_lowess=thr_lowess,
            correction_AIC=correction_AIC,
            )
            vector_AIC =  hcat( vector_AIC,temp_results_1[end])

        vectorized_temp_results = expand_res(
            temp_results_1[1],
            lb_param_array,
            string(well_name),
            label_exp;
            number_of_segment=length(temp_results_1[1])
        )
        if verbose == true
            println("the results are:")
            println(vectorized_temp_results)
        end

        parameter_of_optimization = hcat(parameter_of_optimization, vectorized_temp_results)

    end


    if write_res == true

        CSV.write(
            string(
                path_to_results,
                label_exp,
                "_parameters_model_cpd_nseg_",
                n_max_change_points + 1,
                ".csv",
            ),
            Tables.table(Matrix(parameter_of_optimization)),
        )


    end
    return parameter_of_optimization,vector_AIC




end

export fit_file_ODE
export fit_file_custom_ODE
export ODE_model_selection_file
export selection_ODE_fixed_change_points_file
export segmentation_ODE_file
