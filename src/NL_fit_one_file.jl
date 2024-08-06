using OptimizationMultistartOptimization

"""
    fit_NL_model_file(
    label_exp::String,
    path_to_data::String,
    model::Any,
    u0;
    lb_param::Vector{Float64}=nothing,
    ub_param::Vector{Float64}=nothing,
    path_to_annotation::Any=missing,
    method_of_fitting="NA",
    nrep=100,
    errors_estimation=false,
    path_to_results="NA",
    loss_type="RE",
    smoothing=false,
    type_of_smoothing="lowess",
    verbose=false,
    write_res=false,
    pt_avg=1,
    pt_smooth_derivative=7,
    do_blank_subtraction="avg_blank",
    avg_replicate=false,
    correct_negative="remove",
    thr_negative=0.01,
    multiple_scattering_correction=false,
    method_multiple_scattering_correction="interpolation",
    calibration_OD_curve="NA",
    thr_lowess=0.05,
    penality_CI=8.0,
    size_bootstrap=0.7,
    blank_value=0.0,
    blank_array=[0.0],
    optimizer=BBO_adaptive_de_rand_1_bin_radiuslimited(),
    multistart=false,
    n_restart=50,
    auto_diff_method=nothing,
    cons=nothing,
    opt_params...
    )

This function performs nonlinear (NL) model fitting for a given CSV file using specified models and optimization techniques.

# Arguments:

- `label_exp::String`: The label of the experiment.

- `path_to_data::String`: Path to the CSV file containing the data. The file should be formatted with time in the first column and the corresponding data values in the second column.

- `model::Any`: The nonlinear model to fit. This can be a function or a string representing a hardcoded NL model.

- `u0`: Initial guess for the model parameters.

# Key Arguments:

- `lb_param::Vector{Float64} = nothing`: Vector o

"""
function fit_NL_model_file(
    label_exp::String, #label of the experiment
    path_to_data::String, # path to the folder to analyze
    model::Any, # string of the used model
    u0;
    lb_param::Vector{Float64}=nothing,# array of the array of the lower bound of the parameters
    ub_param::Vector{Float64}=nothing, # array of the array of the upper bound of the parameters
    path_to_annotation::Any = missing,# path to the annotation of the wells
    method_of_fitting="Normal",
    nrep=10,
    errors_estimation=false,
    path_to_results="NA", # path where save results
    loss_type="RE", # string of the type of the used loss
    smoothing=false, # 1 do smoothing of data with rolling average
    type_of_smoothing="lowess",
    verbose=false, # 1 true verbose
    write_res=false, # write results
    pt_avg=1, # number of points to do smoothing average
    pt_smooth_derivative=7, # number of points to do ssmooth_derivative
    do_blank_subtraction="avg_blank", # string on how to use blank (NO,avg_subtraction,time_avg)
    avg_replicate=false, # if true the average between replicates is fitted. If false all replicate are fitted indipendelitly
    correct_negative="remove", # if "thr_correction" it put a thr on the minimum value of the data with blank subracted, if "blank_correction" uses blank distrib to impute negative values
    thr_negative=0.01,  # used only if correct_negative == "thr_correction"
    multiple_scattering_correction=false, # if true uses the given calibration curve to fix the data
    method_multiple_scattering_correction="interpolation",
    calibration_OD_curve="NA",  #  the path to calibration curve to fix the data
    thr_lowess=0.05,
    penality_CI=8.0,
    size_bootstrap=0.7,
    blank_value = 0.0,
    blank_array = [0.0],
    optimizer=BBO_adaptive_de_rand_1_bin_radiuslimited(),
    multistart=false,
    n_restart=50,
    auto_diff_method=nothing,
    cons=nothing,
    opt_params...
)


    if write_res == true
        mkpath(path_to_results)
    end

    fits= ()
    data_to_save = ()
    parameter_of_optimization = initialize_df_results_ode_custom(lb_param)
    errors_of_optimization = initialize_df_results_ode_custom(lb_param)
  
    if typeof(model) == String

        model_string = NL_models[model].name


    else

        model_string = "custom"


    end



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

        data = Matrix(data)

        # defining time steps of the inference
        if method_of_fitting == "Bootstrap"



            temp_results_1 = fit_NL_model_bootstrap(data, # dataset first row times second row OD
                string(well_name), # name of the well
                label_exp, #label of the experiment
                model, #  model to use
                u0;
                lb_param=lb_param, # lower bound param
                ub_param=ub_param, # upper bound param
                nrep=nrep,
                optimizer=optimizer,
                size_bootstrap=size_bootstrap,
                pt_avg=pt_avg, # numebr of the point to generate intial condition
                pt_smooth_derivative=pt_smooth_derivative,
                smoothing=smoothing, # the smoothing is done or not?
                type_of_smoothing=type_of_smoothing,
                type_of_loss=loss_type, # type of used loss
                multiple_scattering_correction=multiple_scattering_correction, # if true uses the given calibration curve to fix the data
                method_multiple_scattering_correction=method_multiple_scattering_correction,
                calibration_OD_curve=calibration_OD_curve,  #  the path to calibration curve to fix the data
                thr_lowess=thr_lowess,
                write_res=write_res,
                penality_CI=penality_CI,
                path_to_results = path_to_results,
                auto_diff_method=auto_diff_method,
                multistart=multistart,
                n_restart=n_restart,
                cons=cons,
                 opt_params...
                )

            temp_mean = temp_results_1[7]
            temp_mean = vcat("mean",temp_mean)
            temp_mean = vcat(string(well_name),temp_mean)

            errors_of_optimization = hcat(errors_of_optimization, temp_mean)

            temp_ci_low = temp_results_1[9]
            temp_ci_low = vcat("lower_CI",temp_ci_low)
            temp_ci_low = vcat(string(well_name),temp_ci_low)

            errors_of_optimization = hcat(errors_of_optimization, temp_ci_low)
            
            temp_ci_up = temp_results_1[10]
            temp_ci_up = vcat("upper_CI",temp_ci_up)
            temp_ci_up = vcat(string(well_name),temp_ci_up)


            errors_of_optimization = hcat(errors_of_optimization, temp_ci_up)



        elseif method_of_fitting == "Morris_sensitivity"


            temp_results_1 = fit_NL_model_with_sensitivity(data, # dataset first row times second row OD
                string(well_name), # name of the well
                label_exp, #label of the experiment
                model, #  model to use
                lb_param, # lower bound param
                ub_param; # upper bound param
                nrep=nrep,
                optimizer=optimizer,
                pt_avg=pt_avg, # number of the point to generate intial condition
                pt_smooth_derivative=pt_smooth_derivative,
                smoothing=smoothing, # the smoothing is done or not?
                type_of_smoothing=type_of_smoothing,
                type_of_loss=type_of_loss, # type of used loss
                multiple_scattering_correction=multiple_scattering_correction, # if true uses the given calibration curve to fix the data
                method_multiple_scattering_correction=method_multiple_scattering_correction,
                calibration_OD_curve=calibration_OD_curve,  #  the path to calibration curve to fix the data
                thr_lowess=thr_lowess,
                write_res=write_res,
                penality_CI=penality_CI,
                auto_diff_method=auto_diff_method,
                multistart=multistart,
                n_restart=n_restart,
                cons=cons,
                opt_params...
            )

 


        else



            temp_results_1 = fit_NL_model_bootstrap(data, # dataset first row times second row OD
            name_well, # name of the well
            label_exp, #label of the experiment
            model_to_test, #  model to use
            u0;
            lb_param=temp_param_lb, # lower bound param
            ub_param=temp_param_ub, # upper bound param
            nrep=nrep,
            optimizer=optimizer,
            size_bootstrap=size_bootstrap,
            pt_avg=pt_avg, # number of the point to generate intial condition
            pt_smooth_derivative=pt_smooth_derivative,
            smoothing=smoothing, # the smoothing is done or not?
            type_of_smoothing=type_of_smoothing,
            type_of_loss=type_of_loss, # type of used loss
            multiple_scattering_correction=multiple_scattering_correction, # if true uses the given calibration curve to fix the data
            method_multiple_scattering_correction=method_multiple_scattering_correction,
            calibration_OD_curve=calibration_OD_curve,  #  the path to calibration curve to fix the data
            thr_lowess=thr_lowess,
            write_res=write_res,
            penality_CI=penality_CI,
            auto_diff_method=auto_diff_method,
            multistart=multistart,
            n_restart=n_restart,
            cons=cons,
             opt_params...
        )





        end


        data = Matrix(data)




        if verbose == true
            println("the results are:")
            println(temp_results_1[2])
        end
        temp_array =  temp_results_1[2]

        parameter_of_optimization = hcat(parameter_of_optimization, temp_array)

        if errors_estimation == true && method_of_fitting != "Bootstrap"


            best_param = temp_results_1[2][3:(end-3)]
            best_param = convert.(Float64, best_param)
            temp_errors_of_optimization = NL_error_blanks(data, # dataset first row times second row OD
                string(well_name), # name of the well
                label_exp, #label of the experiment
                model, #  model to use
                u0,
                blank_array; # upper bound param
                nrep=nrep,
                u0=best_param,# initial guess param
                optimizer=optimizer,
                pt_avg=pt_avg, # numebr of the point to generate intial condition
                pt_smooth_derivative=pt_smooth_derivative,
                smoothing=smoothing, # the smoothing is done or not?
                type_of_smoothing=type_of_smoothing,
                type_of_loss=loss_type, # type of used loss
                multiple_scattering_correction=multiple_scattering_correction, # if true uses the given calibration curve to fix the data
                method_multiple_scattering_correction=method_multiple_scattering_correction,
                calibration_OD_curve=calibration_OD_curve,  #  the path to calibration curve to fix the data
                thr_lowess=thr_lowess,
                penality_CI=penality_CI,
                multistart=multistart,
                n_restart=n_restart,
                auto_diff_method=auto_diff_method,
                cons=cons,
                opt_params...
            )




            temp_mean = temp_errors_of_optimization[5]
            temp_mean = vcat("mean",temp_mean)
            temp_mean = vcat(string(well_name),temp_mean)

            errors_of_optimization = hcat(errors_of_optimization, temp_mean)

            temp_ci_low = temp_errors_of_optimization[7]
            temp_ci_low = vcat("lower_CI",temp_ci_low)
            temp_ci_low = vcat(string(well_name),temp_ci_low)

            errors_of_optimization = hcat(errors_of_optimization, temp_ci_low)
            
            temp_ci_up = temp_errors_of_optimization[8]
            temp_ci_up = vcat("upper_CI",temp_ci_up)
            temp_ci_up = vcat(string(well_name),temp_ci_up)


            errors_of_optimization = hcat(errors_of_optimization, temp_ci_up)



        end
        fits = (fits...,temp_results_1[3] )
        if multiple_scattering_correction == true

            data = correction_OD_multiple_scattering(data, calibration_OD_curve; method=method_multiple_scattering_correction)
        end
        data_to_save = (data_to_save...,data)

    end


    if write_res == true

        CSV.write(
            string(path_to_results, label_exp, "_parameters_", model_string, ".csv"),
            Tables.table(Matrix(parameter_of_optimization)),
        )


    end

    KinBiont_res_one_file = ("NL",parameter_of_optimization,fits,data_to_save)


    return KinBiont_res_one_file, errors_of_optimization




end





"""
    fit_NL_model_selection_file(
    label_exp::String,
    path_to_data::String,
    list_model_function::Any,
    list_u0;
    lb_param_array::Any = nothing,
    ub_param_array::Any = nothing,
    path_to_annotation::Any = missing,
    method_of_fitting="Normal",
    nrep=10,
    path_to_results="NA",
    loss_type="RE",
    smoothing=false,
    type_of_smoothing="lowess",
    verbose=false,
    write_res=false,
    pt_avg=1,
    pt_smooth_derivative=7,
    do_blank_subtraction="avg_blank",
    avg_replicate=false,
    correct_negative="remove",
    thr_negative=0.01,
    multiple_scattering_correction=false,
    method_multiple_scattering_correction="interpolation",
    calibration_OD_curve="NA",
    thr_lowess=0.05,
    beta_smoothing_ms=2.0,
    penality_CI=8.0,
    size_bootstrap=0.7,
    correction_AIC=true,
    blank_value = 0.0,
    blank_array = [0.0],
    optimizer=BBO_adaptive_de_rand_1_bin_radiuslimited(),
    multistart=false,
    n_restart=50,
    auto_diff_method=nothing,
    cons=nothing,
    opt_params...
    )

This function performs nonlinear (NL) model selection for an array of NL models using either AIC or AICc. It applies to a full CSV file of data.

# Arguments:

- `label_exp::String`: The label for the experiment.

- `path_to_data::String`: Path to the CSV file containing the data. The file should be formatted with time in the first column and corresponding data values in subsequent columns.

- `list_model_function::Any`: Array of functions or strings representing the NL models to be tested.

- `list_u0::Vector{Float64}`: Initial guesses for the parameters of each NL model.

# Key Arguments:

- `lb_param_array::Any = nothing`: Array of lower bounds for the parameters of each model. If not provided, no bounds are enforced.

- `ub_param_array::Any = nothing`: Array of upper bounds fo

"""
function fit_NL_model_selection_file(
    label_exp::String, #label of the experiment
    path_to_data::String, # path to the folder to analyze
    list_model_function::Any, #  model to use
    list_u0;# initial guess param
    lb_param_array::Any = nothing, # lower bound param
    ub_param_array::Any = nothing, # upper bound param
    path_to_annotation::Any = missing,# path to the annotation of the wells
    method_of_fitting="Normal",
    nrep=10,
    path_to_results="NA", # path where save results
    loss_type="RE", # string of the type of the used loss
    smoothing=false, # 1 do smoothing of data with rolling average
    type_of_smoothing="lowess",
    verbose=false, # 1 true verbose
    write_res=false, # write results
    pt_avg=1, # number of points to do smoothing average
    pt_smooth_derivative=7, # number of points to do ssmooth_derivative
    do_blank_subtraction="avg_blank", # string on how to use blank (NO,avg_subtraction,time_avg)
    avg_replicate=false, # if true the average between replicates is fitted. If false all replicate are fitted indipendelitly
    correct_negative="remove", # if "thr_correction" it put a thr on the minimum value of the data with blank subracted, if "blank_correction" uses blank distrib to impute negative values
    thr_negative=0.01,  # used only if correct_negative == "thr_correction"
    multiple_scattering_correction=false, # if true uses the given calibration curve to fix the data
    method_multiple_scattering_correction="interpolation",
    calibration_OD_curve="NA",  #  the path to calibration curve to fix the data
    thr_lowess=0.05,
    beta_smoothing_ms=2.0,
    penality_CI=8.0,
    size_bootstrap=0.7,
    correction_AIC=true,
    blank_value = 0.0,
    blank_array = [0.0],
    optimizer=BBO_adaptive_de_rand_1_bin_radiuslimited(),
    multistart=false,
    n_restart=50,
    auto_diff_method=nothing,
    cons=nothing,
    opt_params...
)


    if write_res == true
        mkpath(path_to_results)
    end

    parameter_of_optimization = initialize_res_ms(list_u0)
    names_of_annotated_df,properties_of_annotation,list_of_blank, list_of_discarded = reading_annotation(path_to_annotation)

    fits= ()
    data_to_save = ()
    
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

        data = Matrix(data)

        # defining time steps of the inference


        temp_results_1 = NL_model_selection(data, # dataset first row times second row OD
            string(well_name), # name of the well
            label_exp, #label of the experiment
            list_model_function, #  model to use
            list_u0;
            lb_param_array=lb_param_array, # lower bound param
            ub_param_array=ub_param_array, # upper bound param
            method_of_fitting=method_of_fitting,
            nrep=nrep,
            size_bootstrap=size_bootstrap,
            pt_avg=pt_avg, # numebr of the point to generate intial condition
            pt_smooth_derivative=pt_smooth_derivative,
            smoothing=smoothing, # the smoothing is done or not?
            type_of_smoothing=type_of_smoothing,
            type_of_loss=loss_type, # type of used loss
            multiple_scattering_correction=multiple_scattering_correction, # if true uses the given calibration curve to fix the data
            method_multiple_scattering_correction=method_multiple_scattering_correction,
            calibration_OD_curve=calibration_OD_curve,  #  the path to calibration curve to fix the data
            thr_lowess=thr_lowess,
            write_res=false,
            beta_smoothing_ms=beta_smoothing_ms,
            penality_CI=penality_CI,
            correction_AIC=correction_AIC,
            optimizer=optimizer,
            auto_diff_method=auto_diff_method,
            multistart=multistart,
            n_restart=n_restart,
            cons=cons,
            opt_params...)



        data = Matrix(data)




        if verbose == true
            println("the results are:")
            println(temp_results_1[2])
        end
        parameter_of_optimization = hcat(parameter_of_optimization, temp_results_1[2])
        fits = (fits...,temp_results_1[3] )
        if multiple_scattering_correction == true

            data = correction_OD_multiple_scattering(data, calibration_OD_curve; method=method_multiple_scattering_correction)
        end
        data_to_save = (data_to_save...,data)
    end


    if write_res == true

        CSV.write(
            string(path_to_results, label_exp, "_parameters_nl.csv"),
            Tables.table(Matrix(parameter_of_optimization)),
        )


    end
    KinBiont_res_one_file = ("NL_model_selection",parameter_of_optimization,fits,data_to_save)

    return KinBiont_res_one_file




end



"""
    fit_NL_segmentation_file(
    label_exp::String,
    path_to_data::String,
    list_model_function::Any,
    list_u0,
    n_change_points::Int;
    lb_param_array::Vector{Vector{Float64}}=nothing,
    ub_param_array::Vector{Vector{Float64}}=nothing,
    path_to_annotation::Any = missing,
    method_of_fitting="NA",
    nrep=10,
    path_to_results="NA",
    loss_type="RE",
    smoothing=false,
    type_of_smoothing="lowess",
    verbose=false,
    write_res=false,
    pt_avg=1,
    pt_smooth_derivative=7,
    do_blank_subtraction="avg_blank",
    avg_replicate=false,
    correct_negative="remove",
    thr_negative=0.01,
    multiple_scattering_correction=false,
    method_multiple_scattering_correction="interpolation",
    calibration_OD_curve="NA",
    size_bootstrap=0.7,
    thr_lowess=0.05,
    detect_number_cpd=true,
    type_of_detection="sliding_win",
    type_of_curve="original",
    fixed_cpd=false,
    penality_CI=8.0,
    beta_smoothing_ms=2.0,
    win_size=7,
    n_bins=40,
    correction_AIC=true,
    blank_value = 0.0,
    blank_array = [0.0],
    optimizer=BBO_adaptive_de_rand_1_bin_radiuslimited(),
    auto_diff_method=nothing,
    multistart=false,
    n_restart=50,
    cons=nothing,
    opt_params...
   )

This function performs nonlinear (NL) model selection on a segmented time series using AIC or AICc. It operates on an entire CSV file of data.

# Arguments:

- `label_exp::String`: The label for the experiment.

- `path_to_data::String`: Path to the CSV file containing the data. The file should be formatted with time in the first column and corresponding data values in subsequent columns.

- `list_model_function::Any`: Array of functions or strings representing the NL models to be tested.

- `list_u0::Any`: Initial guesses for the parameters of each NL model.

- `n_change_points::Int`: Maximum number of change points to consider in the segmentation process.

# Key Arguments:

- `lb_param_array::Vector{Vector{Float64}} = nothing`: Array of lower bounds for the parameters of each model. Each entry corresponds to the parameter bounds for a specific model.

- `ub_param_array::Vector{Vector{Float64}} = nothing`: Array of upper bounds for the parameters of each model. Each entry corresponds to the parameter bounds for a specific model.

- `path_to_annotation::Any = missing`: Path to a CSV file with annotation data, if available.

- `method_of_fitting::String = "NA"`: Method for NL fitting. Options include `"Bootstrap"`, `"Normal"`, and `"Morris_sensitivity"`. Default is `"NA"`.

- `nrep::Int = 10`: Number of repetitions for methods like Morris sensitivity or bootstrap. Used only if `method_of_fitting` is `"Bootstrap"` or `"Morris_sensitivity"`.

- `path_to_results::String = "NA"`: Path to the folder where results will be saved.

- `loss_type::String = "RE"`: Type of loss function used. Options include `"RE"` (Residual Error), `"L2"`, `"L2_derivative"`, and `"blank_weighted_L2"`.

- `smoothing::Bool = false`: Whether to apply smoothing to the data.

- `type_of_smoothing::String = "lowess"`: Type of smoothing. Options include `"NO"`, `"rolling_avg"`, and `"lowess"`.

- `pt_avg::Int = 1`: Number of points for rolling average smoothing or initial condition generation.

- `pt_smooth_derivative::Int = 7`: Number of points for evaluating the specific growth rate. If less than 2, uses an interpolation algorithm.

- `do_blank_subtraction::String = "avg_blank"`: Method for blank subtraction. Options include `"NO"`, `"avg_subtraction"`, and `"time_avg"`.

- `avg_replicate::Bool = false`: If true, averages replicates if annotation data is provided.

- `correct_negative::String = "remove"`: Method for treating negative values after blank subtraction. Options include `"thr_correction"`, `"blank_correction"`, and `"remove"`.

- `thr_negative::Float64 = 0.01`: Threshold for correcting negative values if `correct_negative` is `"thr_correction"`.

- `multiple_scattering_correction::Bool = false`: If true, uses a calibration curve to correct data for multiple scattering.

- `method_multiple_scattering_correction::String = "interpolation"`: Method for correcting multiple scattering. Options include `"interpolation"` and `"exp_fit"`.

- `calibration_OD_curve::String = "NA"`: Path to the CSV file containing calibration data, used if `multiple_scattering_correction` is true.

- `thr_lowess::Float64 = 0.05`: Parameter for lowess smoothing.

- `beta_smoothing_ms::Float64 = 2.0`: Penalty parameter for AIC (or AICc) evaluation.

- `penality_CI::Float64 = 8.0`: Penalty parameter for ensuring continuity in segmentation.

- `size_bootstrap::Float64 = 0.7`: Fraction of data used for each bootstrap run, used only if `method_of_fitting` is `"Bootstrap"`.

- `correction_AIC::Bool = true`: If true, performs finite sample correction for AIC.

- `blank_value::Float64 = 0.0`: Average value of blanks used if `do_blank_subtraction` is not `"NO"` and `path_to_annotation` is missing.

- `blank_array::Vector{Float64} = [0.0]`: Array of blank values used if `do_blank_subtraction` is not `"NO"` and `path_to_annotation` is missing.

- `type_of_detection::String = "sliding_win"`: Method for detecting change points. Options include `"sliding_win"` for sliding window approach and `"lsdd"` for least square density difference (LSDD) from `ChangePointDetection.jl`.

- `type_of_curve::String = "original"`: Specifies the curve on which change point detection is performed. Options include `"original"` (the original time series) and `"deriv"` (the specific growth rate time series).

- `fixed_cpd::Bool = false`: If true, returns the fitting using exactly `n_change_points` change points.

- `detect_number_cpd::Bool = true`: If true, all possible combinations of lengths 1, 2, ..., `n_change_points` are tested, and the best combination for AICc is returned.

- `method_peaks_detection::String = "peaks_prominence"`: Method for peak detection on the dissimilarity curve. Options include `"peaks_prominence"` (orders peaks by prominence) and `"thr_scan"` (uses a threshold to select peaks).

- `n_bins::Int = 40`: Number of bins used to generate the threshold for peak detection if `method_peaks_detection` is `"thr_scan"`.

- `win_size::Int = 14`: Size of the window used by the change point detection algorithms.

- `auto_diff_method::Any = nothing`: Differentiation method to be specified if required by the optimizer.

- `cons::Any = nothing`: Constraints for optimization equations.

- `multistart::Bool = false`: If true, performs multistart optimization.

- `n_restart::Int = 50`: Number of restarts for multistart optimization if `multistart` is true.

- `optimizer::Any = BBO_adaptive_de_rand_1_bin_radiuslimited()`: Optimizer used for the fitting process. Default is `BBO_adaptive_de_rand_1_bin_radiuslimited()`.

- `opt_params...`: Additional optional parameters for the optimizer.

# Outputs:

- **Method String**: Description of the fitting method used.

- **Results Matrix**: A matrix where each row contains:
  - `"name of model"`: The name of the model used.
  - `"well"`: The name of the well.
  - `"param_1", "param_2", ..., "param_n"`: Parameters of the selected NL model.
  - `"maximum specific gr using NL"`: Maximum specific growth rate obtained using the NL model.
  - `"maximum specific gr using data"`: Maximum specific growth rate obtained from the data.
  - `"objective function value"`: The value of the objective function (i.e., loss of the solution).

- **Fittings**: Details of the model fittings.

- **Preprocessed Data**: Data that has been preprocessed according to the specified arguments.

- **Change Points Time**: Time of change points detected for each well.

- **AIC (or AICc) Score**: AIC or AICc score for each well.

"""
function fit_NL_segmentation_file(
    label_exp::String, #label of the experiment
    path_to_data::String, # path to the folder to analyze
    list_model_function::Any, #  model to use
    list_u0,# initial guess param
    n_change_points::Int;
    lb_param_array::Vector{Vector{Float64}}=nothing, # lower bound param
    ub_param_array::Vector{Vector{Float64}}=nothing, # upper bound param
    path_to_annotation::Any = missing,# path to the annotation of the wells
    method_of_fitting="NA",
    nrep=10,
    path_to_results="NA", # path where save results
    loss_type="RE", # string of the type of the used loss
    smoothing=false, # 1 do smoothing of data with rolling average
    type_of_smoothing="lowess",
    verbose=false, # 1 true verbose
    write_res=false, # write results
    pt_avg=1, # number of points to do smoothing average
    pt_smooth_derivative=7, # number of points to do ssmooth_derivative
    do_blank_subtraction="avg_blank", # string on how to use blank (NO,avg_subtraction,time_avg)
    avg_replicate=false, # if true the average between replicates is fitted. If false all replicate are fitted indipendelitly
    correct_negative="remove", # if "thr_correction" it put a thr on the minimum value of the data with blank subracted, if "blank_correction" uses blank distrib to impute negative values
    thr_negative=0.01,  # used only if correct_negative == "thr_correction"
    multiple_scattering_correction=false, # if true uses the given calibration curve to fix the data
    method_multiple_scattering_correction="interpolation",
    calibration_OD_curve="NA",  #  the path to calibration curve to fix the data
    size_bootstrap=0.7,
    thr_lowess=0.05,
    detect_number_cpd=true,
    type_of_detection="sliding_win",
    type_of_curve="original",
    fixed_cpd=false,
    penality_CI=8.0,
    beta_smoothing_ms=2.0,
    win_size=7, # number of the point of cpd sliding win
    n_bins=40,
    correction_AIC=true,
    blank_value = 0.0,
    blank_array = [0.0],
    optimizer=BBO_adaptive_de_rand_1_bin_radiuslimited(),
    auto_diff_method=nothing,
    multistart=false,
    n_restart=50,
    cons=nothing,
    opt_params...
)


    if write_res == true
        mkpath(path_to_results)
    end

    parameter_of_optimization = initialize_res_ms(list_u0, number_of_segment=n_change_points)
    names_of_annotated_df,properties_of_annotation,list_of_blank, list_of_discarded = reading_annotation(path_to_annotation)
    fits= ()
    data_to_save = ()
    cps = ()

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

        data = Matrix(data)

        # defining time steps of the inference


        temp_results_1 = segmentation_NL(
            data, # dataset first row times second row OD
            string(well_name), # name of the well
            label_exp, #label of the experiment
            list_model_function, # ode models to use
            list_u0,
            n_change_points;
            lb_param_array=lb_param_array, # lower bound param
            ub_param_array=ub_param_array, # upper bound param
            type_of_loss=loss_type, # type of used loss
            method_of_fitting=method_of_fitting, # selection of sciml integrator
            type_of_detection=type_of_detection,
            type_of_curve=type_of_curve,
            smoothing=smoothing,
            nrep=nrep,
            type_of_smoothing=type_of_smoothing,
            thr_lowess=thr_lowess,
            pt_avg=pt_avg,
            win_size=win_size, # numebr of the point to generate intial condition
            pt_smooth_derivative=pt_smooth_derivative,
            multiple_scattering_correction=multiple_scattering_correction, # if true uses the given calibration curve to fix the data
            method_multiple_scattering_correction=method_multiple_scattering_correction,
            calibration_OD_curve=calibration_OD_curve, #  the path to calibration curve to fix the data
            beta_smoothing_ms=beta_smoothing_ms, #  parameter of the AIC penality
            n_bins=n_bins,
            detect_number_cpd=detect_number_cpd,
            fixed_cpd=fixed_cpd,
            penality_CI=penality_CI,
            size_bootstrap=size_bootstrap,
            correction_AIC=correction_AIC,
            optimizer=optimizer,
            multistart=multistart,
            n_restart=n_restart,
            auto_diff_method=auto_diff_method,
            cons=cons,
            opt_params...
            
        )


        data = Matrix(data)




        if verbose == true
            println("the results are:")
            println(temp_results_1[2])
        end

        results_to_bind = expand_res(
            temp_results_1[2],
            list_u0,
            string(well_name),
            label_exp;
            number_of_segment=length(temp_results_1[2]))

        parameter_of_optimization = hcat(parameter_of_optimization, results_to_bind)
        seg_fit = hcat(temp_results_1[4],temp_results_1[3])
        fits = (fits...,seg_fit )
        if multiple_scattering_correction == true

            data = correction_OD_multiple_scattering(data, calibration_OD_curve; method=method_multiple_scattering_correction)
        end
        data_to_save = (data_to_save...,data)
        cps = (cps...,temp_results_1[5])

    end


    if write_res == true

        CSV.write(
            string(path_to_results, label_exp, "_parameters_ segmentation_NL.csv"),
            Tables.table(Matrix(parameter_of_optimization)),
        )


    end

    KinBiont_res_one_file = ("NL_segmentation", parameter_of_optimization,fits,data_to_save)

    return KinBiont_res_one_file




end

export fit_NL_model_file
export fit_NL_model_selection_file
export fit_NL_segmentation_file
