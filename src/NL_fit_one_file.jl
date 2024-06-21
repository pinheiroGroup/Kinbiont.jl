using OptimizationMultistartOptimization

"""
    fit_NL_model_file(
    label_exp::String,
    path_to_data::String, 
    model::Any,
    lb_param::Vector{Float64},
    ub_param::Vector{Float64};
    path_to_annotation::Any = missing,
    u0=lb_param .+ (ub_param .- lb_param) ./ 2,
    method_of_fitting="MCMC",
    nrep=100,
    errors_estimation=false,
    optimizer=BBO_adaptive_de_rand_1_bin_radiuslimited(), 
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
    correct_negative="thr_correction",
    thr_negative=0.01,  
    multiple_scattering_correction=false, 
    method_multiple_scattering_correction="interpolation",
    calibration_OD_curve="NA",
    thr_lowess=0.05,
    penality_CI=8.0,
    size_bootstrap=0.7,
    blank_value = 0.0,
    blank_array = [0.0],
    )



This function performs NL model selection of one NL model for a full csv file
# Arguments

- `label_exp::String`,  label of the experiment.
- `path_to_data::String`. path to the csv data frame. See documentation for formatting it.
-  `model::Any`:  functions or strings (for harcoded NL model) of the NL models
-  `lb_param::Any`:Array of Lower bounds for the parameters (compatible with the models).
-  `ub_param::Any`:Array of Upper bounds for the parameters (compatible with the models).

# Key Arguments:
- `method_of_fitting="MCMC"`: String, how perform the NL fit. Options "MCMC","Bootstrap","Normal", and "Morris_sensitivity"
- `nrep=100`. Number of MCMC steps.
- `param= lb_param .+ (ub_param.-lb_param)./2`:Vector{Float64}, Initial guess for the model parameters.
- `optimizer =   BBO_adaptive_de_rand_1_bin_radiuslimited()` optimizer from optimizationBBO.
- `type_of_smoothing="rolling_avg"`: String, How to smooth the data, options: "NO" , "rolling avg" rolling average of the data, and "lowess".
- `pt_avg=7`: Number of points to generate the initial condition or do the rolling avg smoothing.
- `smoothing=false`: Whether to apply smoothing to the data or not.
- `type_of_loss:="RE" `: Type of loss function to be used. (options= "RE", "L2", "L2_derivative" and "blank_weighted_L2").
- `pt_smoothing_derivative=7`:Int,  Number of points for evaluation of specific growth rate. If <2 it uses interpolation algorithm otherwise a sliding window approach.
- `calibration_OD_curve="NA"`: String, The path where the .csv calibration data are located, used only if `multiple_scattering_correction=true`.
- `multiple_scattering_correction=false`: Bool, if true uses the given calibration curve to correct the data for muliple scattering.
- `method_multiple_scattering_correction="interpolation"`: String, How perform the inference of multiple scattering curve, options: "interpolation" or   "exp_fit" it uses an exponential fit from "Direct optical density determination of bacterial cultures in microplates for high-throughput screening applications"
-  `thr_lowess=0.05`: Float64, keyword argument of lowess smoothing
- `blank_value = 0.0`: used only if `path_to_annotation = missing`and `blank_subtraction != "NO "`. It is used as average value of the blank.
- `blank_array = [0.0]`:used only if `path_to_annotation = missing`and `blank_subtraction != "NO "`. It is used as array of the blanks values.
-  `correct_negative="thr_correction"`  ;: String, How to treat negative values after blank subtraction. If `"thr_correction"` it put a thr on the minimum value of the data with blank subracted, if `"blank_correction"` uses blank distribution to impute negative values, if `"remove"` the values are just removed..
- `do_blank_subtraction="NO"`: String, how perform the blank subtration, options "NO","avg_subtraction" (subtration of average value of blanks) and "time_avg" (subtration of  time average value of blanks).  
- `penality_CI=2.0`, used only in segementation to force the optimization to respect continuty on bonduar
-  `correction_AIC=true`: Bool, do finite samples correction of AIC.
-  `beta_smoothing_ms=2.0` penality  parameters for AIC (or AICc) evaluation.
-  `size_bootstrap=0.7`: Float, the fraction of data used each Bootstrap run. Used only if method is "Bootstrap"
- `write_res=false`: Bool, write the results in path_to_results folder.
- `path_to_results= "NA"`:String, path to the folder where save the results.
-  `correct_negative="thr_correction"`: # if "thr_correction" it put a thr on the minimum value of the data with blank subracted, if "blank_correction" uses blank distrib to impute negative values.

# Output (if `results_NL_fit =fit_NL_model_file(...)`:


- a matrix with the following contents for each row : `[ "label of exp", "well", "param_1","param_2",..,"param_n","maximum specific gr using model","maximum specific gr using data", "objective function value (i.e. loss of the solution)"]` where ' "param_1","param_2",..,"param_n" ' .

"""
function fit_NL_model_file(
    label_exp::String, #label of the experiment
    path_to_data::String, # path to the folder to analyze
    model::Any, # string of the used model
    u0;
    lb_param::Vector{Float64}=nothing,# array of the array of the lower bound of the parameters
    ub_param::Vector{Float64}=nothing, # array of the array of the upper bound of the parameters
    path_to_annotation::Any = missing,# path to the annotation of the wells
    method_of_fitting="NA",
    nrep=100,
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
    correct_negative="thr_correction", # if "thr_correction" it put a thr on the minimum value of the data with blank subracted, if "blank_correction" uses blank distrib to impute negative values
    thr_negative=0.01,  # used only if correct_negative == "thr_correction"
    multiple_scattering_correction=false, # if true uses the given calibration curve to fix the data
    method_multiple_scattering_correction="interpolation",
    calibration_OD_curve="NA",  #  the path to calibration curve to fix the data
    thr_lowess=0.05,
    penality_CI=8.0,
    size_bootstrap=0.7,
    blank_value = 0.0,
    blank_array = [0.0],
    optimizer=NLopt.LN_BOBYQA(),
    multistart=false,
    n_restart=50,
    auto_diff_method=nothing,
    cons=nothing,
    opt_params...
)


    if write_res == true
        mkpath(path_to_results)
    end

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
                model, # ode model to use
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
                model, # ode model to use
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
            model_to_test, # ode model to use
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

        parameter_of_optimization = hcat(parameter_of_optimization, temp_results_1[2])

        if errors_estimation == true && method_of_fitting != "Bootstrap"


            best_param = temp_results_1[2][3:(end-3)]
            best_param = convert.(Float64, best_param)
            temp_errors_of_optimization = NL_error_blanks(data, # dataset first row times second row OD
                string(well_name), # name of the well
                label_exp, #label of the experiment
                model, # ode model to use
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


    end


    if write_res == true

        CSV.write(
            string(path_to_results, label_exp, "_parameters_", model_string, ".csv"),
            Tables.table(Matrix(parameter_of_optimization)),
        )


    end

    Kimchi_res_one_file = (temp_results_1[1],parameter_of_optimization)


    return Kimchi_res_one_file, errors_of_optimization




end





"""
    function fit_NL_model_selection_file(
    label_exp::String, 
    path_to_data::String, 
    list_model_function::Any,
    list_lb_param::Vector{Float64}, 
    list_ub_param::Vector{Float64}; 
    path_to_annotation::Any = missing,
    method_of_fitting="MCMC",
    nrep=100,
    list_u0=lb_param .+ (ub_param .- lb_param) ./ 2,
    optimizer=BBO_adaptive_de_rand_1_bin_radiuslimited(), 
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
    correct_negative="thr_correction", 
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
    )



This function performs NL model selection of an array of NL models, it uses AIC or AICc depending on user inputs. This is done for a full .csv file.

# Arguments


- `label_exp::String`,  label of the experiment.
- `path_to_data::String`. path to the csv data frame. See documentation for formatting it.
-  `list_model_function::Any`: Array containing functions or strings of the NL models
-  `list_lb_param::Any`:Array of Lower bounds for the parameters (compatible with the models).
-  `list_ub_param::Any`:Array of Upper bounds for the parameters (compatible with the models).

# Key Arguments:
- `method_of_fitting="MCMC"`: String, how perform the NL fit. Options "MCMC","Bootstrap","Normal", and "Morris_sensitivity"
- `nrep=100`. Number of MCMC steps.
- `param= lb_param .+ (ub_param.-lb_param)./2`:Vector{Float64}, Initial guess for the model parameters.
- `optimizer =   BBO_adaptive_de_rand_1_bin_radiuslimited()` optimizer from optimizationBBO.
- `type_of_smoothing="rolling_avg"`: String, How to smooth the data, options: "NO" , "rolling avg" rolling average of the data, and "lowess".
- `pt_avg=7`: Number of points to generate the initial condition or do the rolling avg smoothing.
- `smoothing=false`: Whether to apply smoothing to the data or not.
- `type_of_loss:="RE" `: Type of loss function to be used. (options= "RE", "L2", "L2_derivative" and "blank_weighted_L2").
- `pt_smoothing_derivative=7`:Int,  Number of points for evaluation of specific growth rate. If <2 it uses interpolation algorithm otherwise a sliding window approach.
- `calibration_OD_curve="NA"`: String, The path where the .csv calibration data are located, used only if `multiple_scattering_correction=true`.
- `multiple_scattering_correction=false`: Bool, if true uses the given calibration curve to correct the data for muliple scattering.
- `method_multiple_scattering_correction="interpolation"`: String, How perform the inference of multiple scattering curve, options: "interpolation" or   "exp_fit" it uses an exponential fit from "Direct optical density determination of bacterial cultures in microplates for high-throughput screening applications"
-  `thr_lowess=0.05`: Float64 keyword argument of lowess smoothing
- `penality_CI=8.0`, used only in segementation to force the optimization to respect continuty on bonduar
-  `correction_AIC=true`: Bool, do finite samples correction of AIC.
-  `beta_smoothing_ms=2.0` penality  parameters for AIC (or AICc) evaluation.
-  `size_bootstrap=0.7`: Float, the fraction of data used each Bootstrap run. Used only if method is "Bootstrap"
- `write_res=false`: Bool, write the results in path_to_results folder.
- `blank_value = 0.0`: used only if `path_to_annotation = missing`and `blank_subtraction != "NO "`. It is used as average value of the blank.
- `blank_array = [0.0]`:used only if `path_to_annotation = missing`and `blank_subtraction != "NO "`. It is used as array of the blanks values.
-  `correct_negative="thr_correction"`: String, How to treat negative values after blank subtraction. If `"thr_correction"` it put a thr on the minimum value of the data with blank subracted, if `"blank_correction"` uses blank distribution to impute negative values, if `"remove"` the values are just removed..
- `do_blank_subtraction="NO"`: String, how perform the blank subtration, options "NO","avg_subtraction" (subtration of average value of blanks) and "time_avg" (subtration of  time average value of blanks).  
- `path_to_results= "NA"`:String, path to the folder where save the results.

# Output (if `results_NL_fit =NL_model_selection(...)`:


- a matrix with the following contents for each row : `[ "label of exp", "well", "param_1","param_2",..,"param_n","maximum specific gr using model","maximum specific gr using data", "objective function value (i.e. loss of the solution)"]` where ' "param_1","param_2",..,"param_n" ' .

"""
function fit_NL_model_selection_file(
    label_exp::String, #label of the experiment
    path_to_data::String, # path to the folder to analyze
    list_model_function::Any, # ode model to use
    list_u0;# initial guess param
    lb_param_array::Any = nothing, # lower bound param
    ub_param_array::Any = nothing, # upper bound param
    path_to_annotation::Any = missing,# path to the annotation of the wells
    method_of_fitting="NA",
    nrep=100,
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
    correct_negative="thr_correction", # if "thr_correction" it put a thr on the minimum value of the data with blank subracted, if "blank_correction" uses blank distrib to impute negative values
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
    optimizer=NLopt.LN_BOBYQA(),
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
            list_model_function, # ode model to use
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
            opt_params)



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
    Kimchi_res_one_file = ("NL_model_selection",parameter_of_optimization)

    return parameter_of_optimization




end



"""
    fit_NL_segmentation_file(
    label_exp::String, 
    path_to_data::String, 
    list_model_function::Any,
    list_lb_param::Vector{Vector{Float64}}, 
    list_ub_param::Vector{Vector{Float64}}, 
    n_change_points::Int;
    path_to_annotation::Any = missing,
    method_of_fitting="MCMC",
    nrep=100,
    list_u0=lb_param .+ (ub_param .- lb_param) ./ 2,
    optimizer=BBO_adaptive_de_rand_1_bin_radiuslimited(), 
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
    correct_negative="thr_correction", 
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
    )



This function performs NL model selection  on a segmented time series, it uses AIC or AICc depending on user inputs. This fuction works on an entire csv file.

# Arguments

- `label_exp::String`,  label of the experiment./
- `path_to_data::String`. path to the csv data frame. See documentation for formatting it.
-  `list_model_function::Any`: Array containing functions or strings of the NL models
-  `list_lb_param::Any`:Array of Lower bounds for the parameters (compatible with the models).
-  `list_ub_param::Any`:Array of Upper bounds for the parameters (compatible with the models).
-  `n_max_change_points::Int`: Number of change point used, the results will have different number of cp depending on the values of key argument 'type_of_detection' and 'fixed_cpd'

# Key Arguments:
- `method_of_fitting="MCMC"`: String, how perform the NL fit. Options "MCMC","Bootstrap","Normal", and "Morris_sensitivity"
- `nrep=100`. Number of MCMC steps.
- `param= lb_param .+ (ub_param.-lb_param)./2`:Vector{Float64}, Initial guess for the model parameters.
- `optimizer =   BBO_adaptive_de_rand_1_bin_radiuslimited()` optimizer from optimizationBBO.
- `type_of_smoothing="rolling_avg"`: String, How to smooth the data, options: "NO" , "rolling avg" rolling average of the data, and "lowess".
- `pt_avg=7`: Number of points to generate the initial condition or do the rolling avg smoothing.
- `smoothing=false`: Whether to apply smoothing to the data or not.
- `type_of_loss:="RE" `: Type of loss function to be used. (options= "RE", "L2", "L2_derivative" and "blank_weighted_L2").
- `pt_smoothing_derivative=7`:Int,  Number of points for evaluation of specific growth rate. If <2 it uses interpolation algorithm otherwise a sliding window approach.
- `calibration_OD_curve="NA"`: String, The path where the .csv calibration data are located, used only if `multiple_scattering_correction=true`.
- `multiple_scattering_correction=false`: Bool, if true uses the given calibration curve to correct the data for muliple scattering.
- `method_multiple_scattering_correction="interpolation"`: String, How perform the inference of multiple scattering curve, options: "interpolation" or   "exp_fit" it uses an exponential fit from "Direct optical density determination of bacterial cultures in microplates for high-throughput screening applications"
-  `thr_lowess=0.05`: Float64 keyword argument of lowess smoothing

- `penality_CI=2.0`, used only in segementation to force the optimization to respect continuty on bonduar
-  `correction_AIC=true`: Bool, do finite samples correction of AIC.
-  `beta_smoothing_ms=2.0` penality  parameters for AIC (or AICc) evaluation.
-  `size_bootstrap=0.7`: Float, the fraction of data used each Bootstrap run. Used only if method is "Bootstrap"
- `write_res=false`: Bool, write the results in path_to_results folder.
- `path_to_results= "NA"`:String, path to the folder where save the results.
- 'type_of_detection="slinding_win"': String, algorithm of cpd to use. Options '"slinding_win"' use a slinding window approach, '"lsdd"' uses least square density difference (LSDD) from ChangePointDetection.jl 
- 'type_of_curve="original"': String, on which curve is performed the change point detection algorithm. If '"original"' it use the original time series. With '"deriv"' it use the specific growth rate time series to perform the cdp.
- `method_peaks_detection="peaks_prominence"`: How the peak detection is performed on the dissimilarity curve.  `"peaks_prominence"` orders the peaks by prominence. `thr_scan` uses a threshold to choose the peaks
- `n_bins=40`: Int, used if `method_peaks_detection="thr_scan"` number of bins used to generate the threshold that has n_change_points peaks
- 'detect_number_cpd=true': Bool, if equal to true all the possible combination of lenght 1,2,...,n_change_points are tested and the best for AICc is returned.
- 'fixed_cpd=false': Bool If  true it returns the fitting using top n_change_points.
-  `correct_negative="thr_correction"`: # if "thr_correction" it put a thr on the minimum value of the data with blank subracted, if "blank_correction" uses blank distrib to impute negative values.
-  'win_size=14': Int, size of the windows used by the cdp algorithms

# Output (if `results_NL_fit =fit_NL_segmentation_file(...)`:


- an matrix with the following contents for each row : `[ "name of model", "well", "param_1","param_2",..,"param_n","maximum specific gr using model","maximum specific gr using data", "objective function value (i.e. loss of the solution)" "segment number"]` where ' "param_1","param_2",..,"param_n" ' .

"""
function fit_NL_segmentation_file(
    label_exp::String, #label of the experiment
    path_to_data::String, # path to the folder to analyze
    list_model_function::Any, # ode model to use
    list_u0,# initial guess param
    n_change_points::Int;
    list_lb_param::Vector{Vector{Float64}}=nothing, # lower bound param
    list_ub_param::Vector{Vector{Float64}}=nothing, # upper bound param
    path_to_annotation::Any = missing,# path to the annotation of the wells
    method_of_fitting="NA",
    nrep=100,
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
    correct_negative="thr_correction", # if "thr_correction" it put a thr on the minimum value of the data with blank subracted, if "blank_correction" uses blank distrib to impute negative values
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
    optimizer=NLopt.LN_BOBYQA(),
    auto_diff_method=nothing,
    multistart=false,
    n_restart=50,
    cons=nothing,
    opt_params...
)


    if write_res == true
        mkpath(path_to_results)
    end

    parameter_of_optimization = initialize_res_ms(list_ub_param, number_of_segment=n_change_points)
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

    end


    if write_res == true

        CSV.write(
            string(path_to_results, label_exp, "_parameters_ segmentation_NL.csv"),
            Tables.table(Matrix(parameter_of_optimization)),
        )


    end

    Kimchi_res_one_file = ("NL_segmentation",parameter_of_optimization)

    return Kimchi_res_one_file




end

export fit_NL_model_file
export fit_NL_model_selection_file
export fit_NL_segmentation_file
