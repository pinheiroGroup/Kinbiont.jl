using Optimization

"""
    fit_one_file_Log_Lin(
    label_exp::String, 
    path_to_data::String; 
    path_to_annotation::Any = missing,
    path_to_results="NA",
    write_res=false,
    type_of_smoothing="rolling_avg",
    pt_avg=7,
    pt_smoothing_derivative=7, 
    pt_min_size_of_win=7, 
    type_of_win="maximum", 
    threshold_of_exp=0.9,
    do_blank_subtraction="avg_blank",
    avg_replicate=false, 
    correct_negative="remove", 
    thr_negative=0.01, 
    multiple_scattering_correction=false, 
    method_multiple_scattering_correction="interpolation",
    calibration_OD_curve="NA",
    thr_lowess=0.05, 
    verbose=false,
    blank_value = 0.0,
    blank_array = [0.0],)


This function fits a logarithmic-linear model to a csv file. The function assumes that the first column is the time, see the documentation for example of the data format. It evaluate the specific growht rate, the with a statistical threshold it individuates a exponetial window and perform a-log lin fitting

# Arguments:

- `label_exp::String`: The label of the experiment.
-  `path_to_data::String`: Path to csv file containing the data    - `label_exp::String`: The label of the experiment.

# Key Arguments:

-  `path_to_annotation::Any = missing`: The path to the .csv of annotation .
-   `write_res=false`: Bool, write the results in path_to_results folder.
-  ` path_to_results= "NA"`:String, path to the folder where save the results.
- `average_replicate=false` Bool, perform or not the average of replicates. Works only if an annotation path is provided
-  `type_of_smoothing="rolling_avg"`: String, How to smooth the data, options: `"NO"` , `"rolling avg"` rolling average of the data, and `"lowess"`.
- `pt_avg=7`:Int, The number of points to do rolling average smoothing.
- `pt_smoothing_derivative=7`:Int,  Number of points for evaluation of specific growth rate. If <2 it uses interpolation algorithm otherwise a sliding window approach.
- `pt_min_size_of_win=7`:Int, The minimum size of the exponential windows in the number of smoothed points.
- `type_of_win="maximum"`:String, How the exponential phase window is selected ("maximum" or "global_thr").
- `threshold_of_exp=0.9`:Float, The threshold of the growth rate in quantile to define the exponential windows, a value between 0 and 1.
- `multiple_scattering_correction=false`:Bool, Whether or not correct the data qith a calibration curve.
- `calibration_OD_curve="NA"`: String, The path where the .csv calibration data are located, used only if `multiple_scattering_correction=true`.
- `multiple_scattering_correction=false`: Bool, if true uses the given calibration curve to correct the data for muliple scattering.
- `method_multiple_scattering_correction="interpolation"`: String, How perform the inference of multiple scattering curve, options: "interpolation" or   "exp_fit" it uses an exponential fit from "Direct optical density determination of bacterial cultures in microplates for high-throughput screening applications"
-  `thr_lowess=0.05`: Float64 keyword argument of lowees smoothing.
-  `correct_negative="remove"`: # if "thr_correction" it put a thr on the minimum value of the data with blank subracted, if "blank_correction" uses blank distrib to impute negative values.
- `blank_value = 0.0`: used only if `path_to_annotation = missing`and `do_blank_subtraction != "NO "`. It is used as average value of the blank.
- `blank_array = [0.0]`:used only if `path_to_annotation = missing`and `do_blank_subtraction != "NO "`. It is used as array of the blanks values.
-  `correct_negative="remove"`  ;: String, How to treat negative values after blank subtraction. If `"thr_correction"` it put a thr on the minimum value of the data with blank subracted, if `"blank_correction"` uses blank distribution to impute negative values, if `"remove"` the values are just removed.
-  `thr_negative=0.01`: FLoat: used only if `correct_negative == "thr_correction"` the data under this threshold will be changed to this value.
- `do_blank_subtraction="NO"`: String, how perform the blank subtration, options "NO","avg_subtraction" (subtration of average value of blanks) and "time_avg" (subtration of  time average value of blanks).  
- `start_exp_win_thr=0.05` minimum value (of OD) to consider the start of exp window

# Output:

- a matrix with the following contents in each row : `[label_exp, name_well, start of exp win,  end of exp win,  start of exp win, Maximum specific GR ,specific GR,  2 sigma  CI of GR, doubling time,doubling time - 2 sigma ,doubling time + 2 sigma  , intercept log-lin fitting, 2 sigma intercept ,R^2]`

"""
function fit_one_file_Log_Lin(
    label_exp::String, #label of the experiment
    path_to_data::String; # path to the folder to analyze
    path_to_annotation::Any=missing,# path to the annotation of the wells
    path_to_results="NA",# path where save results
    write_res=false, # write results
    type_of_smoothing="rolling_avg", # option, NO, gaussian, rolling avg
    pt_avg=7, # number of points to do smoothing average
    pt_smoothing_derivative=7, # number of poits to smooth the derivative
    pt_min_size_of_win=7, # minimum size of the exp windows in number of smooted points
    type_of_win="maximum", # how the exp. phase win is selected, "maximum" of "global_thr"
    threshold_of_exp=0.9, # threshold of growth rate in quantile to define the exp windows
    do_blank_subtraction="avg_blank", # string on how to use blank (NO,avg_subtraction,time_avg)
    avg_replicate=false, # if true the average between replicates is fitted. If false all replicate are fitted indipendelitly
    correct_negative="remove", # if "thr_correction" it put a thr on the minimum value of the data with blank subracted, if "blank_correction" uses blank distrib to impute negative values
    thr_negative=0.01, # used only if correct_negative == "thr_correction"
    multiple_scattering_correction=false, # if true uses the given calibration curve to fix the data
    method_multiple_scattering_correction="interpolation",
    calibration_OD_curve="NA", #  the path to calibration curve to fix the data
    thr_lowess=0.05, # keyword argument of lowees smoothing
    verbose=false,
    start_exp_win_thr=0.05,
    blank_value=0.0,
    blank_array=[0.0],)


    fits= ()
    data_to_save = ()
    confidence_bands = ()






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

    if write_res == true
        mkpath(path_to_results)
    end
    names_of_annotated_df, properties_of_annotation, list_of_blank, list_of_discarded = reading_annotation(path_to_annotation)


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
        index_tot = eachindex(data_values)
        index_tot = setdiff(index_tot, index_missing)
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
            start_exp_win_thr=start_exp_win_thr
        )

        if verbose == true
            println("the results are:")
            println(temp_results_1[2])
        end
        fits= (fits...,temp_results_1[3])
        data_to_save = (data_to_save...,temp_results_1[4] )
        confidence_bands = (confidence_bands...,temp_results_1[5] )
        results_Log_Lin = hcat(results_Log_Lin, temp_results_1[2])

    end

    if write_res == true

        CSV.write(
            string(path_to_results, label_exp, "_results.csv"),
            Tables.table(Matrix(results_Log_Lin)),
        )

    end
    Kimchi_res_Log_Lin_files = ("Log-Lin", results_Log_Lin,fits,data_to_save,confidence_bands)


    return Kimchi_res_Log_Lin_files




end




"""
    fit_file_ODE(
    label_exp::String, 
    path_to_data::String, 
    model::String,
    lb_param::Vector{Float64},
    ub_param::Vector{Float64};
    path_to_annotation::Any = missing,
    optimizer=BBO_adaptive_de_rand_1_bin_radiuslimited(), 
    integrator=Tsit5(),
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
    PopulationSize=300,
    maxiters=2000000,
    abstol=0.00001,
    thr_lowess=0.05,
    blank_value = 0.0,
    blank_array = [0.0],)

This function fits a ODE model to a csv file. The function assumes that the first column is the time, see the documentation for example of the data format. It evaluate the specific growht rate, the with a statistical threshold it individuates a exponetial window and perform a-log lin fitting
# Arguments:
- `data::Matrix{Float64}`: The dataset with the growth curve, where the first row represents times, and the second row represents the variable to fit (e.g., OD), see documentation.
- `label_exp::String`: The label of the experiment.
-  `model::String`:String of the ODE to be fitted. See the documentation for the complete list.
- `lb_param::Vector{Float64}`: Lower bounds of the model parameters.
- `ub_param::Vector{Float64}`: Upper bounds of the model parameters.

# Key Arguments:
- `param= lb_param .+ (ub_param.-lb_param)./2`:Vector{Float64}, Initial guess for the model parameters.
- `integrator =Tsit5()' sciML integrator. If using piecewise model please use  'KenCarp4(autodiff=true)'.
- `optimizer = BBO_adaptive_de_rand_1_bin_radiuslimited()` optimizer from optimizationBBO.
- `type_of_loss:="RE" `: Type of loss function to be used. (options= "RE", "L2", "L2_derivative" and "blank_weighted_L2").
- `average_replicate=false` Bool, perform or not the average of replicates. Works only if an annotation path is provided
- `path_to_annotation::Any = missing`: The path to the .csv of annotation .
- `write_res=false`: Bool, write the results in path_to_results folder.
- `path_to_results= "NA"`:String, path to the folder where save the results.
- `type_of_smoothing="rolling_avg"`: String, How to smooth the data, options: `"NO"` , `"rolling avg"` rolling average of the data, and `"lowess"`.
- `pt_avg=7`:Int, The number of points to do rolling average smoothing.
- `pt_smoothing_derivative=7`:Int,  Number of points for evaluation of specific growth rate. If <2 it uses interpolation algorithm otherwise a sliding window approach.
- `pt_min_size_of_win=7`:Int, The minimum size of the exponential windows in the number of smoothed points.
- `type_of_win="maximum"`:String, How the exponential phase window is selected ("maximum" or "global_thr").
- `threshold_of_exp=0.9`:Float, The threshold of the growth rate in quantile to define the exponential windows, a value between 0 and 1.
- `multiple_scattering_correction=false`:Bool, Whether or not correct the data qith a calibration curve.
- `calibration_OD_curve="NA"`: String, The path where the .csv calibration data are located, used only if `multiple_scattering_correction=true`.
- `multiple_scattering_correction=false`: Bool, if true uses the given calibration curve to correct the data for muliple scattering.
- `method_multiple_scattering_correction="interpolation"`: String, How perform the inference of multiple scattering curve, options: "interpolation" or   "exp_fit" it uses an exponential fit from "Direct optical density determination of bacterial cultures in microplates for high-throughput screening applications"
- `thr_lowess=0.05`: Float64 keyword argument of lowees smoothing.
- `PopulationSize =100`: Size of the population of the optimization
- `maxiters=2000000`: stop criterion, the optimization is stopped when the number of iterations is bigger than `maxiters`
- `abstol = 0.00001`: stop criterion, the optimization is stopped when the loss is lesser than `abstol`
-  `correct_negative="remove"`: # if "thr_correction" it put a thr on the minimum value of the data with blank subracted, if "blank_correction" uses blank distrib to impute negative values.
- `blank_value = 0.0`: used only if `path_to_annotation = missing`and `do_blank_subtraction != "NO "`. It is used as average value of the blank.
- `blank_array = [0.0]`:used only if `path_to_annotation = missing`and `do_blank_subtraction != "NO "`. It is used as array of the blanks values.
-  `correct_negative="remove"`  ;: String, How to treat negative values after blank subtraction. If `"thr_correction"` it put a thr on the minimum value of the data with blank subracted, if `"blank_correction"` uses blank distribution to impute negative values, if `"remove"` the values are just removed..
- `do_blank_subtraction="NO"`: String, how perform the blank subtration, options "NO","avg_subtraction" (subtration of average value of blanks) and "time_avg" (subtration of  time average value of blanks).  


# Output:

- an matrix with the following contents for each row :`[] "name of model", "well", "param_1","param_2",..,"param_n","maximum specific gr using ode","maximum specific gr using data", "objective function value (i.e. loss of the solution)"]` where ' "param_1","param_2",..,"param_n" ' are the parameter of the selected ODE as in the documentation.
"""
function fit_file_ODE(
    label_exp::String, #label of the experiment
    path_to_data::String, # path to the folder to analyze
    model::String, # string of the used model
    param;
    path_to_annotation::Any=missing,# path to the annotation of the wells
    integrator=Tsit5(), # selection of sciml integrator
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
    blank_value=0.0,
    blank_array=[0.0],
    multistart=false,
    n_restart=50,
    optimizer=BBO_adaptive_de_rand_1_bin_radiuslimited(),
    auto_diff_method=nothing,
    cons=nothing,
    opt_params...
)


    if write_res == true
        mkpath(path_to_results)
    end

    parameter_of_optimization = initialize_df_results(model)

    fits = ()
    data_to_save = ()

    names_of_annotated_df, properties_of_annotation, list_of_blank, list_of_discarded = reading_annotation(path_to_annotation)


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
        index_tot = eachindex(data_values)
        index_tot = setdiff(index_tot, index_missing)
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
            param; # upper bound param
            integrator=integrator, # selection of sciml integrator
            pt_avg=pt_avg, # numebr of the point to generate intial condition
            pt_smooth_derivative=pt_smooth_derivative,
            smoothing=smoothing, # the smoothing is done or not?
            type_of_loss=loss_type, # type of used loss 
            blank_array=blank_array, # data of all blanks
            multiple_scattering_correction=multiple_scattering_correction, # if true uses the given calibration curve to fix the data
            method_multiple_scattering_correction=method_multiple_scattering_correction,
            calibration_OD_curve=calibration_OD_curve, #  the path to calibration curve to fix the data
            thr_lowess=thr_lowess,
            type_of_smoothing=type_of_smoothing,
            multistart=multistart,
            n_restart=n_restart,
            optimizer=optimizer,
            auto_diff_method=auto_diff_method,
            cons=cons,
            opt_params...
        )


        if verbose == true
            println("the results are:")
            println(temp_results_1[2])
        end

        parameter_of_optimization = hcat(parameter_of_optimization, temp_results_1[2])
        fits = (fits...,temp_results_1[3] )
        data_to_save = (data_to_save...,data)
    


    end


    if write_res == true

        CSV.write(
            string(path_to_results, label_exp, "_parameters_", model, ".csv"),
            Tables.table(Matrix(parameter_of_optimization)),
        )


    end
    Kimchi_res_one_file = ("ODE", parameter_of_optimization,fits,data_to_save)


    return Kimchi_res_one_file




end

"""

    fit_file_custom_ODE(
    label_exp::String, 
    path_to_data::String,
    model::Any, 
    lb_param::Vector{Float64},
    ub_param::Vector{Float64},
    n_equation::Int;
    path_to_annotation::Any = missing,
    optimizer=BBO_adaptive_de_rand_1_bin_radiuslimited(), 
    integrator=Tsit5(),
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
    PopulationSize=300,
    maxiters=2000000,
    abstol=0.00001,
    thr_lowess=0.05,
    blank_value = 0.0,
    blank_array = [0.0],)

This function is designed for fitting an ordinary differential equation (ODE) model to a dataset in a csv file. . It utilizes a customizable ODE model, see documentation on how declare the model

# Arguments:


- `label_exp::String`: The label of the experiment.
-  `path_to_data::String`: Path to csv file containing the data
-  `model::Any`: Function of the ODE to be fitted. See the documentation for examples.
- `lb_param::Vector{Float64}`: Lower bounds of the model parameters.
- `ub_param::Vector{Float64}`: Upper bounds of the model parameters.
- `n_equation::Int`:  number ode in the system

# Key Arguments:

- `param= lb_param .+ (ub_param.-lb_param)./2`:Vector{Float64}, Initial guess for the model parameters.
- `integrator =Tsit5()' sciML integrator. If using piecewise model please use  'KenCarp4(autodiff=true)'.
- `optimizer = BBO_adaptive_de_rand_1_bin_radiuslimited()` optimizer from optimizationBBO.
- `type_of_loss:="RE" `: Type of loss function to be used. (options= "RE", "L2", "L2_derivative" and "blank_weighted_L2").
- `average_replicate=false` Bool, perform or not the average of replicates. Works only if an annotation path is provided
- `path_to_annotation::Any = missing`: The path to the .csv of annotation .
- `write_res=false`: Bool, write the results in path_to_results folder.
- `path_to_results= "NA"`:String, path to the folder where save the results.
- `type_of_smoothing="rolling_avg"`: String, How to smooth the data, options: `"NO"` , `"rolling avg"` rolling average of the data, and `"lowess"`.
- `pt_avg=7`:Int, The number of points to do rolling average smoothing.
- `pt_smoothing_derivative=7`:Int,  Number of points for evaluation of specific growth rate. If <2 it uses interpolation algorithm otherwise a sliding window approach.
- `pt_min_size_of_win=7`:Int, The minimum size of the exponential windows in the number of smoothed points.
- `type_of_win="maximum"`:String, How the exponential phase window is selected ("maximum" or "global_thr").
- `threshold_of_exp=0.9`:Float, The threshold of the growth rate in quantile to define the exponential windows, a value between 0 and 1.
- `multiple_scattering_correction=false`:Bool, Whether or not correct the data qith a calibration curve.
- `calibration_OD_curve="NA"`: String, The path where the .csv calibration data are located, used only if `multiple_scattering_correction=true`.
- `multiple_scattering_correction=false`: Bool, if true uses the given calibration curve to correct the data for muliple scattering.
- `method_multiple_scattering_correction="interpolation"`: String, How perform the inference of multiple scattering curve, options: "interpolation" or   "exp_fit" it uses an exponential fit from "Direct optical density determination of bacterial cultures in microplates for high-throughput screening applications"
- `thr_lowess=0.05`: Float64 keyword argument of lowees smoothing.
- `PopulationSize =100`: Size of the population of the optimization
- `maxiters=2000000`: stop criterion, the optimization is stopped when the number of iterations is bigger than `maxiters`
- `abstol = 0.00001`: stop criterion, the optimization is stopped when the loss is lesser than `abstol`
-  `correct_negative="remove"`: # if "thr_correction" it put a thr on the minimum value of the data with blank subracted, if "blank_correction" uses blank distrib to impute negative values.
- `blank_value = 0.0`: used only if `path_to_annotation = missing`and `blank_subtraction != "NO "`. It is used as average value of the blank.
- `blank_array = [0.0]`:used only if `path_to_annotation = missing`and `blank_subtraction != "NO "`. It is used as array of the blanks values.
-  `correct_negative="remove"`  ;: String, How to treat negative values after blank subtraction. If `"thr_correction"` it put a thr on the minimum value of the data with blank subracted, if `"blank_correction"` uses blank distribution to impute negative values, if `"remove"` the values are just removed..
- `do_blank_subtraction="NO"`: String, how perform the blank subtration, options "NO","avg_subtraction" (subtration of average value of blanks) and "time_avg" (subtration of  time average value of blanks).  


# Output:

- a matrix with the following contents for each row : `[ "name of model", "well", "param_1","param_2",..,"param_n","maximum specific gr using ode","maximum specific gr using data", "objective function value (i.e. loss of the solution)"]` where ' "param_1","param_2",..,"param_n" ' .

"""
function fit_file_custom_ODE(
    label_exp::String, #label of the experiment
    path_to_data::String, # path to the folder to analyze
    model::Any, # string of the used model
    param::Vector{Float64},# array of the array of the lower bound of the parameters
    n_equation::Int;
    path_to_annotation::Any=missing,# path to the annotation of the wells
    integrator=Tsit5(), # selection of sciml integrator
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
    blank_value=0.0,
    blank_array=[0.0],
    multistart=false,
    n_restart=50,
    optimizer=BBO_adaptive_de_rand_1_bin_radiuslimited(),
    auto_diff_method=nothing,
    cons=nothing,
    opt_params...
)


    if write_res == true
        mkpath(path_to_results)
    end

    parameter_of_optimization = initialize_df_results_ode_custom(ub_param)
    fits = ()
    data_to_save = ()

    names_of_annotated_df, properties_of_annotation, list_of_blank, list_of_discarded = reading_annotation(path_to_annotation)


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
        index_tot = eachindex(data_values)
        index_tot = setdiff(index_tot, index_missing)
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
            param, # lower bound param
            n_equation; # number ode in the system
            integrator=integrator, # selection of sciml integrator
            pt_avg=pt_avg, # number of the point to generate intial condition
            pt_smooth_derivative=pt_smooth_derivative,
            smoothing=smoothing, # the smoothing is done or not?
            type_of_loss=loss_type, # type of used loss 
            blank_array=blank_array, # data of all blanks
            multiple_scattering_correction=multiple_scattering_correction, # if true uses the given calibration curve to fix the data
            method_multiple_scattering_correction=method_multiple_scattering_correction,
            calibration_OD_curve=calibration_OD_curve,  #  the path to calibration curve to fix the data
            type_of_smoothing=type_of_smoothing,
            thr_lowess=0.05,
            multistart=multistart,
            n_restart=n_restart,
            optimizer=optimizer,
            auto_diff_method=auto_diff_method,
            cons=cons,
            opt_params...
        )

        well_results = reduce(vcat, temp_results_1[2])

        if verbose == true
            println("the results are:")
            println(well_results)
        end

        parameter_of_optimization = hcat(parameter_of_optimization, well_results)
        fits = (fits...,temp_results_1[3] )
        data_to_save = (data_to_save...,data)
    end


    if write_res == true

        CSV.write(
            string(path_to_results, label_exp, "_parameters_custom_ODE_fitting.csv"),
            Tables.table(Matrix(parameter_of_optimization)),
        )


    end

    Kimchi_res_one_file = ("ODE", parameter_of_optimization,fits,data_to_save)


    return Kimchi_res_one_file




end


"""
    ODE_model_selection_file(
    label_exp::String, 
    path_to_data::String, 
    models_list::Vector{String}, 
    lb_param_array::Any, 
    ub_param_array::Any; 
    path_to_annotation::Any = missing,
    optimizer=BBO_adaptive_de_rand_1_bin_radiuslimited(), 
    integrator=Tsit5(), 
    path_to_results="NA",
    loss_type="L2", 
    smoothing=false,
    type_of_smoothing="lowess",
    beta_smoothing_ms=2.0, 
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
    PopulationSize=300,
    maxiters=2000000,
    abstol=0.00001,
    thr_lowess=0.05,
    correction_AIC=true,
    blank_value = 0.0,
    blank_array = [0.0],)

This function performs model selection  of ODE for a full csv file.

# Arguments:

- `label_exp::String`: The label of the experiment.
-  `path_to_data::String`: Path to csv file containing the data    
-  `models_list::Vector{String}`: list of ODE model used
- `models_list::Vector{String}`: A vector of ODE models to evaluate.
- `lb_param_array::Any`: Lower bounds for the parameters (compatible with the models).
- `ub_param_array::Any`: Upper bounds


# Key Arguments:
- `param= lb_param .+ (ub_param.-lb_param)./2`:Vector{Float64}, Initial guess for the model parameters.
- `integrator =Tsit5()' sciML integrator. If using piecewise model please use  'KenCarp4(autodiff=true)'.
- `optimizer = BBO_adaptive_de_rand_1_bin_radiuslimited()` optimizer from optimizationBBO.
- `type_of_loss:="RE" `: Type of loss function to be used. (options= "RE", "L2", "L2_derivative" and "blank_weighted_L2").
- `average_replicate=false` Bool, perform or not the average of replicates. Works only if an annotation path is provided
- `path_to_annotation::Any = missing`: The path to the .csv of annotation .
- `write_res=false`: Bool, write the results in path_to_results folder.
- `path_to_results= "NA"`:String, path to the folder where save the results.
- `type_of_smoothing="rolling_avg"`: String, How to smooth the data, options: `"NO"` , `"rolling avg"` rolling average of the data, and `"lowess"`.
- `pt_avg=7`:Int, The number of points to do rolling average smoothing.
- `pt_smoothing_derivative=7`:Int,  Number of points for evaluation of specific growth rate. If <2 it uses interpolation algorithm otherwise a sliding window approach.
- `pt_min_size_of_win=7`:Int, The minimum size of the exponential windows in the number of smoothed points.
- `type_of_win="maximum"`:String, How the exponential phase window is selected ("maximum" or "global_thr").
- `threshold_of_exp=0.9`:Float, The threshold of the growth rate in quantile to define the exponential windows, a value between 0 and 1.
- `multiple_scattering_correction=false`:Bool, Whether or not correct the data qith a calibration curve.
- `calibration_OD_curve="NA"`: String, The path where the .csv calibration data are located, used only if `multiple_scattering_correction=true`.
- `multiple_scattering_correction=false`: Bool, if true uses the given calibration curve to correct the data for muliple scattering.
- `method_multiple_scattering_correction="interpolation"`: String, How perform the inference of multiple scattering curve, options: "interpolation" or   "exp_fit" it uses an exponential fit from "Direct optical density determination of bacterial cultures in microplates for high-throughput screening applications"
- `thr_lowess=0.05`: Float64 keyword argument of lowees smoothing.
- `PopulationSize =100`: Size of the population of the optimization
- `maxiters=2000000`: stop criterion, the optimization is stopped when the number of iterations is bigger than `maxiters`
- `abstol = 0.00001`: stop criterion, the optimization is stopped when the loss is lesser than `abstol`
-  `correct_negative="remove"`: # if "thr_correction" it put a thr on the minimum value of the data with blank subracted, if "blank_correction" uses blank distrib to impute negative values.
- `blank_value = 0.0`: used only if `path_to_annotation = missing`and `blank_subtraction != "NO "`. It is used as average value of the blank.
- `blank_array = [0.0]`:used only if `path_to_annotation = missing`and `blank_subtraction != "NO "`. It is used as array of the blanks values.
-  `correct_negative="remove"`  ;: String, How to treat negative values after blank subtraction. If `"thr_correction"` it put a thr on the minimum value of the data with blank subracted, if `"blank_correction"` uses blank distribution to impute negative values, if `"remove"` the values are just removed..
- `do_blank_subtraction="NO"`: String, how perform the blank subtration, options "NO","avg_subtraction" (subtration of average value of blanks) and "time_avg" (subtration of  time average value of blanks).  
-  `correction_AIC=true`: Bool, do finite samples correction of AIC.
-  `beta_smoothing_ms=2.0` penality  parameters for AIC (or AICc) evaluation.

# Output:

- an matrix with the following contents for each row : `[ "name of model", "well", "param_1","param_2",..,"param_n","maximum specific gr using ode","maximum specific gr using data", "objective function value (i.e. loss of the solution)"]` where ' "param_1","param_2",..,"param_n" ' .

"""
function ODE_model_selection_file(
    label_exp::String, #label of the experiment
    path_to_data::String, # path to the folder to analyze
    models_list::Vector{String}, # ode model to use 
    param_array::Any;
    lb_param_array::Any=nothing, # lower bound param
    ub_param_array::Any=nothing, # upper bound param
    path_to_annotation::Any=missing,# path to the annotation of the wells
    integrator=Tsit5(), # selection of sciml integrator
    path_to_results="NA", # path where save results
    loss_type="L2", # string of the type of the used loss
    smoothing=false, # 1 do smoothing of data with rolling average
    type_of_smoothing="lowess",
    beta_smoothing_ms=2.0, # penality for AIC evaluation
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
    correction_AIC=true,
    blank_value=0.0,
    blank_array=[0.0],
    multistart=false,
    n_restart=50,
    optimizer=BBO_adaptive_de_rand_1_bin_radiuslimited(),
    auto_diff_method=nothing,
    cons=nothing,
    opt_params...
)


    if write_res == true
        mkpath(path_to_results)
    end

    parameter_of_optimization = initialize_res_ms(param_array)
    fits = ()
    data_to_save = ()



    names_of_annotated_df, properties_of_annotation, list_of_blank, list_of_discarded = reading_annotation(path_to_annotation)

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
        index_tot = eachindex(data_values)
        index_tot = setdiff(index_tot, index_missing)
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
            param_array;
            lb_param_array=lb_param_array, # lower bound param
            ub_param_array=ub_param_array, # upper bound param
            integrator=integrator, # selection of sciml integrator
            pt_avg=pt_avg, # number of the point to generate intial condition
            beta_smoothing_ms=beta_smoothing_ms, # penality for AIC evaluation
            smoothing=smoothing, # the smoothing is done or not?
            type_of_smoothing=type_of_smoothing,
            thr_lowess=thr_lowess,
            type_of_loss=loss_type, # type of used loss 
            blank_array=blank_array, # data of all blanks
            pt_smooth_derivative=pt_smooth_derivative,
            multiple_scattering_correction=multiple_scattering_correction, # if true uses the given calibration curve to fix the data
            method_multiple_scattering_correction=method_multiple_scattering_correction,
            calibration_OD_curve=calibration_OD_curve, #  the path to calibration curve to fix the data
            verbose=verbose,
            correction_AIC=correction_AIC,
            multistart=multistart,
            n_restart=n_restart,
            optimizer=optimizer,
            auto_diff_method=auto_diff_method,
            cons=cons,
            opt_params...
        )

        vectorized_temp_results = expand_res(
            temp_results_1[end],
            param_array,
            string(well_name),
            label_exp
        )
        if verbose == true
            println("the results are:")
            println(vectorized_temp_results)
        end

        parameter_of_optimization = hcat(parameter_of_optimization, vectorized_temp_results)
        fits = (fits...,temp_results_1[3] )
        data_to_save = (data_to_save...,data)
    end


    if write_res == true

        CSV.write(
            string(path_to_results, label_exp, "_parameters_model_selection.csv"),
            Tables.table(Matrix(parameter_of_optimization)),
        )


    end


    Kimchi_res_one_file = ("ODE_model_selection", parameter_of_optimization,fits,data_to_save)

    return Kimchi_res_one_file


end






"""
    segmentation_ODE_file(
    label_exp::String, 
    path_to_data::String, 
    list_of_models::Vector{String},  
    lb_param_array::Any, 
    ub_param_array::Any,
    n_max_change_points::Int;
    path_to_annotation::Any = missing,
    detect_number_cpd=true,
    fixed_cpd=false,
    optimizer=BBO_adaptive_de_rand_1_bin_radiuslimited(), 
    integrator=Tsit5(), 
    type_of_loss="L2", 
    type_of_detection="sliding_win",
    type_of_curve="original",
    do_blank_subtraction="avg_blank",
    correct_negative="remove",
    thr_negative=0.01,
    pt_avg=1,
    smoothing=true, 
    path_to_results="NA",
    win_size=7,
    pt_smooth_derivative=0,
    beta_smoothing_ms=2.0,
    avg_replicate=false,
    multiple_scattering_correction="false",
    method_multiple_scattering_correction="interpolation",
    calibration_OD_curve="NA",  
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

This function performs model selection for ordinary differential equation (ODE) models while segmenting the time series in various part using change points detection algorithm, for a full csv file.

# Arguments:


- `label_exp::String`: The label of the experiment.
-  `path_to_data::String`: Path to csv file containing the data
- `list_of_models::Vector{String}`: A vector of ODE models to evaluate.
- `list_lb_param::Any`: Lower bounds for the parameters (compatible with the models).
- `list_ub_param::Any`: Upper bounds for the parameters (compatible with the models).
-  `n_max_change_points::Int`: Number of change point used, the results will have different number of cp depending on the values of key argument 'type_of_detection' and 'fixed_cpd'



# Key Arguments:
- `param= lb_param .+ (ub_param.-lb_param)./2`:Vector{Float64}, Initial guess for the model parameters.
- `integrator =Tsit5()' sciML integrator. If using piecewise model please use  'KenCarp4(autodiff=true)'.
- `optimizer = BBO_adaptive_de_rand_1_bin_radiuslimited()` optimizer from optimizationBBO.
- `type_of_loss:="RE" `: Type of loss function to be used. (options= "RE", "L2", "L2_derivative" and "blank_weighted_L2").
- `average_replicate=false` Bool, perform or not the average of replicates. Works only if an annotation path is provided
- `path_to_annotation::Any = missing`: The path to the .csv of annotation .
- `write_res=false`: Bool, write the results in path_to_results folder.
- `path_to_results= "NA"`:String, path to the folder where save the results.
- `type_of_smoothing="rolling_avg"`: String, How to smooth the data, options: `"NO"` , `"rolling avg"` rolling average of the data, and `"lowess"`.
- `pt_avg=7`:Int, The number of points to do rolling average smoothing.
- `pt_smoothing_derivative=7`:Int,  Number of points for evaluation of specific growth rate. If <2 it uses interpolation algorithm otherwise a sliding window approach.
- `pt_min_size_of_win=7`:Int, The minimum size of the exponential windows in the number of smoothed points.
- `type_of_win="maximum"`:String, How the exponential phase window is selected ("maximum" or "global_thr").
- `threshold_of_exp=0.9`:Float, The threshold of the growth rate in quantile to define the exponential windows, a value between 0 and 1.
- `multiple_scattering_correction=false`:Bool, Whether or not correct the data qith a calibration curve.
- `calibration_OD_curve="NA"`: String, The path where the .csv calibration data are located, used only if `multiple_scattering_correction=true`.
- `multiple_scattering_correction=false`: Bool, if true uses the given calibration curve to correct the data for muliple scattering.
- `method_multiple_scattering_correction="interpolation"`: String, How perform the inference of multiple scattering curve, options: "interpolation" or   "exp_fit" it uses an exponential fit from "Direct optical density determination of bacterial cultures in microplates for high-throughput screening applications"
- `thr_lowess=0.05`: Float64 keyword argument of lowees smoothing.
- `PopulationSize =100`: Size of the population of the optimization
- `maxiters=2000000`: stop criterion, the optimization is stopped when the number of iterations is bigger than `maxiters`
- `abstol = 0.00001`: stop criterion, the optimization is stopped when the loss is lesser than `abstol`
-  `correct_negative="remove"`: # if "thr_correction" it put a thr on the minimum value of the data with blank subracted, if "blank_correction" uses blank distrib to impute negative values.
- `blank_value = 0.0`: used only if `path_to_annotation = missing`and `blank_subtraction != "NO "`. It is used as average value of the blank.
- `blank_array = [0.0]`:used only if `path_to_annotation = missing`and `blank_subtraction != "NO "`. It is used as array of the blanks values.
-  `correct_negative="remove"`  ;: String, How to treat negative values after blank subtraction. If `"thr_correction"` it put a thr on the minimum value of the data with blank subracted, if `"blank_correction"` uses blank distribution to impute negative values, if `"remove"` the values are just removed..
- `do_blank_subtraction="NO"`: String, how perform the blank subtration, options "NO","avg_subtraction" (subtration of average value of blanks) and "time_avg" (subtration of  time average value of blanks).  
-  `correction_AIC=true`: Bool, do finite samples correction of AIC.
-  `beta_smoothing_ms=2.0` penality  parameters for AIC (or AICc) evaluation.
- 'type_of_detection="slinding_win"': String, algorithm of cpd to use. Options '"slinding_win"' use a slinding window approach, '"lsdd"' uses least square density difference (LSDD) from ChangePointDetection.jl 
- 'type_of_curve="original"': String, on which curve is performed the change point detection algorithm. If '"original"' it use the original time series. With '"deriv"' it use the specific growth rate time series to perform the cdp.
- `method_peaks_detection="peaks_prominence"`: How the peak detection is performed on the dissimilarity curve.  `"peaks_prominence"` orders the peaks by prominence. `thr_scan` uses a threshold to choose the peaks
- `n_bins=40`: Int, used if `method_peaks_detection="thr_scan"` number of bins used to generate the threshold that has n_change_points peaks
- 'detect_number_cpd=true': Bool, if equal to true all the possible combination of lenght 1,2,...,n_change_points are tested and the best for AICc is returned.
- 'fixed_cpd=false': Bool If  true it returns the fitting using top n_change_points.
-  'win_size=14': Int, size of the windows used by the cdo algorithms
-  'path_to_results="NA"':String, where to save the results.
-  'save_all_model=false': Bool, if true all the tested model are saved.
-  `correction_AIC=true`: Bool, do finite samples correction of AIC.
-  `beta_smoothing_ms=2.0` penality  parameters for AIC (or AICc) evaluation.



# Output:

- an matrix with the following contents for each row : `[ "name of model", "well", "param_1","param_2",..,"param_n","maximum specific gr using ode","maximum specific gr using data", "objective function value (i.e. loss of the solution)" "segment number"]` where ' "param_1","param_2",..,"param_n" ' .
"""
function segmentation_ODE_file(
    label_exp::String, #label of the experiment
    path_to_data::String, # path to the folder to analyze
    list_of_models::Vector{String}, # ode model to use 
    param_array::Any, #  param
    n_max_change_points::Int;
    lb_param_array::Any=nothing, # lower bound param
    ub_param_array::Any=nothing, # upper bound param
    path_to_annotation::Any=missing,# path to the annotation of the wells
    detect_number_cpd=true,
    fixed_cpd=false,
    integrator=Tsit5(), # selection of sciml integrator
    type_of_loss="L2", # type of used loss 
    type_of_detection="sliding_win",
    type_of_curve="original",
    do_blank_subtraction="avg_blank",
    correct_negative="remove",
    thr_negative=0.01,
    pt_avg=3, # number of the point to generate intial condition
    smoothing=true, # the smoothing is done or not?
    path_to_results="NA",
    win_size=7, # numebr of the point to generate intial condition
    pt_smooth_derivative=0,
    beta_smoothing_ms=2.0,
    avg_replicate=false,
    multiple_scattering_correction="false", # if true uses the given calibration curve to fix the data
    method_multiple_scattering_correction="interpolation",
    calibration_OD_curve="NA",  #  the path to calibration curve to fix the data
    write_res=false,
    save_all_model=false,
    method_peaks_detection="peaks_prominence",
    n_bins=40,
    type_of_smoothing="lowess",
    thr_lowess=0.05,
    verbose=false,
    correction_AIC=true,
    blank_value=0.0,
    blank_array=[0.0],
    optimizer=BBO_adaptive_de_rand_1_bin_radiuslimited(),
    multistart=false,
    n_restart=50,
    auto_diff_method=nothing,
    cons=nothing,
    opt_params...
)


    vector_AIC = ["AIC"]

    if write_res == true
        mkpath(path_to_results)
    end

    parameter_of_optimization = initialize_res_ms(param_array, number_of_segment=n_max_change_points + 1)
    fits = ()
    data_to_save = ()
    cps = ()

    names_of_annotated_df, properties_of_annotation, list_of_blank, list_of_discarded = reading_annotation(path_to_annotation)


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
        index_tot = eachindex(data_values)
        index_tot = setdiff(index_tot, index_missing)
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
            param_array,
            n_max_change_points;
            lb_param_array=lb_param_array, # lower bound param
            ub_param_array=ub_param_array, # upper bound param
            detect_number_cpd=detect_number_cpd,
            fixed_cpd=fixed_cpd,
            integrator=integrator, # selection of sciml integrator
            type_of_loss=type_of_loss, # type of used loss
            type_of_detection=type_of_detection,
            type_of_curve=type_of_curve,
            pt_avg=pt_avg, # number of the point to generate intial condition
            smoothing=smoothing, # the smoothing is done or not?
            path_to_results=path_to_results,
            win_size=win_size, # numebr of the point to generate intial condition
            pt_smooth_derivative=pt_smooth_derivative,
            beta_smoothing_ms=beta_smoothing_ms,
            multiple_scattering_correction=multiple_scattering_correction, # if true uses the given calibration curve to fix the data
            method_multiple_scattering_correction=method_multiple_scattering_correction,
            calibration_OD_curve=calibration_OD_curve,  #  the path to calibration curve to fix the data
            save_all_model=save_all_model,
            method_peaks_detection=method_peaks_detection,
            n_bins=n_bins,
            type_of_smoothing=type_of_smoothing,
            thr_lowess=thr_lowess,
            correction_AIC=correction_AIC,
            optimizer=optimizer,
            multistart=multistart,
            n_restart=n_restart,
            auto_diff_method=auto_diff_method,
            cons=cons,
            opt_params...
        )
        vector_AIC = hcat(vector_AIC, temp_results_1[end])

        vectorized_temp_results = expand_res_seg(
            temp_results_1[2],
            param_array,
            string(well_name),
            label_exp;
            number_of_segment=length(temp_results_1[2])
        )
        if verbose == true
            println("the results are:")
            println(vectorized_temp_results)
        end

        parameter_of_optimization = hcat(parameter_of_optimization, vectorized_temp_results)
        seg_fit = hcat(temp_results_1[4],temp_results_1[3])

        fits = (fits...,seg_fit )
        data_to_save = (data_to_save...,data)
        cps = (cps...,temp_results_1[5])
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

    Kimchi_res_one_file_segmentation = ("ODE_segmentation", parameter_of_optimization, fits, data_to_save, cps, vector_AIC)


    return Kimchi_res_one_file_segmentation

end











function segment_gr_analysis_file(
    path_to_data, # dataset first row times second row OD
    label_exp::String; #label of the experiment
    path_to_annotation=missing,
    n_max_change_points=0,
    type_of_smoothing="rolling_avg", # option, NO, gaussian, rolling avg
    pt_avg=7, # number of the point for rolling avg not used in the other cases
    pt_smoothing_derivative=7, # number of poits to smooth the derivative
    multiple_scattering_correction=false, # if true uses the given calibration curve to fix the data
    method_multiple_scattering_correction="interpolation",
    calibration_OD_curve="NA", #  the path to calibration curve to fix the data
    thr_lowess=0.05, # keyword argument of lowees smoothing
    type_of_detection="slinding_win",
    type_of_curve="original",
    win_size=14, #  
    n_bins=40,
    path_to_results="NA",# path where save results
    write_res=false, # write results
    avg_replicate=false,
    blank_value=0.0,
    blank_array=[0.0],
    correct_negative="remove",
    thr_negative=0.01,
    verbose=true,)





    # TEMPORARY results df
    results = ["label_exp", "well_name", "max_gr", "min_gr", "t_of_max", "od_of_max", "max_deriv", "min_deriv", "end_time", "delta_N", "segment_number"]

    if write_res == true
        mkpath(path_to_results)
    end
    names_of_annotated_df, properties_of_annotation, list_of_blank, list_of_discarded = reading_annotation(path_to_annotation)


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
        index_tot = eachindex(data_values)
        index_tot = setdiff(index_tot, index_missing)
        data = Matrix(transpose(hcat(times_data[index_tot], data_values[index_tot])))

        # correcting negative values after blank subtraction
        data = negative_value_correction(data,
            blank_array;
            method=correct_negative,
            thr_negative=thr_negative,)


        data = Matrix(transpose(hcat(data[1, :], data[2, :])))

        # inference


        temp_results_1 = segment_gr_analysis(
            data, # dataset first row times second row OD
            string(well_name), # name of the well
            label_exp; #label of the experiment
            n_max_change_points=n_max_change_points,
            type_of_smoothing=type_of_smoothing, # option, NO, gaussian, rolling avg
            pt_avg=pt_avg, # number of the point for rolling avg not used in the other cases
            pt_smoothing_derivative=pt_smoothing_derivative, # number of poits to smooth the derivative
            multiple_scattering_correction=multiple_scattering_correction, # if true uses the given calibration curve to fix the data
            method_multiple_scattering_correction=method_multiple_scattering_correction,
            calibration_OD_curve=calibration_OD_curve, #  the path to calibration curve to fix the data
            thr_lowess=thr_lowess, # keyword argument of lowees smoothing
            type_of_detection=type_of_detection,
            type_of_curve=type_of_curve,
            win_size=win_size, #  
            n_bins=n_bins,)

        if verbose == true
            println("the results are:")
            println(temp_results_1[2])
        end

        results = hcat(results, temp_results_1[2])

    end

    if write_res == true

        CSV.write(
            string(path_to_results, label_exp, "_results.csv"),
            Tables.table(Matrix(results)),
        )

    end


    Kimchi_res_segment_analysis = ("segment_analysis", results)


    return Kimchi_res_segment_analysis
end


export fit_one_file_Log_Lin
export segment_gr_analysis_file
export fit_file_ODE
export fit_file_custom_ODE
export ODE_model_selection_file
export segmentation_ODE_file
export segment_gr_analysis_file
