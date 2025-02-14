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
    blank_array = [0.0],
    start_exp_win_thr=0.05
    )

Fits a logarithmic-linear model to data from a .csv file. The function assumes that the first column of the file represents time. It evaluates the specific growth rate, identifies an exponential window based on a statistical threshold, and performs a log-linear fitting.

# Arguments:

- `label_exp::String`: The label of the experiment.
- `path_to_data::String`: Path to the .csv file containing the data.
- `path_to_annotation::Any = missing`: Optional path to a .csv file with annotation data.

# Key Arguments:

- `path_to_results="NA"`: Path to the folder where results will be saved.
- `write_res=false`: Boolean flag to indicate whether to write the results to the specified path.
- `type_of_smoothing="rolling_avg"`: Method of data smoothing. Options are `"NO"`, `"rolling_avg"`, or `"lowess"`.
- `pt_avg=7`:  Number of points used in the rolling average smoothing.
- `pt_smoothing_derivative=7`: Number of points for evaluating specific growth rate. If less than 2, uses interpolation; otherwise, a sliding window approach is used.
- `pt_min_size_of_win=7`: Minimum size of the exponential windows in terms of the number of smoothed points.
- `type_of_win="maximum"`: Method for selecting the exponential phase window. Options are "maximum"` or `"global_thr"` "max_with_min_OD".
- `threshold_of_exp=0.9`: Threshold in quantile to define the exponential windows, between 0 and 1.
- `do_blank_subtraction="avg_blank"`: Method for blank subtraction. Options include `"NO"`, `"avg_subtraction"`, and `"time_avg"`.
- `blank_value=0.0`: Average value of the blank, used only if `do_blank_subtraction` is not `"NO"`.
- `blank_array=[0.0]`: Array of blank values, used only if `do_blank_subtraction` is not `"NO"`.
- `correct_negative="remove"`: Method for handling negative values after blank subtraction. Options are `"thr_correction"`, `"blank_correction"`, or `"remove"`.
- `thr_negative=0.01`: Threshold value for correcting negative values if `correct_negative` is `"thr_correction"`.
- `multiple_scattering_correction=false`: Flag indicating whether to correct for multiple scattering.
- `calibration_OD_curve="NA"`: Path to calibration data for multiple scattering correction, used if `multiple_scattering_correction` is true.
- `method_multiple_scattering_correction="interpolation"`: Method for correcting multiple scattering, options include `"interpolation"` or `"exp_fit"`.
- `thr_lowess=0.05`: Threshold for lowess smoothing.
- `verbose=false`: Flag to enable verbose output.
- `start_exp_win_thr=0.05`: Minimum OD value that should be reached to start the exponential window.

# Output:

- A data structure containing:
  1. `method`: Method used for fitting.
  2. A matrix with each row containing:
     - `label_exp`: Experiment label.
     - `name_well`: Name of the well or sample.
     - `start of exp win`: Start of the exponential window.
     - `end of exp win`: End of the exponential window.
     - `Maximum specific GR`: Maximum specific growth rate.
     - `specific GR`: Specific growth rate.
     - `2 sigma CI of GR`: Confidence interval of the growth rate (±2 sigma).
     - `doubling time`: Doubling time.
     - `doubling time - 2 sigma`: Doubling time minus 2 sigma.
     - `doubling time + 2 sigma`: Doubling time plus 2 sigma.
     - `intercept log-lin fitting`: Intercept of the log-linear fitting.
     - `2 sigma intercept`: Confidence interval of the intercept (±2 sigma).
     - `R^2`: Coefficient of determination (R-squared).

"""
function fit_one_file_Log_Lin(
    label_exp::String, 
    path_to_data::String; 
    path_to_annotation::Any=missing,
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
    start_exp_win_thr=0.05,
    blank_value=0.0,
    blank_array=[0.0],)


    fits= ()
    data_to_save = ()
    confidence_bands = ()






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
    Kinbiont_res_Log_Lin_files = ("Log-Lin", results_Log_Lin,fits,data_to_save,confidence_bands)


    return Kinbiont_res_Log_Lin_files




end




"""
    fit_file_ODE(
    label_exp::String,
    path_to_data::String,
    model::String, 
    param;
    path_to_annotation::Any=missing,
    Integration_method=Tsit5(), 
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

This function fits a ODE model to a csv file. The function assumes that the first column is the time, see the documentation for example of the data format. 
# Arguments:
- `label_exp::String`: The label of the experiment.
- `path_to_data::String`: String of the path where data file is located.
- `model::String`:String of the ODE to be fitted. See the documentation for the complete list.
- `param`:Vector{Float64}, Initial guess for the model parameters.


# Key Arguments:
- `Integration_method =Tsit5()' sciML Integration_method. If using piecewise model please use  'KenCarp4(autodiff=true)'.
- `optimizer = BBO_adaptive_de_rand_1_bin_radiuslimited()` optimizer from optimization.jl.
- `type_of_loss:="RE" `: Type of loss function to be used. Some options are "RE", "L2", "L2_derivative" and "blank_weighted_L2"), see documentation for the full list.
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
-  `correct_negative="remove"`: # if "thr_correction" it put a thr on the minimum value of the data with blank subracted, if "blank_correction" uses blank distrib to impute negative values.
- `blank_value = 0.0`: used only if `path_to_annotation = missing`and `do_blank_subtraction != "NO "`. It is used as average value of the blank.
- `blank_array = [0.0]`:used only if `path_to_annotation = missing`and `do_blank_subtraction != "NO "`. It is used as array of the blanks values.
-  `correct_negative="remove"`  ;: String, How to treat negative values after blank subtraction. If `"thr_correction"` it put a thr on the minimum value of the data with blank subracted, if `"blank_correction"` uses blank distribution to impute negative values, if `"remove"` the values are just removed..
- `do_blank_subtraction="NO"`: String, how perform the blank subtration, options "NO","avg_subtraction" (subtration of average value of blanks) and "time_avg" (subtration of  time average value of blanks).  
- `auto_diff_method=nothing`: method of differenzation, to be specified if required by the optimizer.
- `cons=nothing`. Equation constrains for optimization.
- `multistart=false`: use or not multistart optimization. Set to `true` uses Tik-Tak restart (from Benchmarking global optimizers, Arnoud et al 2019).
- `n_restart=50`: number of restart. Used if `multistart = true`.
- `opt_params...` :optional parameters of the required optimizer (e.g., `lb = [0.1, 0.3], ub =[9.0,1.0], maxiters=2000000`).



# Output:

- a data struct containing:
1. method string
2. matrix with the following contents for each row :`[ "name of model", "well", "param_1","param_2",..,"param_n","maximum specific gr using ode","maximum specific gr using data", "objective function value (i.e. loss of the solution)"]` where ' "param_1","param_2",..,"param_n" ' are the parameter of the selected ODE as in the documentation,
3. the fittings
4. the preprocessed data
"""
function fit_file_ODE(
    label_exp::String, #label of the experiment
    path_to_data::String, # path to the folder to analyze
    model::String, # string of the used model
    param;
    path_to_annotation::Any=missing,# path to the annotation of the wells
    Integration_method=Tsit5(), # selection of sciml Integration_method
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
            Integration_method=Integration_method, # selection of sciml Integration_method
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

        if multiple_scattering_correction == true

            data = correction_OD_multiple_scattering(data, calibration_OD_curve; method=method_multiple_scattering_correction)
        end

        data_to_save = (data_to_save...,data)
    


    end


    if write_res == true

        CSV.write(
            string(path_to_results, label_exp, "_parameters_", model, ".csv"),
            Tables.table(Matrix(parameter_of_optimization)),
        )


    end
    Kinbiont_res_one_file = ("ODE", parameter_of_optimization,fits,data_to_save)


    return Kinbiont_res_one_file




end

"""
    fit_file_custom_ODE(
    label_exp::String,
    path_to_data::String, 
    model::Any, 
    param::Vector{Float64},
    n_equation::Int;
    path_to_annotation::Any=missing,
    Integration_method=Tsit5(),
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
    blank_value=0.0,
    blank_array=[0.0],
    multistart=false,
    n_restart=50,
    optimizer=BBO_adaptive_de_rand_1_bin_radiuslimited(),
    auto_diff_method=nothing,
    cons=nothing,
    opt_params...
    )

Fits a customizable Ordinary Differential Equation (ODE) model to a dataset from a .csv file. 
# Arguments:

- `label_exp::String`: The label for the experiment.
- `path_to_data::String`: Path to the .csv file containing the data.
- `model::Any`: Function representing the ODE model to be fitted. See documentation for examples.
- `param::Vector{Float64}`: Initial guesses for the model parameters.
- `n_equation::Int`: Number of ODEs in the system.

# Key Arguments:

- `Integration_method=Tsit5()`: SciML Integration_method to use. For piecewise models, consider using `KenCarp4(autodiff=true)`.
- `optimizer=BBO_adaptive_de_rand_1_bin_radiuslimited()`: Optimizer for parameter fitting from optimization.jl.
- `loss_type="RE"`: Type of loss function. Options include `"RE"`, `"L2"`, `"L2_derivative"`, and `"blank_weighted_L2"`, see documentation for the full list.
- `path_to_annotation::Any=missing`: Path to a .csv file with annotation data. Required if `avg_replicate=true`.
- `path_to_results="NA"`: Path to the folder where results will be saved.
- `write_res=false`: Boolean flag indicating whether to write the results to the specified path.
- `smoothing=false`: Boolean flag for applying data smoothing.
- `type_of_smoothing="lowess"`: Method for smoothing the data. Options are `"NO"`, `"rolling_avg"`, or `"lowess"`.
- `pt_avg=1`: Number of points for rolling average smoothing.
- `pt_smooth_derivative=7`: Number of points for evaluating the specific growth rate. If less than 2, interpolation is used; otherwise, a sliding window approach is applied.
- `do_blank_subtraction="avg_blank"`: Method for blank subtraction. Options include `"NO"`, `"avg_subtraction"`, and `"time_avg"`.
- `blank_value=0.0`: Average value of the blank, used only if `do_blank_subtraction` is not `"NO"`.
- `blank_array=[0.0]`: Array of blank values, used only if `do_blank_subtraction` is not `"NO"`.
- `correct_negative="remove"`: Method for handling negative values after blank subtraction. Options are `"thr_correction"`, `"blank_correction"`, or `"remove"`.
- `thr_negative=0.01`: Threshold value for correcting negative values if `correct_negative` is `"thr_correction"`.
- `multiple_scattering_correction=false`: Flag indicating whether to correct for multiple scattering.
- `calibration_OD_curve="NA"`: Path to the .csv file with calibration data, used only if `multiple_scattering_correction=true`.
- `method_multiple_scattering_correction="interpolation"`: Method for correcting multiple scattering. Options include `"interpolation"` or `"exp_fit"`.
- `thr_lowess=0.05`: Threshold for lowess smoothing.
- `verbose=false`: Boolean flag for verbose output.
- `multistart=false`: Flag to use multistart optimization. Set to `true` uses Tik-Tak restart (from Benchmarking global optimizers, Arnoud et al 2019).
- `n_restart=50`: Number of restarts used if `multistart=true`.
- `auto_diff_method=nothing`: Method for differentiation, if required by the optimizer.
- `cons=nothing`: Constraints for optimization.
- `opt_params...`: Additional optional parameters for the optimizer, such as `lb`, `ub`, and `maxiters`.

# Output:

- A data structure containing:
  1. `method`: Method used for fitting.
  2. A matrix where each row includes:
     - `name of model`: Name of the fitted ODE model.
     - `well`: Identifier for the data well or sample.
     - `param_1`, `param_2`, ..., `param_n`: Parameters of the selected ODE model.
     - `maximum specific gr using ode`: Maximum specific growth rate calculated using the ODE.
     - `maximum specific gr using data`: Maximum specific growth rate calculated from the data.
     - `objective function value`: Loss value indicating the fit quality.
  3. The fittings (model outputs).
  4. The preprocessed data.

"""
function fit_file_custom_ODE(
    label_exp::String, #label of the experiment
    path_to_data::String, # path to the folder to analyze
    model::Any, # string of the used model
    param::Vector{Float64},# array of the array of the lower bound of the parameters
    n_equation::Int;
    path_to_annotation::Any=missing,# path to the annotation of the wells
    Integration_method=Tsit5(), # selection of sciml Integration_method
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
            Integration_method=Integration_method, # selection of sciml Integration_method
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
        if multiple_scattering_correction == true

            data = correction_OD_multiple_scattering(data, calibration_OD_curve; method=method_multiple_scattering_correction)
        end
        data_to_save = (data_to_save...,data)
    end


    if write_res == true

        CSV.write(
            string(path_to_results, label_exp, "_parameters_custom_ODE_fitting.csv"),
            Tables.table(Matrix(parameter_of_optimization)),
        )


    end

    Kinbiont_res_one_file = ("ODE", parameter_of_optimization,fits,data_to_save)


    return Kinbiont_res_one_file




end


"""
    ODE_model_selection_file(
    label_exp::String, 
    path_to_data::String, 
    models_list::Vector{String}, 
    param_array::Any;
    lb_param_array::Any=nothing, 
    ub_param_array::Any=nothing, 
    path_to_annotation::Any=missing,
    Integration_method=Tsit5(), 
    path_to_results="NA", 
    loss_type="L2", 
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

This function performs model selection for ordinary differential equations (ODEs) on a dataset provided in a .csv  file. 

# Arguments:

- `label_exp::String`: The label for the experiment.
- `path_to_data::String`: Path to the .csv file containing the data.
- `models_list::Vector{String}`: List of ODE models to evaluate. See documentation for the full list.
- `param_array::Any`: Initial guesses for the parameters of the models. It can be a matrix where each row corresponds to initial guesses for the models in `models_list`.

# Key Arguments:

- `lb_param_array::Any=nothing`: Lower bounds for the model parameters. Must be compatible with the models.
- `ub_param_array::Any=nothing`: Upper bounds for the model parameters. Must be compatible with the models.
- `path_to_annotation::Any=missing`: Path to the .csv  file with annotation data. Required if `avg_replicate=true`.
- `Integration_method=Tsit5()`: SciML Integration_method to use. For piecewise models, consider using `KenCarp4(autodiff=true)`.
- `optimizer=BBO_adaptive_de_rand_1_bin_radiuslimited()`: Optimizer used for parameter fitting from optimization.jl.
- `loss_type="L2"`: Type of loss function. Options include `"L2"`, `"RE"`, `"L2_derivative"`, and `"blank_weighted_L2"`, see documentation for the full list.
- `path_to_results="NA"`: Path to the folder where results will be saved.
- `write_res=false`: Boolean flag indicating whether to write the results to the specified path.
- `type_of_smoothing="lowess"`: Method for smoothing the data. Options are `"NO"`, `"rolling_avg"`, or `"lowess"`.
- `pt_avg=7`: Number of points for rolling average smoothing.
- `pt_smooth_derivative=7`: Number of points for evaluating the specific growth rate. If less than 2, interpolation is used; otherwise, a sliding window approach is applied.
- `do_blank_subtraction="avg_blank"`: Method for blank subtraction. Options include `"NO"`, `"avg_subtraction"`, and `"time_avg"`.
- `blank_value=0.0`: Average value of the blank, used only if `do_blank_subtraction` is not `"NO"`.
- `blank_array=[0.0]`: Array of blank values, used only if `do_blank_subtraction` is not `"NO"`.
- `correct_negative="remove"`: Method for handling negative values after blank subtraction. Options are `"thr_correction"`, `"blank_correction"`, or `"remove"`.
- `thr_negative=0.01`: Threshold value for correcting negative values if `correct_negative` is `"thr_correction"`.
- `multiple_scattering_correction=false`: Flag indicating whether to correct for multiple scattering.
- `calibration_OD_curve="NA"`: Path to the .csv file with calibration data, used only if `multiple_scattering_correction=true`.
- `method_multiple_scattering_correction="interpolation"`: Method for correcting multiple scattering. Options include `"interpolation"` or `"exp_fit"`.
- `thr_lowess=0.05`: Threshold for lowess smoothing.
- `correction_AIC=true`: Boolean flag indicating whether to apply finite sample correction to the Akaike Information Criterion (AIC).
- `beta_smoothing_ms=2.0`: Penalty parameter for AIC (or AICc) evaluation.
- `verbose=false`: Boolean flag for verbose output.
- `multistart=false`: Flag to use multistart optimization.  Set to `true` uses Tik-Tak restart (from Benchmarking global optimizers, Arnoud et al 2019).
- `n_restart=50`: Number of restarts used if `multistart=true`.
- `auto_diff_method=nothing`: Method for differentiation, if required by the optimizer.
- `cons=nothing`: Constraints for optimization.
- `opt_params...`: Additional optional parameters for the optimizer, such as `lb`, `ub`, and `maxiters`.

# Output:

- A data structure containing:
  1. `method`: The method used for model selection.
  2. A matrix where each row includes:
     - `name of model`: Name of the ODE model.
     - `well`: Identifier for the data well or sample.
     - `param_1`, `param_2`, ..., `param_n`: Parameters of the selected ODE model.
     - `maximum specific gr using ode`: Maximum specific growth rate calculated using the ODE.
     - `maximum specific gr using data`: Maximum specific growth rate calculated from the data.
     - `objective function value`: Loss value indicating the fit quality.
  3. The fittings (model outputs).
  4. The preprocessed data.

"""
function ODE_model_selection_file(
    label_exp::String, #label of the experiment
    path_to_data::String, # path to the folder to analyze
    models_list::Vector{String}, # ode model to use 
    param_array::Any;
    lb_param_array::Any=nothing, # lower bound param
    ub_param_array::Any=nothing, # upper bound param
    path_to_annotation::Any=missing,# path to the annotation of the wells
    Integration_method=Tsit5(), # selection of sciml Integration_method
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
            Integration_method=Integration_method, # selection of sciml Integration_method
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
        if multiple_scattering_correction == true

            data = correction_OD_multiple_scattering(data, calibration_OD_curve; method=method_multiple_scattering_correction)
        end
        data_to_save = (data_to_save...,data)
    end


    if write_res == true

        CSV.write(
            string(path_to_results, label_exp, "_parameters_model_selection.csv"),
            Tables.table(Matrix(parameter_of_optimization)),
        )


    end


    Kinbiont_res_one_file = ("ODE_model_selection", parameter_of_optimization,fits,data_to_save)

    return Kinbiont_res_one_file


end






"""
    segmentation_ODE_file(
    label_exp::String,
    path_to_data::String, 
    list_of_models::Vector{String}, 
    param_array::Any, 
    n_max_change_points::Int;
    lb_param_array::Any=nothing, 
    ub_param_array::Any=nothing, 
    path_to_annotation::Any=missing,
    detect_number_cpd=true,
    fixed_cpd=false,
    Integration_method=Tsit5(), 
    type_of_loss="L2", 
    type_of_detection="sliding_win",
    type_of_curve="original",
    do_blank_subtraction="avg_blank",
    correct_negative="remove",
    thr_negative=0.01,
    pt_avg=3, 
    smoothing=true, 
    path_to_results="NA",
    win_size=7, 
    pt_smooth_derivative=0,
    beta_smoothing_ms=2.0,
    avg_replicate=false,
    multiple_scattering_correction=false, 
    method_multiple_scattering_correction="interpolation",
    calibration_OD_curve="NA",  
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

This function performs model selection for ordinary differential equation (ODE) models while segmenting the time series data using change point detection algorithms.

# Arguments:

- `label_exp::String`: The label for the experiment.
- `path_to_data::String`: Path to the .csv file containing the data.
- `list_of_models::Vector{String}`: List of ODE models to evaluate. See documentation for the full list.
- `param_array::Any`: Initial guesses for the parameters of the models.
- `n_max_change_points::Int`: Maximum number of change points to use. The number of detected change points may vary based on the `type_of_detection` and `fixed_cpd` settings.

# Key Arguments:

- `lb_param_array::Any=nothing`: Lower bounds for the model parameters. Must be compatible with the models.
- `ub_param_array::Any=nothing`: Upper bounds for the model parameters. Must be compatible with the models.
- `path_to_annotation::Any=missing`: Path to the .csv file with annotation data. Required if `avg_replicate=true`.
- `detect_number_cpd=true`: Boolean flag indicating whether to detect the optimal number of change points. If `true`, all combinations from length 1 to `n_max_change_points` are tested, and the best is selected based on AICc.
- `fixed_cpd=false`: Boolean flag indicating whether to use a fixed number of change points. If `true`, the fitting will use exactly `n_max_change_points`.
- `Integration_method=Tsit5()`: SciML Integration_method to use. For piecewise models, consider using `KenCarp4(autodiff=true)`.
- `optimizer=BBO_adaptive_de_rand_1_bin_radiuslimited()`: Optimizer used for parameter fitting from optimization.jl.
- `type_of_loss="L2"`: Type of loss function. Options include `"L2"`, `"RE"`, `"L2_derivative"`, and `"blank_weighted_L2"`, see documentation for the full list.
- `path_to_results="NA"`: Path to the folder where results will be saved.
- `write_res=false`: Boolean flag indicating whether to write the results to the specified path.
- `type_of_smoothing="lowess"`: Method for smoothing the data. Options are `"NO"`, `"rolling_avg"`, or `"lowess"`.
- `pt_avg=3`: Number of points for rolling average smoothing.
- `pt_smooth_derivative=0`: Number of points for evaluating the specific growth rate. If less than 2, interpolation is used; otherwise, a sliding window approach is applied.
- `win_size=7`: Size of the window used by the change point detection algorithms.
- `type_of_detection="sliding_win"`: Change point detection algorithm to use. Options include `"sliding_win"` and `"lsdd"` (Least Squares Density Difference from ChangePointDetection.jl).
- `type_of_curve="original"`: Curve on which change point detection is performed. Options are `"original"` (original time series) or `"deriv"` (specific growth rate time series).
- `method_peaks_detection="peaks_prominence"`: Method for detecting peaks in the dissimilarity curve. Options are `"peaks_prominence"` (orders peaks by prominence) and `"thr_scan"` (uses a threshold to choose peaks).
- `n_bins=40`: Number of bins used to generate the threshold for peak detection if `method_peaks_detection="thr_scan"`.
- `smoothing=true`: Boolean flag to enable or disable smoothing.
- `multiple_scattering_correction=false`: Boolean flag indicating whether to correct for multiple scattering.
- `calibration_OD_curve="NA"`: Path to the .csv file with calibration data, used only if `multiple_scattering_correction=true`.
- `method_multiple_scattering_correction="interpolation"`: Method for correcting multiple scattering. Options include `"interpolation"` or `"exp_fit"`.
- `thr_lowess=0.05`: Threshold for lowess smoothing.
- `correct_negative="remove"`: Method for handling negative values after blank subtraction. Options include `"thr_correction"`, `"blank_correction"`, and `"remove"`.
- `blank_value=0.0`: Average value of the blank, used only if `do_blank_subtraction` is not `"NO"`.
- `blank_array=[0.0]`: Array of blank values, used only if `do_blank_subtraction` is not `"NO"`.
- `do_blank_subtraction="avg_blank"`: Method for blank subtraction. Options include `"NO"`, `"avg_subtraction"`, and `"time_avg"`.
- `correction_AIC=true`: Boolean flag indicating whether to apply finite sample correction to AIC (or AICc).
- `beta_smoothing_ms=2.0`: Penalty parameter for AIC (or AICc) evaluation.
- `auto_diff_method=nothing`: Method for differentiation, if required by the optimizer.
- `cons=nothing`: Constraints for optimization.
- `multistart=false`: Boolean flag to use multistart optimization. Set to `true` uses Tik-Tak restart (from Benchmarking global optimizers, Arnoud et al 2019).
- `n_restart=50`: Number of restarts used if `multistart=true`.
- `opt_params...`: Additional optional parameters for the optimizer, such as `lb`, `ub`, and `maxiters`.

# Output:

- A data structure containing:
  1. `method`: The method used for model selection.
  2. A matrix where each row includes:
     - `name of model`: Name of the ODE model.
     - `well`: Identifier for the data well or sample.
     - `param_1`, `param_2`, ..., `param_n`: Parameters of the selected ODE model.
     - `maximum specific gr using ode`: Maximum specific growth rate calculated using the ODE.
     - `maximum specific gr using data`: Maximum specific growth rate calculated from the data.
     - `objective function value`: Loss value indicating the fit quality.
  3. The fittings (model outputs).
  4. The preprocessed data.
  5. The change points detected in the data for each well.
  6. The AIC (or AICc) score for each well.

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
    Integration_method=Tsit5(), # selection of sciml Integration_method
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
            Integration_method=Integration_method, # selection of sciml Integration_method
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
        if multiple_scattering_correction == true

            data = correction_OD_multiple_scattering(data, calibration_OD_curve; method=method_multiple_scattering_correction)
        end
        data_to_save = (data_to_save...,data)
        cps = (cps...,temp_results_1[5])
        end


    if write_res == true

        CSV.write(
            string(
                path_to_results,
                label_exp,
                "_parameters_model_cpd_ncp_",
                n_max_change_points,
                ".csv",
            ),
            Tables.table(Matrix(parameter_of_optimization)),
        )


    end

    Kinbiont_res_one_file_segmentation = ("ODE_segmentation", parameter_of_optimization, fits, data_to_save, cps, vector_AIC)


    return Kinbiont_res_one_file_segmentation

end








"""

    segment_gr_analysis_file(
    path_to_data::String,
    label_exp::String; # Label of the experiment
    path_to_annotation=missing,
    n_max_change_points=0,
    type_of_smoothing="rolling_avg",
    pt_avg=7,
    pt_smoothing_derivative=7,
    multiple_scattering_correction=false,
    method_multiple_scattering_correction="interpolation",
    calibration_OD_curve="NA",
    thr_lowess=0.05,
    type_of_detection="slinding_win",
    type_of_curve="original",
    do_blank_subtraction="avg_blank",
    win_size=14,
    n_bins=40,
    path_to_results="NA",
    write_res=false,
    avg_replicate=false,
    blank_value=0.0,
    blank_array=[0.0],
    correct_negative="remove",
    thr_negative=0.01,
    verbose=false,
    )

This function analyzes change points in growth curve data from a .csv file and evaluates various metrics for each segment. 

# Arguments:

- `path_to_data::String`: Path to the .csv file containing the data. The file should have time values in the first row and OD values in the second row.
- `label_exp::String`: Label or identifier for the experiment.

# Key Arguments:

- `path_to_annotation::Any = missing`: Optional path to a .csv file containing annotations. If provided, it will be used for averaging replicates.
- `n_max_change_points::Int`: Maximum number of change points to consider. The actual number of change points in the results may vary based on the detection method and other parameters.
- `type_of_smoothing="rolling_avg"`: Method for smoothing the data. Options include `"NO"` (no smoothing), `"rolling_avg"` (rolling average), and `"lowess"` (locally weighted scatterplot smoothing).
- `pt_avg=7`: Size of the rolling average window for smoothing, applicable if `type_of_smoothing` is `"rolling_avg"`.
- `pt_smoothing_derivative=7`: Number of points used for the evaluation of specific growth rate. If less than 2, interpolation is used; otherwise, a sliding window approach is applied.
- `do_blank_subtraction="avg_blank"`: Method for blank subtraction. Options include `"NO"` (no blank subtraction), `"avg_blank"` (subtract average blank value), and `"time_avg"` (subtract time-averaged blank value).
- `blank_value=0.0`: Used as the average blank value if `path_to_annotation = missing` and `do_blank_subtraction != "NO"`.
- `blank_array=[0.0]`: Array of blank values used if `path_to_annotation = missing` and `do_blank_subtraction != "NO"`.
- `correct_negative="remove"`: Method for handling negative values after blank subtraction. Options include `"thr_correction"` (threshold correction), `"blank_correction"` (impute negative values based on blank distribution), and `"remove"` (remove negative values).
- `thr_negative=0.01`: Threshold used for correcting negative values if `correct_negative` is `"thr_correction"`. Values below this threshold will be adjusted to this value.
- `multiple_scattering_correction=false`: Whether to correct data for multiple scattering using a calibration curve.
- `calibration_OD_curve="NA"`: Path to the .csv file containing calibration data for multiple scattering correction.
- `method_multiple_scattering_correction="interpolation"`: Method for performing multiple scattering correction. Options include `"interpolation"` and `"exp_fit"`.
- `thr_lowess=0.05`: Parameter for lowess smoothing if `type_of_smoothing` is `"lowess"`.
- `type_of_detection="slinding_win"`: Method for detecting change points. Options include `"slinding_win"` (sliding window approach) and `"lsdd"` (least squares density difference).
- `type_of_curve="original"`: Specifies the curve used for change point detection. Options include `"original"` (raw time series) and `"deriv"` (specific growth rate time series).
- `method_peaks_detection="peaks_prominence"`: Method for detecting peaks in the dissimilarity curve. Options include `"peaks_prominence"` (orders peaks by prominence) and `"thr_scan"` (uses a threshold to select peaks).
- `n_bins=40`: Number of bins used to generate thresholds if `method_peaks_detection` is `"thr_scan"`.
- `win_size=14`: Size of the window used by change point detection algorithms.
- `path_to_results="NA"`: Path to the folder where results will be saved.
- `write_res=false`: Whether to write results to the specified folder.
- `avg_replicate=false`: Whether to average replicates if an annotation path is provided.
- `verbose=false`: Whether to print additional information during processing.

# Output:

If `res = segment_gr_analysis_file(...)`:

- `res[1]`: A string describing the method used for segmentation and analysis.
- `res[2]`: A matrix with rows containing: `[label_exp, name_well, max_specific_gr, min_specific_gr, t_of_max, od_of_max, max_deriv, min_deriv, start_OD_of_segment, delta_OD, segment_number]`.
- `res[3]`: List of change points for each well or sample.
- `res[4]`: Preprocessed data including smoothed values and corrected data.

"""
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
    do_blank_subtraction="avg_blank",
    win_size=14, #  
    n_bins=40,
    path_to_results="NA",# path where save results
    write_res=false, # write results
    avg_replicate=false,
    blank_value=0.0,
    blank_array=[0.0],
    correct_negative="remove",
    thr_negative=0.01,
    verbose=false,)





    # TEMPORARY results df
    results = ["label_exp", "well_name", "max_gr", "min_gr", "t_of_max", "od_of_max", "max_deriv", "min_deriv", "end_time", "delta_N", "segment_number"]

    if write_res == true
        mkpath(path_to_results)
    end
    names_of_annotated_df, properties_of_annotation, list_of_blank, list_of_discarded = reading_annotation(path_to_annotation)
    cp_list = ()
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
        cp_list = (cp_list...,temp_results_1[3])
        data_to_save = (data_to_save...,temp_results_1[4])

    end

    if write_res == true

        CSV.write(
            string(path_to_results,label_exp, "_results.csv"),
            Tables.table(Matrix(results)),
        )

    end


    Kinbiont_res_one_file = ("segment_analysis", results, cp_list, data_to_save)


    return Kinbiont_res_one_file
end


export fit_one_file_Log_Lin
export segment_gr_analysis_file
export fit_file_ODE
export fit_file_custom_ODE
export ODE_model_selection_file
export segmentation_ODE_file
export segment_gr_analysis_file
