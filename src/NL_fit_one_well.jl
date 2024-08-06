


"""
    fit_NL_model(
    data::Matrix{Float64},
    name_well::String,
    label_exp::String,
    model_function::Any,
    u0;
    pt_avg=1,
    pt_smooth_derivative=7,
    smoothing=false,
    type_of_smoothing="rolling_avg",
    type_of_loss="RE",
    multiple_scattering_correction=false,
    method_multiple_scattering_correction="interpolation",
    calibration_OD_curve="NA",
    thr_lowess=0.05,
    penality_CI=3.0,
    optimizer=BBO_adaptive_de_rand_1_bin_radiuslimited(),
    multistart=false,
    n_restart=50,
    auto_diff_method=nothing,
    cons=nothing,
    opt_params...
    )

This function fits a nonlinear function to the time series input data of a single well.

# Arguments:

- `data::Matrix{Float64}`: The dataset containing the growth curve. The first row should represent time values, and the second row should represent the variable to fit (e.g., optical density). Refer to the documentation for proper formatting.

- `name_well::String`: The name of the well being analyzed.

- `label_exp::String`: The label for the experiment to identify the results.

- `model_function::Any`: The nonlinear model function to be used for fitting. This can be a custom function or a predefined model specified as a string (see documentation for available models).

- `u0::Vector{Float64}`: Initial guess for the model parameters.

# Key Arguments:

- `optimizer::Any = BBO_adaptive_de_rand_1_bin_radiuslimited()`: Optimizer for parameter estimation, from the BBO optimization library.

- `type_of_smoothing::String = "rolling_avg"`: Method for smoothing the data. Options include `"NO"` (no smoothing), `"rolling_avg"` (rolling average), and `"lowess"` (locally weighted scatterplot smoothing).

- `pt_avg::Int = 7`: Size of the rolling average window for smoothing if `type_of_smoothing` is `"rolling_avg"`.

- `smoothing::Bool = false`: Flag to apply data smoothing. Set to `true` to smooth the data; `false` to skip.

- `type_of_loss::String = "RE"`: Type of loss function used for optimization. Options include `"RE"` (relative error), `"L2"` (L2 norm), `"L2_derivative"`, and `"blank_weighted_L2"`.

- `blank_array::Vector{Float64} = zeros(100)`: Array containing data of blanks for correction.

- `pt_smooth_derivative::Int = 7`: Number of points for evaluating the specific growth rate. Uses interpolation if less than 2; otherwise, a sliding window approach is applied.

- `multiple_scattering_correction::Bool = false`: Flag to perform multiple scattering correction. Set to `true` to apply correction, requiring a calibration curve.

- `calibration_OD_curve::String = "NA"`: Path to a CSV file containing calibration data for optical density, used if `multiple_scattering_correction` is `true`.

- `method_multiple_scattering_correction::String = "interpolation"`: Method for performing multiple scattering correction. Options are `"interpolation"` or `"exp_fit"` (based on Meyers et al., 2018).

- `thr_lowess::Float64 = 0.05`: Threshold parameter for lowess smoothing.

- `penality_CI::Float64 = 3.0`: (Deprecated) Penalty for enforcing continuity at segment boundaries. Consider removing or updating.

- `auto_diff_method::Any = nothing`: Differentiation method for the optimizer, if required.

- `cons::Any = nothing`: Constraints for optimization.

- `multistart::Bool = false`: Flag to enable multistart optimization. Set to `true` to use multiple starting points.

- `n_restart::Int = 50`: Number of restarts for multistart optimization, used if `multistart` is `true`.

- `opt_params...`: Additional optional parameters for the optimizer (e.g., `lb` (lower bounds), `ub` (upper bounds), `maxiters`).

# Output:

- A data structure containing:
  1. `method`: A string describing the method used for fitting.
  2. Parameters array: `["name of model", "well", "param_1", "param_2", ..., "param_n", "maximum specific GR using NL", "maximum specific GR using data", "objective function value (i.e., loss of the solution)"]`, where `"param_1", "param_2", ..., "param_n"` are the parameters of the fitted model.
  3. The numerical solution of the fitted nonlinear model.
  4. Time coordinates corresponding to the fitted nonlinear model.

"""
function fit_NL_model(data::Matrix{Float64}, # dataset first row times second row OD
    name_well::String, # name of the well
    label_exp::String, #label of the experiment
    model_function::Any, #  model to use
    u0;# initial guess param    
    pt_avg=1, # number of the point to generate intial condition/
    pt_smooth_derivative=7,
    smoothing=false, # the smoothing is done or not?
    type_of_smoothing="rolling_avg",
    type_of_loss="RE", # type of used loss
    multiple_scattering_correction=false, # if true uses the given calibration curve to fix the data
    method_multiple_scattering_correction="interpolation",
    calibration_OD_curve="NA",  #  the path to calibration curve to fix the data
    thr_lowess=0.05,
    penality_CI=3.0,
    optimizer=BBO_adaptive_de_rand_1_bin_radiuslimited(),
    multistart=false,
    n_restart=50,
    auto_diff_method=nothing,
    cons=nothing,
    opt_params...
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

        model_string = KinBiont.NL_models[model_function].name
        model_function = KinBiont.NL_models[model_string].func


    else

        model_string = "custom"


    end

    loss_function =
        select_loss_function_NL(type_of_loss, data, penality_CI, model_function)


    # Solve the optimization problem
    sol = KinBiontSolve_NL(loss_function,
        u0,
        data;
        opt = optimizer,
        auto_diff_method=auto_diff_method,
        multistart = multistart,
        n_restart=n_restart,
        cons=cons,
        opt_params...)

    # evaluate the fitted  model
    fitted_model = model_function(sol, data[1, :])

    sol_fin, index_not_zero = remove_negative_value(fitted_model)

    data_th = transpose(hcat(data[1, index_not_zero], sol_fin))


    # max empirical gr
    if length(data[2, :]) < pt_smooth_derivative + 2
        max_em_gr = missing
        max_th_gr =missing
    else
        max_em_gr = maximum(specific_gr_evaluation(data, pt_smooth_derivative))
        max_th_gr = maximum(specific_gr_evaluation(Matrix(data_th), pt_smooth_derivative))

    end
    loss_value = sol.objective


    res_param = [[label_exp,name_well,model_string], [sol[1:end]], [max_th_gr, max_em_gr, loss_value]]


    res_param = reduce(vcat, reduce(vcat, res_param))

    KinBiont_res_one_well = ("NL", res_param, fitted_model, data[1, :])


    return KinBiont_res_one_well
end





"""
    fit_NL_model_with_sensitivity(
    data::Matrix{Float64},
    name_well::String,
    label_exp::String,
    model_function::Any,
    lb_param::Vector{Float64},
    ub_param::Vector{Float64};
    nrep=9,
    pt_avg=1,
    pt_smooth_derivative=7,
    smoothing=false,
    type_of_smoothing="rolling_avg",
    type_of_loss="RE",
    multiple_scattering_correction=false,
    method_multiple_scattering_correction="interpolation",
    calibration_OD_curve="NA",
    thr_lowess=0.05,
    write_res=false,
    penality_CI=3.0,
    optimizer=BBO_adaptive_de_rand_1_bin_radiuslimited(),
    multistart=false,
    n_restart=50,
    auto_diff_method=nothing,
    cons=nothing,
    opt_params...
    )

This function performs Morris sensitivity analysis on the nonlinear fit optimization. It evaluates the sensitivity of the fit parameters to variations in the initial guess, which is useful for quality checks of nonlinear model fits.

# Arguments:

- `data::Matrix{Float64}`: The growth curve data matrix. Time values are in the first row and the observable values (e.g., optical density) are in the second row.

- `name_well::String`: The name of the well for which the model is being fitted.

- `label_exp::String`: Label for the experiment, used for identifying results.

- `model_function::Any`: The nonlinear model function to be used for fitting. This can be a custom function or a predefined model name.

- `lb_param::Vector{Float64}`: Lower bounds for the model parameters. Defines the hyperspace for sensitivity analysis.

- `ub_param::Vector{Float64}`: Upper bounds for the model parameters. Defines the hyperspace for sensitivity analysis.

# Key Arguments:

- `nrep::Int = 100`: Number of steps for the Morris sensitivity analysis. Determines the number of sampling points.

- `optimizer::Any = BBO_adaptive_de_rand_1_bin_radiuslimited()`: Optimizer for parameter estimation, from the BBO optimization library.

- `type_of_smoothing::String = "rolling_avg"`: Method for smoothing the data. Options include `"NO"` (no smoothing), `"rolling_avg"` (rolling average), and `"lowess"` (locally weighted scatterplot smoothing).

- `pt_avg::Int = 7`: Size of the rolling average window for smoothing if `type_of_smoothing` is `"rolling_avg"`.

- `smoothing::Bool = false`: Flag to apply data smoothing. Set to `true` to enable smoothing; `false` to skip.

- `type_of_loss::String = "RE"`: Type of loss function used for optimization. Options include `"RE"` (relative error), `"L2"` (L2 norm), `"L2_derivative"`, and `"blank_weighted_L2"`.

- `blank_array::Vector{Float64} = zeros(100)`: Array containing data of blanks for correction.

- `pt_smooth_derivative::Int = 7`: Number of points for evaluating the specific growth rate. Uses interpolation if less than 2; otherwise, a sliding window approach is applied.

- `multiple_scattering_correction::Bool = false`: Flag to perform multiple scattering correction. Set to `true` to apply correction, requiring a calibration curve.

- `calibration_OD_curve::String = "NA"`: Path to a CSV file containing calibration data for optical density, used if `multiple_scattering_correction` is `true`.

- `method_multiple_scattering_correction::String = "interpolation"`: Method for performing multiple scattering correction. Options are `"interpolation"` or `"exp_fit"` (based on Meyers et al., 2018).

- `thr_lowess::Float64 = 0.05`: Threshold parameter for lowess smoothing.

- `write_res::Bool = false`: Flag to write results to a file. Set to `true` to enable file writing.

- `penality_CI::Float64 = 3.0`: (Deprecated) Penalty for enforcing continuity at segment boundaries. Consider removing or updating.

- `auto_diff_method::Any = nothing`: Differentiation method for the optimizer, if required.

- `cons::Any = nothing`: Constraints for optimization.

- `multistart::Bool = false`: Flag to enable multistart optimization. Set to `true` to use multiple starting points.

- `n_restart::Int = 50`: Number of restarts for multistart optimization, used if `multistart` is `true`.

- `opt_params...`: Additional optional parameters for the optimizer (e.g., `lb` (lower bounds), `ub` (upper bounds), `maxiters`).

# Output:

- A data structure containing:
  1. `method`: A string describing the method used for fitting and sensitivity analysis.
  2. Parameters array: `["name of model", "well", "param_1", "param_2", ..., "param_n", "maximum specific GR using NL", "maximum specific GR using data", "objective function value (i.e., loss of the solution)"]`, where `"param_1", "param_2", ..., "param_n"` are the parameters of the fitted model.
  3. The numerical solution of the fitted nonlinear model.
  4. Time coordinates corresponding to the fitted nonlinear model.
  5. Final parameters from the sensitivity analysis.
  6. Sensitivity analysis results including mean parameter values from the analysis.
  7. Standard deviation of the parameters from the sensitivity analysis.


"""
function fit_NL_model_with_sensitivity(data::Matrix{Float64}, # dataset first row times second row OD
    name_well::String, # name of the well
    label_exp::String, #label of the experiment
    model_function::Any, #  model to use
    lb_param::Vector{Float64}, # lower bound param
    ub_param::Vector{Float64}; # upper bound param
    nrep=9,
    pt_avg=1, # number of the point to generate intial condition
    pt_smooth_derivative=7,
    smoothing=false, # the smoothing is done or not?
    type_of_smoothing="rolling_avg",
    type_of_loss="RE", # type of used loss
    multiple_scattering_correction=false, # if true uses the given calibration curve to fix the data
    method_multiple_scattering_correction="interpolation",
    calibration_OD_curve="NA",  #  the path to calibration curve to fix the data
    thr_lowess=0.05,
    write_res=false,
    penality_CI=3.0,
    optimizer=BBO_adaptive_de_rand_1_bin_radiuslimited(),
    multistart=false,
    n_restart=50,
    auto_diff_method=nothing,
    cons=nothing,
    opt_params...
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

    if typeof(model_function) == String
        model_string = NL_models[model_function].name
        model_function = NL_models[model_string].func

    else

        model_string = "custom"


    end
    # Define the optimization problem LOSS


    loss_function =
        select_loss_function_NL(type_of_loss, data, penality_CI, model_function)


    if length(data[2, :]) < pt_smooth_derivative + 2
        max_em_gr = missing

    else
        max_em_gr = maximum(specific_gr_evaluation(data, pt_smooth_derivative))

    end
    fin_param = initialize_df_results_ode_custom(lb_param)

    param_combination =
        generation_of_combination_of_IC_morris(lb_param, ub_param, nrep)

    for i = 1:size(param_combination)[2]
        u0 = param_combination[:, i]



        sol = KinBiontSolve_NL(loss_function,
            u0,
            data;
            opt = optimizer,
            auto_diff_method=auto_diff_method,
            multistart = multistart,
            n_restart=n_restart,
            cons=cons,
            opt_params...)
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


        res_param = [[label_exp,name_well, model_string], [sol[1:end]], [max_th_gr, max_em_gr, loss_value]]
        res_param = reduce(vcat, reduce(vcat, res_param))



        fin_param = hcat(fin_param, res_param)


    end


    index_best = findmin(fin_param[end, 2:end])[2]


    best_res_param = fin_param[:, index_best+1]

    best_fitted_model = model_function(best_res_param[3:(end-3)], data[1, :])

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


    KinBiont_res_sensitivity_NL = ("NL_sensitivity", best_res_param, best_fitted_model, data[1, :], Matrix(param_combination))

    return KinBiont_res_sensitivity_NL
end





"""
    function fit_NL_model_bootstrap(
    data::Matrix{Float64},
    name_well::String,
    label_exp::String,
    model_function::Any,
    u0;
    lb_param=nothing,
    ub_param=nothing,
    nrep=100,
    size_bootstrap=0.7,
    pt_avg=1,
    pt_smooth_derivative=7,
    smoothing=false,
    type_of_smoothing="rolling_avg",
    type_of_loss="RE",
    multiple_scattering_correction=false,
    method_multiple_scattering_correction="interpolation",
    calibration_OD_curve="NA",
    thr_lowess=0.05,
    write_res=false,
    penality_CI=3.0,
    path_to_results="NA",
    optimizer=BBO_adaptive_de_rand_1_bin_radiuslimited(),
    multistart=false,
    n_restart=50,
    auto_diff_method=nothing,
    cons=nothing,
    opt_params...
    )

This function performs nonlinear (NL) fitting of the growth curve data using a bootstrap approach to evaluate confidence intervals and mitigate issues with poor initializations. 

# Arguments:

- `data::Matrix{Float64}`: The growth curve data matrix, where the first row contains time values, and the second row contains the observable values (e.g., optical density (OD)).

- `name_well::String`: The name of the well for which the model is being fitted.

- `label_exp::String`: The label for the experiment, used for identifying the results.

- `model_function::Any`: The nonlinear model function to be used for fitting. This can be a custom function or a predefined model name.

- `u0`: Initial guess for the model parameters.

# Key Arguments:

- `lb_param::Vector{Float64} = nothing`: Lower bounds for the model parameters. Defines the parameter space.

- `ub_param::Vector{Float64} = nothing`: Upper bounds for the model parameters. Defines the parameter space.

- `size_bootstrap::Float64 = 0.7`: Fraction of the data used for each bootstrap iteration.

- `nrep::Int = 100`: Number of bootstrap iterations to perform.

- `pt_avg::Int = 7`: Number of points to use for initial conditions or rolling average smoothing.

- `pt_smooth_derivative::Int = 7`: Number of points for evaluating the specific growth rate. Uses interpolation if less than 2; otherwise, applies a sliding window approach.

- `type_of_smoothing::String = "rolling_avg"`: Method for smoothing the data. Options include `"NO"` (no smoothing), `"rolling_avg"` (rolling average), and `"lowess"` (locally weighted scatterplot smoothing).

- `smoothing::Bool = false`: Flag to apply data smoothing. Set to `true` to enable smoothing; `false` to skip.

- `type_of_loss::String = "RE"`: Type of loss function used for optimization. Options include `"RE"` (relative error), `"L2"` (L2 norm), `"L2_derivative"`, and `"blank_weighted_L2"`.

- `blank_array::Vector{Float64} = zeros(100)`: Array containing data of blanks for correction.

- `calibration_OD_curve::String = "NA"`: Path to the CSV file containing calibration data for optical density, used if `multiple_scattering_correction` is `true`.

- `multiple_scattering_correction::Bool = false`: Flag to apply multiple scattering correction using the given calibration curve.

- `method_multiple_scattering_correction::String = "interpolation"`: Method for performing multiple scattering correction. Options are `"interpolation"` or `"exp_fit"` (based on Meyers et al., 2018).

- `thr_lowess::Float64 = 0.05`: Threshold parameter for lowess smoothing.

- `penality_CI::Float64 = 3.0`: (Deprecated) Penalty for enforcing continuity at segment boundaries. Consider removing or updating.

- `auto_diff_method::Any = nothing`: Differentiation method for the optimizer, if required.

- `cons::Any = nothing`: Constraints for optimization.

- `multistart::Bool = false`: Flag to enable multistart optimization. Set to `true` to use multiple starting points.

- `n_restart::Int = 50`: Number of restarts for multistart optimization, used if `multistart` is `true`.

- `optimizer::Any = BBO_adaptive_de_rand_1_bin_radiuslimited()`: Optimizer for parameter estimation, from the BBO optimization library.

- `opt_params...`: Additional optional parameters for the optimizer (e.g., `lb` (lower bounds), `ub` (upper bounds), `maxiters`).

- `write_res::Bool = false`: Flag to write results to a file. Set to `true` to enable file writing.

- `path_to_results::String = "NA"`: Path to the folder where results will be saved if `write_res` is `true`.

# Output:

- A data structure containing:
  1. `method`: A string describing the method used, including details of the bootstrap approach.
  2. Parameters array: `["name of model", "well", "param_1", "param_2", ..., "param_n", "maximum specific GR using NL", "maximum specific GR using data", "objective function value (i.e., loss of the solution)"]`, where `"param_1", "param_2", ..., "param_n"` are the parameters of the fitted model.
  3. The numerical solution of the fitted nonlinear model.
  4. Time coordinates corresponding to the fitted nonlinear model.
  5. Final parameters from the bootstrap analysis.
  6. Mean parameters from the bootstrap analysis.
  7. Standard deviation of parameters from the bootstrap analysis.
  8. (Optional) Results saved in the specified `path_to_results` folder if `write_res` is `true`.

"""
function fit_NL_model_bootstrap(data::Matrix{Float64}, # dataset first row times second row OD
    name_well::String, # name of the well
    label_exp::String, #label of the experiment
    model_function::Any, #  model to use
    u0;# initial guess param
    lb_param=nothing,
    ub_param=nothing,
    nrep=100,
    size_bootstrap=0.7,
    pt_avg=1, # number of the point to generate intial condition
    pt_smooth_derivative=7,
    smoothing=false, # the smoothing is done or not?
    type_of_smoothing="rolling_avg",
    type_of_loss="RE", # type of used loss
    multiple_scattering_correction=false, # if true uses the given calibration curve to fix the data
    method_multiple_scattering_correction="interpolation",
    calibration_OD_curve="NA",  #  the path to calibration curve to fix the data
    thr_lowess=0.05,
    write_res=false,
    penality_CI=3.0,
    path_to_results="NA",
    optimizer=BBO_adaptive_de_rand_1_bin_radiuslimited(),
    multistart=false,
    n_restart=50,
    auto_diff_method=nothing,
    cons=nothing,
    opt_params...
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

    if typeof(model_function) == String
        model_string = NL_models[model_function].name
        model_function = NL_models[model_string].func

    else

        model_string = "custom"


    end
    # Define the optimization problem LOSS



    if length(data[2, :]) > pt_smooth_derivative + 2
        max_em_gr = missing

    else
        max_em_gr = maximum(specific_gr_evaluation(data, pt_smooth_derivative))

    end

    fin_param = initialize_df_results_ode_custom(u0)


    for i = 1:nrep
        idxs = rand(1:1:size(data[2, :], 1), convert.(Int, floor(length(data[2, :]) * size_bootstrap)))
        times_boot = data[1, idxs]
        data_boot = data[2, idxs]

        idxs_2 = sortperm(times_boot)

        times_boot = times_boot[idxs_2]
        data_boot = data_boot[idxs_2]

        data_to_fit = Matrix(transpose(hcat(times_boot, data_boot)))


        loss_function =
            select_loss_function_NL(type_of_loss, data_to_fit, penality_CI, model_function)

        sol = KinBiontSolve_NL(loss_function,
            u0,
            data_to_fit;
           opt =  optimizer,
            auto_diff_method=auto_diff_method,
            cons=cons,
            multistart = multistart,
            n_restart=n_restart,
            opt_params...)
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


        res_param = [[label_exp,name_well, model_string], [sol[1:end]], [max_th_gr, max_em_gr, loss_value]]
        res_param = reduce(vcat, reduce(vcat, res_param))



        fin_param = hcat(fin_param, res_param)


    end


    index_best = findmin(fin_param[end, 2:end])[2]

    best_res_param = fin_param[:, index_best+1]
    best_fitted_model = model_function(best_res_param[3:(end-3)], data[1, :])

    if write_res == true
        mkpath(path_to_results)
        CSV.write(
            string(path_to_results, label_exp, "_", model_string, "_results_bootstrap.csv"),
            Tables.table(Matrix(fin_param)),
        )

    end

    quantile_loss = quantile(fin_param[end, 2:end], 0.90)
    index_loss = findall(fin_param[end, 2:end] .< quantile_loss)
    new_param_fin = fin_param[:, (index_loss.+1)]
    mean_param = [mean(new_param_fin[i, 2:end]) for i in 3:axes(new_param_fin)[1][end]]
    sd_param = [std(new_param_fin[i, 2:end]) for i in 3:axes(new_param_fin)[1][end]]


    KinBiont_res_bootstrap_NL = ("NL_bootstrap", best_res_param, best_fitted_model, data[1, :], fin_param, new_param_fin, mean_param, sd_param)


    return KinBiont_res_bootstrap_NL
end




"""
    function NL_error_blanks(
    data::Matrix{Float64},
    name_well::String,
    label_exp::String,
    model_function::Any,
    u0,
    blank_array::Vector{Float64};
    nrep=100,
    pt_avg=1,
    pt_smooth_derivative=7,
    smoothing=false,
    type_of_smoothing="rolling_avg",
    type_of_loss="RE",
    multiple_scattering_correction=false,
    method_multiple_scattering_correction="interpolation",
    calibration_OD_curve="NA",
    thr_lowess=0.05,
    write_res=false,
    penality_CI=3.0,
    optimizer=BBO_adaptive_de_rand_1_bin_radiuslimited(),
    multistart=false,
    n_restart=50,
    auto_diff_method=nothing,
    cons=nothing,
    opt_params...
    )

This function performs nonlinear (NL) fitting of the growth curve data while accounting for blank data noise. It executes `nrep` iterations to estimate the posterior distribution of the parameters.

# Arguments:

- `data::Matrix{Float64}`: The dataset with the growth curve. The first row contains time values, and the second row contains optical density (OD) measurements.

- `name_well::String`: The name of the well for which the model is being fitted.

- `label_exp::String`: The label for the experiment, used to identify the results.

- `model_function::Any`: The nonlinear model function to be used for fitting. This can be a custom function or one of the predefined models.

- `u0`: Initial guess for the model parameters.

- `blank_array::Vector{Float64}`: Array of blank data used to model the noise in the measurements.

- `nrep::Int = 100`: Number of iterations for estimating the posterior distribution of the parameters.

- `pt_avg::Int = 1`: Number of points used for generating initial conditions or performing smoothing.

- `pt_smooth_derivative::Int = 7`: Number of points for evaluating the specific growth rate. Uses interpolation if less than 2; otherwise, applies a sliding window approach.

- `smoothing::Bool = false`: Whether to apply data smoothing. Set to `true` to enable smoothing; `false` to skip.

- `type_of_smoothing::String = "rolling_avg"`: Method for smoothing the data. Options include `"NO"` (no smoothing), `"rolling_avg"` (rolling average), and `"lowess"` (locally weighted scatterplot smoothing).

- `type_of_loss::String = "RE"`: Type of loss function used for optimization. Options include `"RE"` (relative error), `"L2"` (L2 norm), `"L2_derivative"`, and `"blank_weighted_L2"`.

- `multiple_scattering_correction::Bool = false`: Whether to apply multiple scattering correction using the given calibration curve.

- `method_multiple_scattering_correction::String = "interpolation"`: Method for performing multiple scattering correction. Options are `"interpolation"` or `"exp_fit"` (based on Meyers et al., 2018).

- `calibration_OD_curve::String = "NA"`: Path to the CSV file containing calibration data for optical density, used if `multiple_scattering_correction` is `true`.

- `thr_lowess::Float64 = 0.05`: Threshold parameter for lowess smoothing.

- `penality_CI::Float64 = 3.0`: (Deprecated) Penalty for enforcing continuity at segment boundaries. Consider removing or updating.

- `auto_diff_method::Any = nothing`: Differentiation method for the optimizer, if required.

- `cons::Any = nothing`: Constraints for optimization.

- `multistart::Bool = false`: Flag to enable or disable multistart optimization. Set to `true` to use multiple starting points.

- `n_restart::Int = 50`: Number of restarts for multistart optimization, used if `multistart` is `true`.

- `opt_params...`: Additional optional parameters for the optimizer (e.g., `lb` (lower bounds), `ub` (upper bounds), `maxiters`).

- `write_res::Bool = false`: Flag to write results to a file. Set to `true` to enable file writing.

- `path_to_results::String = "NA"`: Path to the folder where results will be saved if `write_res` is `true`.

# Output:

- A data structure containing:
  1. `method`: A string describing the method used, including details of the blank noise approach.
  2. Parameters array: `["name of model", "well", "param_1", "param_2", ..., "param_n", "maximum specific GR using NL", "maximum specific GR using data", "objective function value (i.e., loss of the solution)"]`, where `"param_1", "param_2", ..., "param_n"` are the parameters of the fitted model.
  3. The numerical solution of the fitted nonlinear model.
  4. Time coordinates corresponding to the fitted nonlinear model.
  5. Estimated parameters from the iterations considering blank noise.
  6. Posterior distribution summary (mean and standard deviation of parameters) if applicable.
  7. (Optional) Results saved in the specified `path_to_results` folder if `write_res` is `true`.


"""
function NL_error_blanks(data::Matrix{Float64}, # dataset first row times second row OD
    name_well::String, # name of the well
    label_exp::String, #label of the experiment
    model_function::Any, #  model to use
    u0,# initial guess param
    blank_array::Vector{Float64}; # upper bound param
    nrep=100,
    pt_avg=1, # number of the point to generate intial condition
    pt_smooth_derivative=7,
    smoothing=false, # the smoothing is done or not?
    type_of_smoothing="rolling_avg",
    type_of_loss="RE", # type of used loss
    multiple_scattering_correction=false, # if true uses the given calibration curve to fix the data
    method_multiple_scattering_correction="interpolation",
    calibration_OD_curve="NA",  #  the path to calibration curve to fix the data
    thr_lowess=0.05,
    write_res=false,
    penality_CI=3.0,
    optimizer=BBO_adaptive_de_rand_1_bin_radiuslimited(),
    multistart=false,
    n_restart=50,
    auto_diff_method=nothing,
    cons=nothing,
    opt_params...
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





    if length(data[2, :]) > pt_smooth_derivative + 2
        max_em_gr = missing

    else
        max_em_gr = maximum(specific_gr_evaluation(data, pt_smooth_derivative))

    end
    fin_param = initialize_df_results_ode_custom(u0)


    for i = 1:nrep

        data[2, :] = data[2, :] .+ rand(Normal(0.0, std(blank_array)), length(data[2, :]))
        sol_fin, index_not_zero = remove_negative_value(data[2, :])

        data = transpose(hcat(data[1, index_not_zero], sol_fin))




        loss_function =
            select_loss_function_NL(type_of_loss, data, penality_CI, model_function)

        sol = KinBiontSolve_NL(loss_function,
            u0,
            data;
           opt = optimizer,
            auto_diff_method=auto_diff_method,
            multistart=multistart,
            n_restart=n_restart,
            cons=cons,
            opt_params...)
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


        res_param = [[label_res,name_well, model_string], [sol[1:end]], [max_th_gr, max_em_gr, loss_value]]
        res_param = reduce(vcat, reduce(vcat, res_param))



        fin_param = hcat(fin_param, res_param)


    end


    index_best = findmin(fin_param[end, 2:end])[2]

    best_res_param = fin_param[:, index_best+1]
    best_fitted_model = model_function(best_res_param[3:(end-3)], data[1, :])

    if write_res == true
        mkpath(path_to_results)
        CSV.write(
            string(path_to_results, label_exp, "_", model_string, "_results_bootstrap.csv"),
            Tables.table(Matrix(fin_param)),
        )

    end


    quantile_loss = quantile(fin_param[end, 2:end], 0.90)
    index_loss = findall(fin_param[end, 2:end] .< quantile_loss)
    new_param_fin = fin_param[:, (index_loss.+1)]
    mean_param = [mean(new_param_fin[i, 2:end]) for i in 3:axes(new_param_fin)[1][end]]
    sd_param = [std(new_param_fin[i, 2:end]) for i in 3:axes(new_param_fin)[1][end]]
    CI_param_low = [quantile(new_param_fin[i, 2:end], 0.1) for i in 3:axes(new_param_fin)[1][end]]
    CI_param_up = [quantile(new_param_fin[i, 2:end], 0.9) for i in 3:axes(new_param_fin)[1][end]]

    return best_res_param, best_fitted_model, fin_param, new_param_fin, mean_param, sd_param, CI_param_low, CI_param_up
end












"""
    NL_model_selection(
    data::Matrix{Float64},
    name_well::String,
    label_exp::String,
    list_model_function::Any,
    list_u0;
    lb_param_array::Any = nothing,
    ub_param_array::Any = nothing,
    method_of_fitting="NA",
    nrep=100,
    size_bootstrap=0.7,
    pt_avg=1,
    pt_smooth_derivative=7,
    smoothing=false,
    type_of_smoothing="rolling_avg",
    type_of_loss="RE",
    multiple_scattering_correction=false,
    method_multiple_scattering_correction="interpolation",
    calibration_OD_curve="NA",
    thr_lowess=0.05,
    write_res=false,
    beta_smoothing_ms=2.0,
    penality_CI=8.0,
    correction_AIC=false,
    optimizer=BBO_adaptive_de_rand_1_bin_radiuslimited(),
    auto_diff_method=nothing,
    multistart=false,
    n_restart=50,
    cons=nothing,
    opt_params...
    )

This function performs nonlinear (NL) model selection from an array of NL models using AIC or AICc, depending on user inputs. It performs `nrep` iterations to estimate the posterior distribution of parameters. It uses the blank distribution as noise.

# Arguments:

- `data::Matrix{Float64}`: The dataset with the growth curve. The first row represents times, and the second row represents the variable to fit (e.g., OD).

- `name_well::String`: Name of the well.

- `label_exp::String`: Label of the experiment.

- `list_model_function::Any`: Array containing functions or strings representing the NL models.

- `list_u0`: Array of initial guesses for the parameters for each model.

# Key Arguments:

- `lb_param_array::Any = nothing`: Array of lower bounds for the parameters, compatible with the models.

- `ub_param_array::Any = nothing`: Array of upper bounds for the parameters, compatible with the models.

- `method_of_fitting::String = "Normal"`: Method of performing the NL fit. Options are `"Bootstrap"`, `"Normal"`, and `"Morris_sensitivity"`.

- `nrep::Int = 100`: Number of iterations for estimating the posterior distribution.

- `size_bootstrap::Float = 0.7`: Fraction of data used in each bootstrap run, applicable if `method_of_fitting` is `"Bootstrap"`.

- `pt_avg::Int = 1`: Number of points to generate initial conditions or perform smoothing.

- `pt_smooth_derivative::Int = 7`: Number of points for evaluating the specific growth rate. Uses interpolation if less than 2; otherwise, applies a sliding window approach.

- `smoothing::Bool = false`: Whether to apply smoothing to the data. Set to `true` to enable smoothing; `false` to skip.

- `type_of_smoothing::String = "rolling_avg"`: Method for smoothing the data. Options are `"NO"`, `"rolling_avg"`, and `"lowess"`.

- `type_of_loss::String = "RE"`: Type of loss function used for optimization. Options include `"RE"`, `"L2"`, `"L2_derivative"`, and `"blank_weighted_L2"`.

- `blank_array::Vector{Float64} = zeros(100)`: Data of all blanks in a single array.

- `pt_smoothing_derivative::Int = 7`: Number of points for evaluating the specific growth rate. Uses interpolation if less than 2; otherwise, uses a sliding window approach.

- `calibration_OD_curve::String = "NA"`: Path to the CSV file containing calibration data for optical density, used if `multiple_scattering_correction` is `true`.

- `multiple_scattering_correction::Bool = false`: Whether to apply multiple scattering correction using the given calibration curve.

- `method_multiple_scattering_correction::String = "interpolation"`: Method for performing multiple scattering correction. Options are `"interpolation"` or `"exp_fit"`.

- `thr_lowess::Float64 = 0.05`: Threshold parameter for lowess smoothing.

- `beta_smoothing_ms::Float64 = 2.0`: Penalty parameter for AIC (or AICc) evaluation.

- `penality_CI::Float64 = 8.0`: Penalty for enforcing continuity at segment boundaries.

- `correction_AIC::Bool = false`: Whether to apply finite sample correction to AIC.

- `auto_diff_method::Any = nothing`: Differentiation method for the optimizer, if required.

- `multistart::Bool = false`: Flag to enable or disable multistart optimization. Set to `true` to use multiple starting points.

- `n_restart::Int = 50`: Number of restarts for multistart optimization, used if `multistart` is `true`.

- `opt_params...`: Additional optional parameters for the optimizer (e.g., `lb` (lower bounds), `ub` (upper bounds), `maxiters`).

- `write_res::Bool = false`: Flag to write results to a file. Set to `true` to enable file writing.

- `path_to_results::String = "NA"`: Path to the folder where results will be saved if `write_res` is `true`.

# Output:

- A data structure containing:
  1. `method`: A string describing the method used for model selection.
  2. Parameters array: `["name of model", "well", "param_1", "param_2", ..., "param_n", "maximum specific GR using NL", "maximum specific GR using data", "objective function value (i.e., loss of the solution)"]`, where `"param_1", "param_2", ..., "param_n"` are the parameters of the fitted model.
  3. The numerical solution of the fitted nonlinear model.
  4. Time coordinates corresponding to the fitted nonlinear model.
  5. The AIC values for all models.
  6. The loss value of the top model.





"""
function NL_model_selection(data::Matrix{Float64}, # dataset first row times second row OD
    name_well::String, # name of the well
    label_exp::String, #label of the experiment
    list_model_function::Any, #  model to use
    list_u0;
    lb_param_array::Any = nothing, # lower bound param
    ub_param_array::Any = nothing, # upper bound param
    method_of_fitting="NA",
    nrep=100,
    size_bootstrap=0.7,
    pt_avg=1, # number of the point to generate intial condition
    pt_smooth_derivative=7,
    smoothing=false, # the smoothing is done or not?
    type_of_smoothing="rolling_avg",
    type_of_loss="RE", # type of used loss
    multiple_scattering_correction=false, # if true uses the given calibration curve to fix the data
    method_multiple_scattering_correction="interpolation",
    calibration_OD_curve="NA",  #  the path to calibration curve to fix the data
    thr_lowess=0.05,
    write_res=false,
    beta_smoothing_ms=2.0,
    penality_CI=8.0,
    correction_AIC=false,
    optimizer=BBO_adaptive_de_rand_1_bin_radiuslimited(),
    auto_diff_method=nothing,
    multistart=false,
    n_restart=50,
    cons=nothing,
    opt_params...
)
    score_res = ["AIC"]
    top_score = 10^20
    top_model = Vector{Any}
    top_fitted_sol = Vector{Any}
    top_loss = 10^20
    opt_param_temp = copy(opt_params)
    temp_param_ub = nothing
    temp_param_lb = nothing


    for mm in 1:eachindex(list_model_function)[end]

        model_to_test = list_model_function[mm]

        if !isnothing(lb_param_array)
            temp_param_lb = lb_param_array[mm]
            temp_param_ub = ub_param_array[mm]


            opt_param_temp = (opt_params...,
                lb=temp_param_lb,
                ub=temp_param_ub,
                )
        end

        u0 = list_u0[mm]

        if method_of_fitting == "Bootstrap"



            temp_res = fit_NL_model_bootstrap(data, # dataset first row times second row OD
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
                opt_param_temp...
            )

            n_param = length(u0)
            temp_AIC = AICc_evaluation2(n_param, beta_smoothing_ms, data[2, :], temp_res[2][end], correction=correction_AIC)

            temp = [model_to_test, temp_res[end], temp_AIC]
            score_res = hcat(score_res, temp_AIC)
            if mm == 1

                top_score = copy(temp_AIC)
                top_model = copy(temp_res[2])
                top_fitted_sol = copy(temp_res[3])
                top_loss = copy(temp_res[2][end])

            elseif top_score > temp_AIC
                top_score = copy(temp_AIC)
                top_model = copy(temp_res[2])
                top_fitted_sol = copy(temp_res[3])
                top_loss = copy(temp_res[2][end])

            end


        elseif method_of_fitting == "Morris_sensitivity"


            temp_res = fit_NL_model_with_sensitivity(data, # dataset first row times second row OD
                name_well, # name of the well
                label_exp, #label of the experiment
                model_to_test, #  model to use
                temp_param_lb, # lower bound param
                temp_param_ub; # upper bound param
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
                opt_param_temp...
            )

            n_param = length(temp_param_lb)

            temp_AIC = AICc_evaluation2(n_param, beta_smoothing_ms, data[2, :], temp_res[2][end], correction=correction_AIC)

            score_res = hcat(score_res, temp_AIC)

            if mm == 1

                top_score = copy(temp_AIC)
                top_model = copy(temp_res[2])
                top_fitted_sol = copy(temp_res[3])
                top_loss = copy(temp_res[2][end])


            elseif top_score > temp_AIC
                top_score = copy(temp_AIC)
                top_model = copy(temp_res[2])
                top_fitted_sol = copy(temp_res[3])
                top_loss = copy(temp_res[2][end])

            end

    

        else



            temp_res = fit_NL_model(data, # dataset first row times second row OD
                name_well, # name of the well
                label_exp, #label of the experiment
                model_to_test, #  model to use
                u0;
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
                penality_CI=penality_CI,
                auto_diff_method=auto_diff_method,
                multistart=multistart,
                n_restart=n_restart,
                cons=cons,
                opt_param_temp...
            )


            n_param = length(u0)

            temp_AIC = AICc_evaluation2(n_param, beta_smoothing_ms, data[2, :], temp_res[2][end], correction=correction_AIC)

            score_res = hcat(score_res, temp_AIC)
            if mm == 1

                top_score = copy(temp_AIC)
                top_model = copy(temp_res[2])
                top_fitted_sol = copy(temp_res[3])
                top_loss = copy(temp_res[2][end])


            elseif top_score > temp_AIC

                top_score = copy(temp_AIC)
                top_model = copy(temp_res[2])
                top_fitted_sol = copy(temp_res[3])
                top_loss = copy(temp_res[2][end])

            end



        end


    end
    KinBiont_res_NL_model_selection = ("NL_model_selection", top_model, top_fitted_sol, data[1, :], score_res, top_loss)
    return KinBiont_res_NL_model_selection
end

"""
    selection_NL_fixed_interval(
    data_testing::Matrix{Float64},
    name_well::String,
    label_exp::String,
    list_of_models::Vector{String},
    list_u0,
    intervals_changepoints::Any;
    lb_param_array::Any = nothing,
    ub_param_array::Any = nothing,
    type_of_loss="L2",
    method_of_fitting="NA",
    smoothing=false,
    size_bootstrap=0.7,
    nrep=100,
    type_of_smoothing="lowess",
    thr_lowess=0.05,
    pt_avg=1,
    pt_smooth_derivative=0,
    multiple_scattering_correction=false,
    method_multiple_scattering_correction="interpolation",
    calibration_OD_curve="NA",
    beta_smoothing_ms=2.0,
    penality_CI=8.0,
    correction_AIC=true,
    optimizer=BBO_adaptive_de_rand_1_bin_radiuslimited(),
    multistart=false,
    n_restart=50,
    auto_diff_method=nothing,
    cons=nothing,
    opt_params...
    )

This function fits a segmented nonlinear (NL) model to a curve, using specified change points.

# Arguments:

- `data_testing::Matrix{Float64}`: The dataset with the growth curve, where the first row represents times, and the second row represents the variable to fit (e.g., OD).

- `name_well::String`: Name of the well.

- `label_exp::String`: Label of the experiment.

- `list_of_models::Vector{String}`: Array containing functions or strings representing the NL models.

- `list_u0`: Initial guesses for the parameters for each model.

- `intervals_changepoints::Any`: Array containing the change points for segmentation, e.g., `[0.0, 10.0, 30.0]`.

# Key Arguments:

- `lb_param_array::Any = nothing`: Array of lower bounds for the parameters, compatible with the models.

- `ub_param_array::Any = nothing`: Array of upper bounds for the parameters, compatible with the models.

- `type_of_loss::String = "L2"`: Type of loss function used. Options include `"L2"`, `"RE"`, `"L2_derivative"`, and `"blank_weighted_L2"`.

- `method_of_fitting::String = "Normal"`: Method for performing the NL fit. Options are `"Normal"`, `"Bootstrap"`, and `"Morris_sensitivity"`.

- `nrep::Int = 100`: Number of iterations for Bootstrap or Morris sensitivity analysis.

- `smoothing::Bool = false`: Whether to apply smoothing to the data.

- `type_of_smoothing::String = "lowess"`: Method for smoothing the data. Options include `"NO"`, `"rolling_avg"`, and `"lowess"`.

- `thr_lowess::Float64 = 0.05`: Threshold parameter for lowess smoothing.

- `pt_avg::Int = 1`: Number of points for generating initial conditions or performing smoothing.

- `pt_smooth_derivative::Int = 0`: Number of points for evaluating the specific growth rate. Uses interpolation if less than 2; otherwise, applies a sliding window approach.

- `multiple_scattering_correction::Bool = false`: Whether to apply multiple scattering correction using the given calibration curve.

- `method_multiple_scattering_correction::String = "interpolation"`: Method for performing multiple scattering correction. Options are `"interpolation"` or `"exp_fit"`.

- `calibration_OD_curve::String = "NA"`: Path to the CSV file with calibration data for optical density, used if `multiple_scattering_correction` is `true`.

- `beta_smoothing_ms::Float64 = 2.0`: Penalty parameter for AIC (or AICc) evaluation.

- `penality_CI::Float64 = 8.0`: Penalty for enforcing continuity at segment boundaries.

- `correction_AIC::Bool = true`: Whether to apply finite sample correction to AIC.

- `size_bootstrap::Float = 0.7`: Fraction of data used in each bootstrap run, applicable if `method_of_fitting` is `"Bootstrap"`.

- `auto_diff_method::Any = nothing`: Differentiation method for the optimizer, if required.

- `cons::Any = nothing`: Constraints for optimization.

- `multistart::Bool = false`: Whether to use multistart optimization. Set to `true` to use multiple starting points.

- `n_restart::Int = 50`: Number of restarts for multistart optimization, used if `multistart` is `true`.

- `opt_params...`: Additional optional parameters for the optimizer (e.g., `lb` (lower bounds), `ub` (upper bounds), `maxiters`).

# Output:

If `results_NL_fit = selection_NL_fixed_interval(...)`:

- `results_NL_fit[1]`: Parameters of each segment.

- `results_NL_fit[2]`: Numerical solutions of the fit.

- `results_NL_fit[3]`: Time coordinates of the numerical solutions of the fit.

- `results_NL_fit[4]`: Loss of the best model.

"""
function selection_NL_fixed_interval(
    data_testing::Matrix{Float64}, # dataset first row times second row OD
    name_well::String, # name of the well
    label_exp::String, #label of the experiment
    list_of_models::Vector{String}, #  models to use
    list_u0,
    intervals_changepoints::Any;
    lb_param_array::Any = nothing, # lower bound param
    ub_param_array::Any = nothing, # upper bound param  
    type_of_loss="L2", # type of used loss
    method_of_fitting="NA", # selection of sciml integrator
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
    penality_CI=8.0,
    correction_AIC=true,
    optimizer=BBO_adaptive_de_rand_1_bin_radiuslimited(),
    multistart=false,
    n_restart=50,
    auto_diff_method=nothing,
    cons=nothing,
    opt_params...
)
    interval_changepoints = copy(intervals_changepoints)
    interval_changepoints = push!(interval_changepoints, data_testing[1, 1])

    interval_changepoints = push!(interval_changepoints, data_testing[1, end])

    interval_changepoints = sort(interval_changepoints)
    bc = [data_testing[1, end], data_testing[2, end]]
    param_out = Vector{Vector{Any}}()
    composed_sol = Type{Any}
    composed_time = Type{Any}
    loss_to_use = ""
    sum_of_loss = 0.0
    for i = (length(interval_changepoints)):-1:2

        if i == 2 && i != (length(interval_changepoints))
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
            list_of_models, #  model to use
            list_u0;# initial guess param
            lb_param_array =lb_param_array, # lower bound param
            ub_param_array =ub_param_array,# upper bound param
            method_of_fitting=method_of_fitting,
            nrep=nrep,
            size_bootstrap=size_bootstrap,
            pt_avg=pt_avg, # number of the point to generate intial condition
            pt_smooth_derivative=pt_smooth_derivative,
            smoothing=smoothing, # the smoothing is done or not?
            type_of_smoothing=type_of_smoothing,
            type_of_loss=loss_to_use, # type of used loss
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
            opt_params...
        )

        # param of the best model
        temp_res_win = model_selection_results[2]
        sum_of_loss = sum_of_loss + model_selection_results[2][end]


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
    composed_time, composed_sol = remove_replicate_data(composed_time, composed_sol)

    return param_out, composed_sol, composed_time, sum_of_loss

end

"""
    segmentation_NL(
    data_testing::Matrix{Float64},
    name_well::String,
    label_exp::String,
    list_of_models::Any,
    list_u0,
    n_change_points::Int;
    lb_param_array::Any=nothing,
    ub_param_array::Any=nothing,
    type_of_loss="L2_fixed_CI",
    method_of_fitting="NA",
    type_of_detection="sliding_win",
    type_of_curve="original",
    smoothing=false,
    nrep=100,
    type_of_smoothing="lowess",
    thr_lowess=0.05,
    pt_avg=1,
    win_size=7,
    pt_smooth_derivative=0,
    multiple_scattering_correction=false,
    method_multiple_scattering_correction="interpolation",
    calibration_OD_curve="NA",
    beta_smoothing_ms=2.0,
    method_peaks_detection="peaks_prominence",
    n_bins=40,
    detect_number_cpd=false,
    fixed_cpd=false,
    penality_CI=8.0,
    size_bootstrap=0.7,
    correction_AIC=true,
    optimizer=BBO_adaptive_de_rand_1_bin_radiuslimited(),
    multistart=false,
    n_restart=50,
    auto_diff_method=nothing,
    cons=nothing,
    opt_params...
    )

This function performs model selection for nonlinear (NL) models while segmenting the time series using change point detection algorithms.

# Arguments:

- `data_testing::Matrix{Float64}`: The dataset with the growth curve, where the first row represents times, and the second row represents the variable to fit (e.g., OD).

- `name_well::String`: Name of the well.

- `label_exp::String`: Label of the experiment.

- `list_of_models::Any`: Array containing functions or strings representing the NL models.

- `list_u0`: Initial guesses for the parameters of each model.

- `n_change_points::Int`: Number of change points to use. The results will vary based on the `type_of_detection` and `fixed_cpd` arguments.

# Key Arguments:

- `lb_param_array::Any = nothing`: Array of lower bounds for the parameters, compatible with the models.

- `ub_param_array::Any = nothing`: Array of upper bounds for the parameters, compatible with the models.

- `type_of_loss::String = "L2_fixed_CI"`: Type of loss function used. Options include `"L2_fixed_CI"`, `"RE"`, `"L2"`, `"L2_derivative"`, and `"blank_weighted_L2"`.

- `method_of_fitting::String = "Normal"`: Method for performing the NL fit. Options are `"Normal"`, `"Bootstrap"`, and `"Morris_sensitivity"`.

- `type_of_detection::String = "sliding_win"`: Change point detection algorithm. Options include `"sliding_win"` and `"lsdd"` (Least Squares Density Difference).

- `type_of_curve::String = "original"`: Curve used for change point detection. Options are `"original"` (original time series) and `"deriv"` (specific growth rate time series).

- `smoothing::Bool = false`: Whether to apply smoothing to the data.

- `nrep::Int = 100`: Number of iterations for Bootstrap or Morris sensitivity analysis.

- `type_of_smoothing::String = "lowess"`: Method for smoothing the data. Options include `"NO"`, `"rolling_avg"`, and `"lowess"`.

- `thr_lowess::Float64 = 0.05`: Threshold parameter for lowess smoothing.

- `pt_avg::Int = 1`: Number of points for generating initial conditions or performing smoothing.

- `win_size::Int = 7`: Size of the window used by the change point detection algorithms.

- `pt_smooth_derivative::Int = 0`: Number of points for evaluating the specific growth rate. Uses interpolation if less than 2; otherwise, applies a sliding window approach.

- `multiple_scattering_correction::Bool = false`: Whether to apply multiple scattering correction using the given calibration curve.

- `method_multiple_scattering_correction::String = "interpolation"`: Method for performing multiple scattering correction. Options are `"interpolation"` and `"exp_fit"`.

- `calibration_OD_curve::String = "NA"`: Path to the CSV file with calibration data for optical density, used if `multiple_scattering_correction` is `true`.

- `beta_smoothing_ms::Float64 = 2.0`: Penalty parameter for AIC (or AICc) evaluation.

- `method_peaks_detection::String = "peaks_prominence"`: Method for peak detection on the dissimilarity curve. Options include `"peaks_prominence"` (orders peaks by prominence) and `"thr_scan"` (uses a threshold).

- `n_bins::Int = 40`: Number of bins used to generate the threshold in `thr_scan` method for peak detection.

- `detect_number_cpd::Bool = false`: Whether to test all possible combinations of change points up to `n_change_points` and return the best based on AICc.

- `fixed_cpd::Bool = false`: If `true`, performs fitting using only the top `n_change_points` and tests a single combination.

- `penality_CI::Float64 = 8.0`: Penalty for enforcing continuity at segment boundaries.

- `size_bootstrap::Float = 0.7`: Fraction of data used in each bootstrap run, applicable if `method_of_fitting` is `"Bootstrap"`.

- `correction_AIC::Bool = true`: Whether to apply finite sample correction to AIC.

- `auto_diff_method::Any = nothing`: Differentiation method for the optimizer, if required.

- `cons::Any = nothing`: Constraints for optimization.

- `multistart::Bool = false`: Whether to use multistart optimization.

- `n_restart::Int = 50`: Number of restarts for multistart optimization, used if `multistart` is `true`.

- `opt_params...`: Additional optional parameters for the optimizer (e.g., `lb` (lower bounds), `ub` (upper bounds), `maxiters`).

# Output:

- A data struct containing:

1. Method string.

2. Matrix with each row containing: `["name of model", "well", "param_1", "param_2", ..., "param_n", "maximum specific gr using NL", "maximum specific gr using data", "objective function value (i.e. loss of the solution)"]`, where `"param_1", "param_2", ..., "param_n"` are the parameters of the selected NL model.

3. The y-coordinates of the NL fit.

4. The x-coordinates of the NL fit.

5. The change point intervals.


"""
function segmentation_NL(
    data_testing::Matrix{Float64}, # dataset first row times second row OD
    name_well::String, # name of the well
    label_exp::String, #label of the experiment
    list_of_models::Any, #  models to use
    list_u0,# initial guess param
    n_change_points::Int;
    lb_param_array::Any=nothing, # lower bound param
    ub_param_array::Any=nothing, # upper bound param   
    type_of_loss="L2_fixed_CI", # type of used loss
    method_of_fitting="NA", # selection of sciml integrator
    type_of_detection="sliding_win",
    type_of_curve="original",
    smoothing=false,
    nrep=100,
    type_of_smoothing="lowess",
    thr_lowess=0.05,
    pt_avg=1,
    win_size=7, # number of the point to generate intial condition
    pt_smooth_derivative=0,
    multiple_scattering_correction=false, # if true uses the given calibration curve to fix the data
    method_multiple_scattering_correction="interpolation",
    calibration_OD_curve="NA", #  the path to calibration curve to fix the data
    beta_smoothing_ms=2.0, #  parameter of the AIC penality
    method_peaks_detection="peaks_prominence",
    n_bins=40,
    detect_number_cpd=false,
    fixed_cpd=false,
    penality_CI=8.0,
    size_bootstrap=0.7,
    correction_AIC=true,
    optimizer=BBO_adaptive_de_rand_1_bin_radiuslimited(),
    multistart=false,
    n_restart=50,
    auto_diff_method=nothing,
    cons=nothing,
    opt_params...
)

    top_aicc = 10^20
    top_param = Vector{Any}
    top_fit = Vector{Any}
    top_time = Vector{Any}
    top_intervals = Vector{Any}


    if multiple_scattering_correction == true
        data_testing = correction_OD_multiple_scattering(data_testing, calibration_OD_curve; method=method_multiple_scattering_correction)
    end
    data_testing_1 = copy(data_testing)
    if smoothing == true
        data_testing_1 = smoothing_data(
            data_testing;
            method=type_of_smoothing,
            pt_avg=pt_avg,
            thr_lowess=thr_lowess
        )
    end




    if detect_number_cpd == true

        list_change_points_dev = cpd_local_detection(
            data_testing_1,
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
            data_testing_1,
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
            data_testing_1,
            n_change_points + 2;
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
            list_of_models, #  models to use
            list_u0,
            cpd_temp;
            lb_param_array =lb_param_array, # lower bound param
            ub_param_array =ub_param_array, # upper bound param
            type_of_loss=type_of_loss, # type of used loss
            optimizer=optimizer, # selection of optimization method
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
            penality_CI=penality_CI,
            correction_AIC=correction_AIC,
            auto_diff_method=auto_diff_method,
            multistart=multistart,
            n_restart=n_restart,
            cons=cons,
            opt_params...
        )


        n_param_full_model = sum([
            length(res_this_combination[1][kk][4:(end-3)]) for
            kk = 1:length(res_this_combination[1])
        ])

        n_param_full_model = n_param_full_model + length(cpd_temp)



        AICc_full_model = AICc_evaluation2(n_param_full_model, beta_smoothing_ms, data_testing[2, :], res_this_combination[end], correction=correction_AIC)

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
    KinBiont_res_segmentation_NL = ("NL_segmentation", top_param, top_fit, top_time, sort(top_intervals))


    return KinBiont_res_segmentation_NL
end

export fit_NL_model
export fit_NL_model_with_sensitivity
export fit_NL_model_bootstrap
export NL_error_blanks
export NL_model_selection
export selection_NL_fixed_interval
export segmentation_NL
