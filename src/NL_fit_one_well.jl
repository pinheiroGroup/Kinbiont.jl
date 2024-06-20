


"""
    fit_NL_model(
    data::Matrix{Float64}, 
    name_well::String,
    label_exp::String, 
    model_function::Any, 
    lb_param::Vector{Float64}, 
    ub_param::Vector{Float64}, 
    u0=lb_param .+ (ub_param .- lb_param)./ 2,
    optmizer=  NLopt.LN_BOBYQA(),    
    pt_avg=1, 
    pt_smooth_derivative=7,
    smoothing=false, 
    type_of_smoothing="rolling_avg",
    type_of_loss="RE", 
    multiple_scattering_correction=false, 
    method_multiple_scattering_correction="interpolation",
    calibration_OD_curve="NA",  
    PopulationSize=300,
    maxiters=2000000,
    abstol=0.00001,
    thr_lowess=0.05,
    penality_CI=3.0
    )


This function fits a nonlinear function to the time series input data of a single well. Method (Xx).

# Arguments:

- `data::Matrix{Float64}`: The dataset with the growth curve, where the first row represents times, and the second row represents the variable to fit (e.g., OD), see documentation.
- `name_well::String`: Name of the well.
- `label_exp::String`: Label of the experiment.
-  `model_function::String`: The model to use, here put the non linear function desired (see documentations for examples) or the string of one of the hard-coded NL models
- `lb_param::Vector{Float64}`: Lower bounds of the model parameters.
- `ub_param::Vector{Float64}`: Upper bounds of the model parameters.

# Key Arguments:

- `param=lb_param .+ (ub_param.-lb_param)./2`: Vector{Float64}. Used as the default initial guess for the model parameters.
- `optmizer=  NLopt.LN_BOBYQA()`: Optimizer from optimizationBBO.
- `type_of_smoothing="rolling_avg"`: String. Method of choice to smooth the data. Options: "NO", "rolling_avg" (rolling average of the data), and "lowess".
- `pt_avg=7`: Int. Size of the rolling average window smoothing. 
- `smoothing=false`: Bool. Options: "true" to smooth the data, or "false" not to.
- `type_of_loss:="RE" `: Type of loss function to be used. Options = "RE" (relative error), "L2" (L2 norm), "L2_derivative" (Xx) and "blank_weighted_L2" (Xx).
- `blank_array=zeros(100)`: Data of all blanks in single array.
- `pt_smoothing_derivative=7`:Int. Number of points for evaluation of specific growth rate. If <2 it uses interpolation algorithm otherwise a sliding window approach.
- `calibration_OD_curve="NA"`: String, The path where the .csv calibration data are located, used only if `multiple_scattering_correction=true`.
- `multiple_scattering_correction=false`: Bool. Options: "true" to perform the multiple scattering correction (requires a callibration curve) or "false" not to. 
- `thr_lowess=0.05`: Float64. Argument of the lowess smoothing.
- `PopulationSize=100`: Size of the population of the optimization (Xx)
- `maxiters=2000000`: stop criterion, the optimization stops when the number of iterations is bigger than `maxiters`.
- `abstol=0.00001`: stop criterion, the optimization stops when the loss is smaller than `abstol`.
- `penality_CI=2.0`: Float64. Used only in segementation to enforce the continuity at the boundary between segments. (We should delete this line!)



# Output (if `results_NL_fit =fit_NL_model(...)`):

- `results_NL_fit[1]`. An array with the following contents: 

`["name of model", "well", "param_1", "param_2",.., "param_n", "maximum specific gr using ode", "maximum specific gr using data", "objective function value (i.e. loss of the solution)"]`, 
where `"param_1","param_2",..,"param_n"` are the parameters of the selected model as in the documentation.

- `results_NL_fit[2]`. The numerical solution of the fitted ODE. (What is the thing with ODE int he output of the NL function?)
"""
function fit_NL_model(data::Matrix{Float64}, # dataset first row times second row OD
    name_well::String, # name of the well
    label_exp::String, #label of the experiment
    model_function::Any, # ode model to use
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
    optmizer=NLopt.LN_BOBYQA(),
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

    loss_function =
        select_loss_function_NL(type_of_loss, data, penality_CI, model_function)


    # Solve the optimization problem
    sol = KimchiSolve_NL(loss_function,
        u0,
        data;
        opt = optmizer,
        auto_diff_method=auto_diff_method,
        multistart = multistart,
        n_restart=n_restart,
        cons=cons,
        opt_params...)

    # evaluate the fitted  model
    fitted_model = model_function(sol, data[1, :])

    sol_fin, index_not_zero = remove_negative_value(fitted_model)

    data_th = transpose(hcat(data[1, index_not_zero], sol_fin))

    max_th_gr = maximum(specific_gr_evaluation(Matrix(data_th), pt_smooth_derivative))

    # max empirical gr
    if length(data[2, :]) > pt_smooth_derivative + 2
        max_em_gr = missing

    else
        max_em_gr = maximum(specific_gr_evaluation(data, pt_smooth_derivative))

    end
    loss_value = sol.objective


    res_param = [[name_well, model_string], [sol[1:end]], [max_th_gr, max_em_gr, loss_value]]


    res_param = reduce(vcat, reduce(vcat, res_param))


    Kimchi_res_one_well = ("NL", res_param, fitted_model, data[1, :])


    return Kimchi_res_one_well
end





"""
    fit_NL_model_with_sensitivity(
    data::Matrix{Float64}, 
    name_well::String,
    label_exp::String, 
    model_function::Any, 
    lb_param::Vector{Float64}, 
    ub_param::Vector{Float64};
    nrep=100,
    u0=lb_param .+ (ub_param .- lb_param) ./ 2,
    optmizer=  NLopt.LN_BOBYQA(),
    pt_avg=1, 
    pt_smooth_derivative=7,
    smoothing=false, 
    type_of_smoothing="rolling_avg",
    type_of_loss="RE", 
    multiple_scattering_correction=false, 
    method_multiple_scattering_correction="interpolation",
    calibration_OD_curve="NA",  
    PopulationSize=300,
    maxiters=2000000,
    abstol=0.00001,
    thr_lowess=0.05,
    penality_CI=3.0,
    )


This function performs the Morris sensitivity analysis for the non-linear fit optimization, which assesses the sensitivity of the fit parameters to variations of the initial guess (suitable for quality checks of nonlinear model fits. See https://docs.sciml.ai/GlobalSensitivity/stable/methods/morris/). 

# Arguments:

- `data::Matrix{Float64}`: The growth curve data. Time values are in the first row and the fit observable (e.g., OD) is in the second row, see documentation.
- `name_well::String`: Name of the well.
- `label_exp::String`: Label of the experiment.
- `model_function::String`: Name of the model of choice (see documentations for examples and list of hard-coded non-linear models).
- `lb_param::Vector{Float64}`: Lower bounds of the model parameters.
- `ub_param::Vector{Float64}`: Upper bounds of the model parameters.



# Key Arguments:
- `nrep=100`.  Number of steps for the Morris sensitivity analysis.
- `param=lb_param .+ (ub_param.-lb_param)./2`:Vector{Float64}. Initial guess for the model parameters.
- `optmizer=  NLopt.LN_BOBYQA()`: Optimizer from optimizationBBO.
- `type_of_smoothing="rolling_avg"`: String. Method of choice to smooth the data. Options: "NO", "rolling_avg" (rolling average of the data), and "lowess".
- `pt_avg=7`: Int. Size of the rolling average window smoothing. 
- `pt_smoothing_derivative=7`: Int. Number of points for evaluation of specific growth rate. If <2 it uses interpolation algorithm otherwise a sliding window approach.
- `smoothing=false`: Bool. Options: "true" to smooth the data, or "false" not to.
- `type_of_loss:="RE" `: Type of loss function to be used. Options = "RE" (relative error), "L2" (L2 norm), "L2_derivative" (Xx) and "blank_weighted_L2" (Xx).
- `blank_array=zeros(100)`: Data of all blanks in a single array.
- `pt_smoothing_derivative=7`: Int. Number of points for evaluation of the specific growth rate. If <2 it uses an interpolation algorithm. Otherwise, it uses a sliding window approach.
- `calibration_OD_curve="NA"`: String. The path where the .csv calibration data are located, used only if `multiple_scattering_correction=true`.
- `multiple_scattering_correction=false`: Bool. Options: "true" to perform the multiple scattering correction (requires a callibration curve) or "false" not to. 
- `method_multiple_scattering_correction="interpolation"`: String. Method of choice to perform the multiple scattering curve inference. Options: '"interpolation"' or '"exp_fit"' (adapted from Meyers, A., Furtmann, C., & Jose, J., *Enzyme and microbial technology*, 118, 1-5., 2018). 
- `thr_lowess=0.05`: Float64. Argument of the lowess smoothing.
- `PopulationSize=100`: Size of the population of the optimization (Xx).
- `maxiters=2000000`: Stop criterion, the optimization stops when the number of iterations is bigger than `maxiters`.
- `abstol=0.00001`: Stop criterion, the optimization stops when the loss is smaller than `abstol`.
- `penality_CI=2.0`: Float64. Used only in segementation to enforce the continuity at the boundary between segments. (We should delete this line!)


# Output (if `results_NL_fit =fit_NL_model_with_sensitivity(...)`):

- `results_NL_fit[1]`. Nn array with the following contents: 
`["name of model", "well", "param_1", "param_2", .., "param_n", "maximum specific gr using ode", "maximum specific gr using data", "objective function value (i.e. loss of the solution)"]`,
 where `"param_1", "param_2", .., "param_n"` are the parameters of the optimal selected model fit as in the documentation. 

 - `results_NL_fit[2]`: Optimal fit parameters. 

 - `results_NL_fit[3]`: the results of the fit for any combination tested. (What does this mean?)

 - `results_NL_fit[4]`: the results of the fit for any combination tested. (????)
"""
function fit_NL_model_with_sensitivity(data::Matrix{Float64}, # dataset first row times second row OD
    name_well::String, # name of the well
    label_exp::String, #label of the experiment
    model_function::Any, # ode model to use
    lb_param::Vector{Float64}, # lower bound param
    ub_param::Vector{Float64}; # upper bound param
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
    optmizer=NLopt.LN_BOBYQA(),
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


    if length(data[2, :]) > pt_smooth_derivative + 2
        max_em_gr = missing

    else
        max_em_gr = maximum(specific_gr_evaluation(data, pt_smooth_derivative))

    end
    fin_param = initialize_df_results_ode_custom(lb_param)
    param_combination =
        generation_of_combination_of_IC_morris(lb_param, ub_param, nrep)

    for i = 1:size(param_combination)[2]
        u0 = param_combination[:, i]



        sol = KimchiSolve_NL(loss_function,
            u0,
            data;
            opt = optmizer,
            auto_diff_method=auto_diff_method,
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


        res_param = [[name_well, model_string], [sol[1:end]], [max_th_gr, max_em_gr, loss_value]]
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


    Kimchi_res_sensitivity_NL = ("NL_sensitivity", best_res_param, best_fitted_model, data[1, :], Matrix(param_combination))

    return Kimchi_res_sensitivity_NL
end





"""
fit_NL_model_bootstrap(
    data::Matrix{Float64}, 
    name_well::String,
    label_exp::String, 
    model_function::Any, 
    lb_param::Vector{Float64}, 
    ub_param::Vector{Float64};
    nrep=100,
    u0=lb_param .+ (ub_param .- lb_param) ./ 2,
    optmizer=  NLopt.LN_BOBYQA(),
    pt_avg=1, 
    size_bootstrap=0.7,
    pt_smooth_derivative=7,
    smoothing=false, 
    type_of_smoothing="rolling_avg",
    type_of_loss="RE", 
    multiple_scattering_correction=false, 
    method_multiple_scattering_correction="interpolation",
    calibration_OD_curve="NA",  
    PopulationSize=300,
    maxiters=2000000,
    abstol=0.00001,
    thr_lowess=0.05,
    penality_CI=3.0,
    )


This function performs NL fitting. It perform nrep iterations of Bootstrap to evaluate the confidence intervals and  avoid bad initializations.

- `data::Matrix{Float64}`: The dataset with the growth curve, where the first row represents times, and the second row represents the variable to fit (e.g., OD), see documentation.
- `name_well::String`: Name of the well.
- `label_exp::String`: Label of the experiment.
-  `model_function::String`: The model to use, here put the non linear function desired (see documentations for examples) or the string of one of the hard-coded NL models
- `lb_param::Vector{Float64}`: Lower bounds of the model parameters.
- `ub_param::Vector{Float64}`: Upper bounds of the model parameters.

# Key Arguments:
-  `size_bootstrap=0.7`: Float, the fraction of data used each Bootstrap run
- `nrep=100`. Number of MCMC steps.
- `param= lb_param .+ (ub_param.-lb_param)./2`:Vector{Float64}, Initial guess for the model parameters.
- `optmizer =     NLopt.LN_BOBYQA()` optimizer from optimizationBBO.
- `type_of_smoothing="rolling_avg"`: String, How to smooth the data, options: "NO" , "rolling avg" rolling average of the data, and "lowess".
- `pt_avg=7`: Number of points to generate the initial condition or do the rolling avg smoothing.
- `smoothing=false`: Whether to apply smoothing to the data or not.
- `type_of_loss:="RE" `: Type of loss function to be used. (options= "RE", "L2", "L2_derivative" and "blank_weighted_L2").
- `blank_array=zeros(100)`: Data of all blanks in single array.
- `pt_smoothing_derivative=7`:Int,  Number of points for evaluation of specific growth rate. If <2 it uses interpolation algorithm otherwise a sliding window approach.
- `calibration_OD_curve="NA"`: String, The path where the .csv calibration data are located, used only if `multiple_scattering_correction=true`.
- `multiple_scattering_correction=false`: Bool, if true uses the given calibration curve to correct the data for muliple scattering.
- `method_multiple_scattering_correction="interpolation"`: String, How perform the inference of multiple scattering curve, options: "interpolation" or   "exp_fit" it uses an exponential fit from "Direct optical density determination of bacterial cultures in microplates for high-throughput screening applications"
-  `thr_lowess=0.05`: Float64 keyword argument of lowees smoothing
- ` PopulationSize =100`: Size of the population of the optimization
- ` maxiters=2000000`: stop criterion, the optimization is stopped when the number of iterations is bigger than `maxiters`
- `abstol = 0.00001`: stop criterion, the optimization is stopped when the loss is lesser than `abstol`
- `penality_CI=2.0`, used only in segementation to force the optimization to respect continuty on bonduar

# Output (if `results_NL_fit =fit_NL_model_bootstrap(...)`:

- `results_NL_fit[1]` an array with the following contents: `["name of model", "well", "param_1","param_2",..,"param_n","maximum specific gr using ode","maximum specific gr using data", "objective function value (i.e. loss of the solution)"]` where ' "param_1","param_2",..,"param_n" ' are the parameter of the selected model as in the documentation. This for the fit with less loss.
- `results_NL_fit[2]`: the array of the best fit 
- `results_NL_fit[3]`:parameters 
-`results_NL_fit[4]`: parameters
-`results_NL_fit[5]`: mean best parameters
-`results_NL_fit[6]`: std best parameters
-`results_NL_fit[7]`:CI lower bound
-`results_NL_fit[8]`:CI upper bound

"""
function fit_NL_model_bootstrap(data::Matrix{Float64}, # dataset first row times second row OD
    name_well::String, # name of the well
    label_exp::String, #label of the experiment
    model_function::Any, # ode model to use
    u0;# initial guess param
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
    optmizer=NLopt.LN_BOBYQA(),
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
            select_loss_function_NL(type_of_loss, data_to_fit, penality_CI, model_function)

        sol = KimchiSolve_NL(loss_function,
            u0,
            data_to_fit;
           opt =  optmizer,
            auto_diff_method=auto_diff_method,
            cons=cons,
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


        res_param = [[name_well, model_string], [sol[1:end]], [max_th_gr, max_em_gr, loss_value]]
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

    Kimchi_res_bootstrap_NL = ("NL_bootstrap", best_res_param, best_fitted_model, data[1, :], fin_param, new_param_fin, mean_param, sd_param, CI_param_low, CI_param_up)


    return Kimchi_res_bootstrap_NL
end




"""
    NL_error_blanks(data::Matrix{Float64}, 
    name_well::String, 
    label_exp::String, 
    model_function::Any,
    lb_param::Vector{Float64},
    ub_param::Vector{Float64},
    blank_array::Vector{Float64}; 
    nrep=100,
    u0=lb_param .+ (ub_param .- lb_param) ./ 2,
optmizer=  NLopt.LN_BOBYQA(),    pt_avg=1, 
    pt_smooth_derivative=7,
    smoothing=false,
    type_of_smoothing="rolling_avg",
    type_of_loss="RE", 
    multiple_scattering_correction=false,
    method_multiple_scattering_correction="interpolation",
    calibration_OD_curve="NA",  
    PopulationSize=300,
    maxiters=2000000,
    abstol=0.00001,
    thr_lowess=0.05,
    write_res=false,
    penality_CI=3.0
    )



This function performs NL fitting. It perform nrep iterations to estimate the posterior distribuition of parameters fitting. It uses the blank distribution as noise.

- `data::Matrix{Float64}`: The dataset with the growth curve, where the first row represents times, and the second row represents the variable to fit (e.g., OD), see documentation.
- `name_well::String`: Name of the well.
- `label_exp::String`: Label of the experiment.
-  `model_function::String`: The model to use, here put the non linear function desired (see documentations for examples) or the string of one of the hard-coded NL models
- `lb_param::Vector{Float64}`: Lower bounds of the model parameters.
- `ub_param::Vector{Float64}`: Upper bounds of the model parameters.

# Key Arguments:
- `nrep=100`. Number of MCMC steps.
- `param= lb_param .+ (ub_param.-lb_param)./2`:Vector{Float64}, Initial guess for the model parameters.
- `optmizer =     NLopt.LN_BOBYQA()` optimizer from optimizationBBO.
- `type_of_smoothing="rolling_avg"`: String, How to smooth the data, options: "NO" , "rolling avg" rolling average of the data, and "lowess".
- `pt_avg=7`: Number of points to generate the initial condition or do the rolling avg smoothing.
- `smoothing=false`: Whether to apply smoothing to the data or not.
- `type_of_loss:="RE" `: Type of loss function to be used. (options= "RE", "L2", "L2_derivative" and "blank_weighted_L2").
- `blank_array=zeros(100)`: Data of all blanks in single array.
- `pt_smoothing_derivative=7`:Int,  Number of points for evaluation of specific growth rate. If <2 it uses interpolation algorithm otherwise a sliding window approach.
- `calibration_OD_curve="NA"`: String, The path where the .csv calibration data are located, used only if `multiple_scattering_correction=true`.
- `multiple_scattering_correction=false`: Bool, if true uses the given calibration curve to correct the data for muliple scattering.
- `method_multiple_scattering_correction="interpolation"`: String, How perform the inference of multiple scattering curve, options: "interpolation" or   "exp_fit" it uses an exponential fit from "Direct optical density determination of bacterial cultures in microplates for high-throughput screening applications"
-  `thr_lowess=0.05`: Float64 keyword argument of lowees smoothing
- ` PopulationSize =100`: Size of the population of the optimization
- ` maxiters=2000000`: stop criterion, the optimization is stopped when the number of iterations is bigger than `maxiters`
- `abstol = 0.00001`: stop criterion, the optimization is stopped when the loss is lesser than `abstol`
- `penality_CI=2.0`, used only in segementation to force the optimization to respect continuty on bonduar

# Output (if `results_NL_fit =NL_error_blanks(...)`:

- `results_NL_fit[1]` an array with the following contents: `["name of model", "well", "param_1","param_2",..,"param_n","maximum specific gr using ode","maximum specific gr using data", "objective function value (i.e. loss of the solution)"]` where ' "param_1","param_2",..,"param_n" ' are the parameter of the selected model as in the documentation. This for the fit with less loss.
- `results_NL_fit[2]`: the array of the best fit 
- `results_NL_fit[3]`:parameters 
-`results_NL_fit[4]`: parameters
-`results_NL_fit[5]`: mean best parameters
-`results_NL_fit[6]`: std best parameters
-`results_NL_fit[7]`:CI lower bound
-`results_NL_fit[8]`:CI upper bound

"""
function NL_error_blanks(data::Matrix{Float64}, # dataset first row times second row OD
    name_well::String, # name of the well
    label_exp::String, #label of the experiment
    model_function::Any, # ode model to use
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
    optmizer=NLopt.LN_BOBYQA(),
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
    fin_param = initialize_df_results_ode_custom(lb_param)


    for i = 1:nrep

        data[2, :] = data[2, :] .+ rand(Normal(0.0, std(blank_array)), length(data[2, :]))
        sol_fin, index_not_zero = remove_negative_value(data[2, :])

        data = transpose(hcat(data[1, index_not_zero], sol_fin))




        loss_function =
            select_loss_function_NL(type_of_loss, data, penality_CI, model_function)

        sol = KimchiSolve_NL(loss_function,
            u0,
            data;
           opt = optmizer,
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


        res_param = [[name_well, model_string], [sol[1:end]], [max_th_gr, max_em_gr, loss_value]]
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
    NL_model_selection(data::Matrix{Float64}, 
    name_well::String,
    label_exp::String, 
    list_model_function::Any, 
    list_lb_param::Any, 
    list_ub_param::Any; 
    method_of_fitting="MCMC",
    nrep=100,
    list_u0=list_lb_param .+ (list_ub_param .- list_lb_param) ./ 2,
optmizer=  NLopt.LN_BOBYQA(),    size_bootstrap=0.7,
    pt_avg=1,
    pt_smooth_derivative=7,
    smoothing=false, 
    type_of_smoothing="rolling_avg",
    type_of_loss="RE", 
    multiple_scattering_correction=false, 
    method_multiple_scattering_correction="interpolation",
    calibration_OD_curve="NA",  
    PopulationSize=300,
    maxiters=2000000,
    abstol=0.00001,
    thr_lowess=0.05,
    write_res=false,
    beta_param=2.0,
    penality_CI=8.0,
    correction_AIC=false,
    )



This function performs NL model selection of an array of NL models, it uses AIC or AICc depending on user inputs. It perform nrep iterations to estimate the posterior distribuition of parameters fitting. It uses the blank distribution as noise.

- `data::Matrix{Float64}`: The dataset with the growth curve, where the first row represents times, and the second row represents the variable to fit (e.g., OD), see documentation.
- `name_well::String`: Name of the well.
- `label_exp::String`: Label of the experiment.
-  `list_model_function::Any`: Array containing functions or strings of the NL models
-  `list_lb_param::Any`:Array of Lower bounds for the parameters (compatible with the models).
-  `list_ub_param::Any`:Array of Upper bounds for the parameters (compatible with the models).

# Key Arguments:
- `method_of_fitting="MCMC"`: String, how perform the NL fit. Options "MCMC","Bootstrap","Normal", and "Morris_sensitivity"
- `nrep=100`. Number of MCMC steps.
- `param= lb_param .+ (ub_param.-lb_param)./2`:Vector{Float64}, Initial guess for the model parameters.
- `optmizer =     NLopt.LN_BOBYQA()` optimizer from optimizationBBO.
- `type_of_smoothing="rolling_avg"`: String, How to smooth the data, options: "NO" , "rolling avg" rolling average of the data, and "lowess".
- `pt_avg=7`: Number of points to generate the initial condition or do the rolling avg smoothing.
- `smoothing=false`: Whether to apply smoothing to the data or not.
- `type_of_loss:="RE" `: Type of loss function to be used. (options= "RE", "L2", "L2_derivative" and "blank_weighted_L2").
- `blank_array=zeros(100)`: Data of all blanks in single array.
- `pt_smoothing_derivative=7`:Int,  Number of points for evaluation of specific growth rate. If <2 it uses interpolation algorithm otherwise a sliding window approach.
- `calibration_OD_curve="NA"`: String, The path where the .csv calibration data are located, used only if `multiple_scattering_correction=true`.
- `multiple_scattering_correction=false`: Bool, if true uses the given calibration curve to correct the data for muliple scattering.
- `method_multiple_scattering_correction="interpolation"`: String, How perform the inference of multiple scattering curve, options: "interpolation" or   "exp_fit" it uses an exponential fit from "Direct optical density determination of bacterial cultures in microplates for high-throughput screening applications"
-  `thr_lowess=0.05`: Float64 keyword argument of lowees smoothing
- ` PopulationSize =100`: Size of the population of the optimization
- ` maxiters=2000000`: stop criterion, the optimization is stopped when the number of iterations is bigger than `maxiters`
- `abstol = 0.00001`: stop criterion, the optimization is stopped when the loss is lesser than `abstol`
- `penality_CI=2.0`, used only in segementation to force the optimization to respect continuty on bonduar
-  `correction_AIC=true`: Bool, do finite samples correction of AIC.
-  `beta_param=2.0` penality  parameters for AIC (or AICc) evaluation.
-  `size_bootstrap=0.7`: Float, the fraction of data used each Bootstrap run. Used only if method is "Bootstrap"

# Output (if `results_NL_fit =NL_model_selection(...)`:

- `results_NL_fit[1]` an array with the following contents: `["name of model", "well", "param_1","param_2",..,"param_n","maximum specific gr using ode","maximum specific gr using data", "objective function value (i.e. loss of the solution)"]` where ' "param_1","param_2",..,"param_n" ' are the parameter of the selected model as in the documentation. This for the fit with less loss.
- `results_NL_fit[2]`: an array with the following contents: `["name of model", "well", "param_1","param_2",..,"param_n","maximum specific gr using ode","maximum specific gr using data", "objective function value (i.e. loss of the solution)"]` where ' "param_1","param_2",..,"param_n" ' are the parameter of the selected model as in the documentation. This for the fit with less loss.
- `results_NL_fit[3]`: the array of best fit 
-`results_NL_fit[4]`: scores of the models
-`results_NL_fit[5]`: the loss of the best loss

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
    beta_param=2.0,
    penality_CI=8.0,
    correction_AIC=false,
    optmizer=NLopt.LN_BOBYQA(),
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
    for mm in 1:eachindex(list_model_function)[end]

        model_to_test = list_model_function[mm]

        if !isnothing(lb_param_array)
            temp_param_lb = lb_param_array[mm]
            temp_param_ub = ub_param_array[mm]

            opt_params = (opt_params...,
                lb=temp_param_lb,
                ub=temp_param_ub,
                )
        end

        u0 = list_u0[mm]

        if method_of_fitting == "Bootstrap"



            temp_res = fit_NL_model_bootstrap(data, # dataset first row times second row OD
                name_well, # name of the well
                label_exp, #label of the experiment
                model_to_test, # ode model to use
                u0;
                lb_param=temp_param_lb, # lower bound param
                ub_param=temp_param_ub, # upper bound param
                nrep=nrep,
                optmizer=optmizer,
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

            n_param = length(u0)
            temp_AIC = AICc_evaluation2(n_param, beta_param, data[2, :], temp_res[2][end], correction=correction_AIC)

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
                model_to_test, # ode model to use
                temp_param_lb, # lower bound param
                temp_param_ub; # upper bound param
                nrep=nrep,
                optmizer=optmizer,
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

            n_param = length(temp_param_lb)

            temp_AIC = AICc_evaluation2(n_param, beta_param, data[2, :], temp_res[2][end], correction=correction_AIC)

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
                model_to_test, # ode model to use
                u0;
                optmizer=optmizer,
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
                opt_params...
            )


            n_param = length(u0)

            temp_AIC = AICc_evaluation2(n_param, beta_param, data[2, :], temp_res[2][end], correction=correction_AIC)

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
    Kimchi_res_NL_model_selection = ("NL_model_selection", top_model, top_fitted_sol, data[1, :], score_res, top_loss)
    return Kimchi_res_NL_model_selection
end

"""
    selection_NL_fixed_interval(
    data_testing::Matrix{Float64},
    name_well::String, 
    label_exp::String, 
    list_of_models::Vector{String}, 
    list_lb_param::Any, 
    list_ub_param::Any, 
    intervals_changepoints::Any;
    list_u0=list_lb_param .+ (list_ub_param .- list_lb_param) ./ 2,
    type_of_loss="L2", 
    optmizer=  NLopt.LN_BOBYQA(), 
    method_of_fitting="MCMC",
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
    PopulationSize=300,
    maxiters=2000000,
    abstol=0.000000001,
    penality_CI=8.0,
    correction_AIC=true,
    )



    This function performs a fitting of a segmented NL on one curve. For this function the user must supply the change points.

- `data::Matrix{Float64}`: The dataset with the growth curve, where the first row represents times, and the second row represents the variable to fit (e.g., OD), see documentation.
- `name_well::String`: Name of the well.
- `label_exp::String`: Label of the experiment.
-  `list_model_function::Any`: Array containing functions or strings of the NL models
-  `list_lb_param::Any`:Array of Lower bounds for the parameters (compatible with the models).
-  `list_ub_param::Any`:Array of Upper bounds for the parameters (compatible with the models).
- `models_list::Vector{String}`: Array of  models to evaluate.
- `intervals_changepoints::Any`: the array containings the change point list, e.g., [0.0 10.0 30.0] 


# Key Arguments:
- `method_of_fitting="MCMC"`: String, how perform the NL fit. Options "MCMC","Bootstrap","Normal", and "Morris_sensitivity"
- `nrep=100`. Number of MCMC steps.
- `param= lb_param .+ (ub_param.-lb_param)./2`:Vector{Float64}, Initial guess for the model parameters.
- `optmizer =     NLopt.LN_BOBYQA()` optimizer from optimizationBBO.
- `type_of_smoothing="rolling_avg"`: String, How to smooth the data, options: "NO" , "rolling avg" rolling average of the data, and "lowess".
- `pt_avg=7`: Number of points to generate the initial condition or do the rolling avg smoothing.
- `smoothing=false`: Whether to apply smoothing to the data or not.
- `type_of_loss:="RE" `: Type of loss function to be used. (options= "RE", "L2", "L2_derivative" and "blank_weighted_L2").
- `blank_array=zeros(100)`: Data of all blanks in single array.
- `pt_smoothing_derivative=7`:Int,  Number of points for evaluation of specific growth rate. If <2 it uses interpolation algorithm otherwise a sliding window approach.
- `calibration_OD_curve="NA"`: String, The path where the .csv calibration data are located, used only if `multiple_scattering_correction=true`.
- `multiple_scattering_correction=false`: Bool, if true uses the given calibration curve to correct the data for muliple scattering.
- `method_multiple_scattering_correction="interpolation"`: String, How perform the inference of multiple scattering curve, options: "interpolation" or   "exp_fit" it uses an exponential fit from "Direct optical density determination of bacterial cultures in microplates for high-throughput screening applications"
-  `thr_lowess=0.05`: Float64 keyword argument of lowees smoothing
- ` PopulationSize =100`: Size of the population of the optimization
- ` maxiters=2000000`: stop criterion, the optimization is stopped when the number of iterations is bigger than `maxiters`
- `abstol = 0.00001`: stop criterion, the optimization is stopped when the loss is lesser than `abstol`
- `penality_CI=2.0`, used only in segementation to force the optimization to respect continuty on bonduar
-  `correction_AIC=true`: Bool, do finite samples correction of AIC.
-  `beta_param=2.0` penality  parameters for AIC (or AICc) evaluation.
-  `size_bootstrap=0.7`: Float, the fraction of data used each Bootstrap run. Used only if method is "Bootstrap"

# Output (if `results_NL_fit =selection_NL_fixed_interval(...)`:


- `results_NL_fit[1]`. Parameters of each segment
- `results_NL_fit[2]`. The numerical solutions of the fit
- `results_NL_fit[3]`.  The time of the numerical solutions of the fit
- `results_NL_fit[4]`. the loss of the best loss

"""
function selection_NL_fixed_interval(
    data_testing::Matrix{Float64}, # dataset first row times second row OD
    name_well::String, # name of the well
    label_exp::String, #label of the experiment
    list_of_models::Vector{String}, # ode models to use
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
    optmizer=NLopt.LN_BOBYQA(),
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
            list_of_models, # ode model to use
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
            beta_param=beta_smoothing_ms,
            penality_CI=penality_CI,
            correction_AIC=correction_AIC,
            optmizer=optmizer,
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
    selection_NL_max_change_points(
    data_testing::Matrix{Float64},
    name_well::String, 
    label_exp::String,
    list_of_models::Any, 
    list_lb_param::Any, 
    list_ub_param::Any, 
    n_change_points::Int;
    list_u0=list_lb_param .+ (list_ub_param .- list_lb_param) ./ 2,
    type_of_loss="L2_fixed_CI", 
optmizer=  NLopt.LN_BOBYQA(),    method_of_fitting="MCMC", 
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
    PopulationSize=300,
    maxiters=2000000,
    abstol=0.000000001,
    detect_number_cpd=false,
    fixed_cpd=false,
    penality_CI=8.0,
    size_bootstrap=0.7,
    correction_AIC=true
    )


This function performs model selection for NL models while segmenting the time series in various part using change points detection algorithm.

- `data::Matrix{Float64}`: The dataset with the growth curve, where the first row represents times, and the second row represents the variable to fit (e.g., OD), see documentation.
- `name_well::String`: Name of the well.
- `label_exp::String`: Label of the experiment.
-  `list_model_function::Any`: Array containing functions or strings of the NL models
-  `list_lb_param::Any`:Array of Lower bounds for the parameters (compatible with the models).
-  `list_ub_param::Any`:Array of Upper bounds for the parameters (compatible with the models).
- `models_list::Vector{String}`: Array of  models to evaluate.
- `intervals_changepoints::Any`: the array containings the change point list, e.g., [0.0 10.0 30.0] 


# Key Arguments:
- `method_of_fitting="MCMC"`: String, how perform the NL fit. Options "MCMC","Bootstrap","Normal", and "Morris_sensitivity"
- `nrep=100`. Number of MCMC steps.
- `param= lb_param .+ (ub_param.-lb_param)./2`:Vector{Float64}, Initial guess for the model parameters.
- `optmizer =     NLopt.LN_BOBYQA()` optimizer from optimizationBBO.
- `type_of_smoothing="rolling_avg"`: String, How to smooth the data, options: "NO" , "rolling avg" rolling average of the data, and "lowess".
- `pt_avg=7`: Number of points to generate the initial condition or do the rolling avg smoothing.
- `smoothing=false`: Whether to apply smoothing to the data or not.
- `type_of_loss:="RE" `: Type of loss function to be used. (options= "RE", "L2", "L2_derivative" and "blank_weighted_L2").
- `blank_array=zeros(100)`: Data of all blanks in single array.
- `pt_smoothing_derivative=7`:Int,  Number of points for evaluation of specific growth rate. If <2 it uses interpolation algorithm otherwise a sliding window approach.
- `calibration_OD_curve="NA"`: String, The path where the .csv calibration data are located, used only if `multiple_scattering_correction=true`.
- `multiple_scattering_correction=false`: Bool, if true uses the given calibration curve to correct the data for muliple scattering.
- `method_multiple_scattering_correction="interpolation"`: String, How perform the inference of multiple scattering curve, options: "interpolation" or   "exp_fit" it uses an exponential fit from "Direct optical density determination of bacterial cultures in microplates for high-throughput screening applications"
-  `thr_lowess=0.05`: Float64 keyword argument of lowees smoothing
- ` PopulationSize =100`: Size of the population of the optimization
- ` maxiters=2000000`: stop criterion, the optimization is stopped when the number of iterations is bigger than `maxiters`
- `abstol = 0.00001`: stop criterion, the optimization is stopped when the loss is lesser than `abstol`
- `penality_CI=2.0`, used only in segementation to force the optimization to respect continuty on bonduar
-  `correction_AIC=true`: Bool, do finite samples correction of AIC.
-  `beta_param=2.0` penality  parameters for AIC (or AICc) evaluation.
-  `size_bootstrap=0.7`: Float, the fraction of data used each Bootstrap run. Used only if method is "Bootstrap"

# Output (if `results_NL_fit =selection_NL_fixed_interval(...)`:

- `results_NL_fit[1]` an array with the following the parameters of each segment
- `results_NL_fit[2]`: the list of used change points
- `results_NL_fit[3]`:the numerical best solution
-`results_NL_fit[4]`: the time of the numerical best solution

"""
function selection_NL_max_change_points(
    data_testing::Matrix{Float64}, # dataset first row times second row OD
    name_well::String, # name of the well
    label_exp::String, #label of the experiment
    list_of_models::Any, # ode models to use
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
    optmizer=NLopt.LN_BOBYQA(),
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
            list_of_models, # ode models to use
            list_u0,
            cpd_temp;
            lb_param_array =lb_param_array, # lower bound param
            ub_param_array =ub_param_array, # upper bound param
            type_of_loss=type_of_loss, # type of used loss
            optmizer=optmizer, # selection of optimization method
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
            length(res_this_combination[1][kk][3:(end-4)]) for
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
    Kimchi_res_segmentation_NL = ("NL_segmentation", top_param, top_fit, top_time, sort(top_intervals))


    return Kimchi_res_segmentation_NL
end

export fit_NL_model
export fit_NL_model_with_sensitivity
export fit_NL_model_MCMC_intialization
export fit_NL_model_bootstrap
export NL_error_blanks
export NL_model_selection
export selection_NL_fixed_interval
export selection_NL_max_change_points
