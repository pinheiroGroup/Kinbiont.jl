using Distributions
using StatsBase
using Optimization
using OptimizationBBO
using OptimizationMultistartOptimization
#######################################################################
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

This function fits a nonlinear (NL) function to the time series input data of a single well. 

# Arguments:

- `data::Matrix{Float64}`: The dataset with the growth curve. The first row represents times, and the second row represents the variable to fit (e.g., optical density). Ensure the data is formatted correctly as a 2-row matrix where each column corresponds to a time-point and its associated value.

- `name_well::String`: The name of the well for which the model is being fitted.

- `label_exp::String`: Label for the experiment, used for identification in outputs and results.

- `model_function::Any`: The nonlinear model to use for fitting. This can be a function or a string representing one of the predefined NL models.

- `u0::Any`: Initial guess for the model parameters. This should match the expected input format for the model being used.

# Key Arguments:

- `optimizer::Any = BBO_adaptive_de_rand_1_bin_radiuslimited()`: Optimizer used for parameter estimation. Default is `BBO_adaptive_de_rand_1_bin_radiuslimited()` from `optimizationBBO`.

- `type_of_smoothing::String = "rolling_avg"`: Method for smoothing the data. Options include `"NO"` (no smoothing), `"rolling_avg"` (rolling average), and `"lowess"` (locally weighted scatterplot smoothing).

- `pt_avg::Int = 7`: Number of points used for the rolling average window if `type_of_smoothing` is `"rolling_avg"`.

- `smoothing::Bool = false`: Whether to apply smoothing to the data. Set to `true` to enable smoothing.

- `type_of_loss::String = "RE"`: Type of loss function to be used for fitting. Options include `"RE"` (relative error), `"L2"` (L2 norm), `"L2_derivative"`, and `"blank_weighted_L2"`.

- `blank_array::Vector{Float64} = zeros(100)`: Array containing blank data values, used for blank correction if needed.

- `pt_smooth_derivative::Int = 7`: Number of points for evaluating the specific growth rate. If less than 2, an interpolation algorithm is used; otherwise, a sliding window approach is used.

- `calibration_OD_curve::String = "NA"`: Path to a CSV file with calibration data for optical density, used only if `multiple_scattering_correction` is true.

- `multiple_scattering_correction::Bool = false`: Whether to apply multiple scattering correction. Set to `true` to use the calibration curve provided in `calibration_OD_curve`.

- `thr_lowess::Float64 = 0.05`: Threshold parameter for the lowess smoothing method.

- `penality_CI::Float64 = 3.0`: Penalty parameter for ensuring continuity at boundaries in segmentation. (Note: Consider removing this parameter if it is not applicable.)

- `auto_diff_method::Any = nothing`: Differentiation method for the optimizer, to be specified if required.

- `cons::Any = nothing`: Constraints for the optimization process.

- `multistart::Bool = false`: Whether to use multistart optimization. Set to `true` to enable multistart.

- `n_restart::Int = 50`: Number of restarts for multistart optimization, used if `multistart` is true.

- `opt_params...`: Additional optional parameters for the optimizer, such as `lb` (lower bounds), `ub` (upper bounds), and `maxiters` (maximum iterations).

# Output:

- **Method String**: A string indicating the method used for the fitting process.

- **Results Array**: An array where each entry contains:
  - `"name of model"`: The name of the model used.
  - `"well"`: The name of the well.
  - `"param_1", "param_2", ..., "param_n"`: Parameters of the fitted NL model.
  - `"maximum specific gr using NL"`: Maximum specific growth rate obtained using the NL model.
  - `"maximum specific gr using data"`: Maximum specific growth rate obtained from the data.
  - `"objective function value"`: Value of the objective function (i.e., the loss function value).

- **Fitted NL Data**: Numerical solution of the fitted NL model.

- **Fitted Time Coordinates**: Time coordinates corresponding to the fitted NL data.


"""
function fitting_one_well_Log_Lin(
    data::Matrix{Float64}, # dataset first row times second row OD
    name_well::String, # name of the well
    label_exp::String; #label of the experiment
    type_of_smoothing="rolling_avg", # option, NO, gaussian, rolling avg
    pt_avg=7, # number of the point for rolling avg not used in the other cases
    pt_smoothing_derivative=7, # number of poits to smooth the derivative
    pt_min_size_of_win=7, # minimum size of the exp windows in number of smooted points
    type_of_win="maximum", # how the exp. phase win is selected, "maximum" of "global_thr"
    threshold_of_exp=0.9, # threshold of growth rate in quantile to define the exp windows
    multiple_scattering_correction=false, # if true uses the given calibration curve to fix the data
    method_multiple_scattering_correction="interpolation",
    calibration_OD_curve="NA", #  the path to calibration curve to fix the data
    thr_lowess=0.05, # keyword argument of lowees smoothing
    start_exp_win_thr=0.05, # minimum value to consider the start of exp window
)
    if multiple_scattering_correction == true

        data = correction_OD_multiple_scattering(data, calibration_OD_curve; method=method_multiple_scattering_correction)
    end


    data_smooted = smoothing_data(
        data;
        method=type_of_smoothing,
        pt_avg=pt_avg,
        thr_lowess=thr_lowess
    )


    # local fitting generating of specific growth rate (derivative)
    specific_gr = specific_gr_evaluation(data_smooted, pt_smoothing_derivative)
    specific_gr_times = [
        (data_smooted[1, r] + data_smooted[1, (r+pt_smoothing_derivative)]) / 2 for
        r = 1:1:(eachindex(data_smooted[2, :])[end].-pt_smoothing_derivative)
    ]

    # selecting the max
    gr_max = maximum(specific_gr)
    index_of_max = findall(x -> x == gr_max, specific_gr)[1]

    # selecting the windows
    # lower threshold of the exp phase using quantile
    lb_of_distib = quantile(specific_gr, threshold_of_exp)

    # searching for t_start_exp
    t_start = 0.0

    if type_of_win == "maximum"
        for yy = 1:(index_of_max-2)
            if specific_gr[index_of_max-yy] >= lb_of_distib
                t_start = copy(specific_gr_times[index_of_max-yy])
            else
                if specific_gr[(index_of_max-yy-1)] < lb_of_distib
                    break
                end

            end
        end

        # index of start
        index_of_t_start = findfirst(x -> x > t_start, data_smooted[1, :])[1]

        # searching t_end of the exp phase
        t_end = specific_gr_times[end]

        for yy = index_of_max:(eachindex(specific_gr)[end]-1)
            if specific_gr[yy] >= lb_of_distib
                t_end = copy(specific_gr_times[yy])
            else
                if specific_gr[(yy+1)] < lb_of_distib
                    break
                end
            end
        end

        index_of_t_end = findfirst(x -> x > t_end, data_smooted[1, :])[1]
    end

    # selection of exp win with a global thr on the growht rate
    if type_of_win == "global_thr"

        index_od_over_thr = findfirst(data_smooted[2, :] .> start_exp_win_thr)

        if isnothing(index_od_over_thr) == false && length(specific_gr[index_od_over_thr:end]) > 2

            lb_of_distib = quantile(specific_gr[index_od_over_thr:end], threshold_of_exp)

            index_of_max = argmax(specific_gr[index_od_over_thr:end])[1] + index_od_over_thr - 1

            index_gr_max = findlast(x -> x > lb_of_distib, specific_gr[index_od_over_thr:end])[1] + index_od_over_thr - 1

            index_gr_min = findfirst(x -> x > lb_of_distib, specific_gr[index_od_over_thr:index_of_max])[1] + index_od_over_thr - 1




            t_start = specific_gr_times[index_gr_min]
            t_end = specific_gr_times[index_gr_max]

            index_of_t_start = findfirst(x -> x > t_start, data_smooted[1, :])[1]
            index_of_t_end = findall(x -> x > t_end, data_smooted[1, :])[1]


        else
            # the conditions are not satisfied i fix t_end and  t_start to be discarded later and have all results to missing
            t_end = 100
            t_start = 200
            index_of_t_start = 1
            index_of_t_end = index_of_t_start + 2 * pt_min_size_of_win

        end


    end

    # checking the minimum size of the window before fitting
    if (index_of_t_end - index_of_t_start) < pt_min_size_of_win
        index_of_t_start = convert(Int, index_of_max - floor(pt_min_size_of_win / 2))
        index_of_t_end = convert(Int, index_of_max + floor(pt_min_size_of_win / 2))

        if index_of_t_start < 1
            index_of_t_start =1 
        end

        if index_of_t_end > length(data_smooted[1, :])
            index_of_t_end = length(data_smooted[1, :]) 
        end
    end

    if t_end > t_start
        # fitting data
        data_to_fit_times = data_smooted[1, index_of_t_start:index_of_t_end]
        data_to_fit_values = log.(data_smooted[2, index_of_t_start:index_of_t_end])

        N = length(data_to_fit_times)
        M = [ones(N) data_to_fit_times]
        (coeff_1, coeff_2) = M \ data_to_fit_values
        mean_x = mean(data_to_fit_times)

        sigma_a = sigma_b = r = zeros(N)
        Theoretical_fitting = coeff_1 .+ data_to_fit_times .* coeff_2
        
        Cantrell_errors = sqrt(sum((data_to_fit_values - coeff_2 * data_to_fit_times .- coeff_1) .^ 2) / (N - 2))  # goodness of fit
        sigma_b = sqrt(1 / sum((data_to_fit_times .- mean_x) .^ 2))
        sigma_a = Cantrell_errors * sqrt(1 / N + mean_x^2 * sigma_b^2)
        sigma_b *= Cantrell_errors
        # Pearson's correlation coefficient
        rho = cov(data_to_fit_times, data_to_fit_values) / sqrt(var(data_to_fit_times) * var(data_to_fit_values))
        d = TDist(N - 2)     # t-Student distribution with N-2 degrees of freedom
        cf = quantile(d, 0.975)  # correction factor for 95% confidence intervals (two-tailed distribution)
        confidence_band = cf * Cantrell_errors * sqrt.(1 / N .+ (data_to_fit_times .- mean(data_to_fit_times)) .^ 2 / var(data_to_fit_times) / (N - 1))


        # storing results

        results_lin_log_fit = [
            label_exp,
            name_well,
            data_to_fit_times[1],
            data_to_fit_times[end],
            specific_gr_times[index_of_max],
            gr_max,
            coeff_2,
            sigma_b,
            log(2) / (coeff_2),
            log(2) / (coeff_2 - sigma_b),
            log(2) / (coeff_2 + sigma_b),
            coeff_1,
            sigma_a,
            rho,
        ]

    else

        results_lin_log_fit = [
            label_exp,
            name_well,
            missing,
            missing,
            missing,
            missing,
            missing,
            missing,
            missing,
            missing,
            missing,
            missing,
            missing,
            missing,
        ]

    end


    KinBiont_res_one_well_log_lin = ("Log-lin", results_lin_log_fit, hcat(data_to_fit_times, data_to_fit_values), data_smooted, confidence_band)

    return KinBiont_res_one_well_log_lin
end


"""
    fitting_one_well_ODE_constrained(
    data::Matrix{Float64},
    name_well::String,
    label_exp::String,
    model::String,
    param;
    integrator=Tsit5(),
    pt_avg=1,
    pt_smooth_derivative=7,
    smoothing=false,
    type_of_smoothing="rolling_avg",
    type_of_loss="RE",
    blank_array=zeros(100),
    multiple_scattering_correction=false,
    method_multiple_scattering_correction="interpolation",
    calibration_OD_curve="NA",
    thr_lowess=0.05,
    multistart=false,
    n_restart=50,
    optimizer=BBO_adaptive_de_rand_1_bin_radiuslimited(),
    auto_diff_method=nothing,
    cons=nothing,
    opt_params...
    )

This function uses an ordinary differential equation (ODE) model to fit the data from a single well. It estimates the model parameters while applying specified constraints and optimizations.

# Arguments:

- `data::Matrix{Float64}`: The growth curve data. The first row contains time values, and the second row contains the observable values (e.g., OD).
- `model::String`: The ODE model to be used for fitting.
- `name_well::String`: The name of the well.
- `label_exp::String`: The label for the experiment.
- `param`: Initial guess for the model parameters, provided as a vector of `Float64`.

# Key Arguments:

- `integrator=Tsit5()`: SciML integrator used for solving the ODE. For piecewise models, consider using `KenCarp4(autodiff=true)`.
- `optimizer=BBO_adaptive_de_rand_1_bin_radiuslimited()`: Optimizer for parameter estimation, from the BBO optimization library.
- `type_of_smoothing="rolling_avg"`: Method for smoothing the data. Options include `"NO"`, `"rolling_avg"` (rolling average), and `"lowess"`.
- `pt_avg=7`: Size of the rolling average window for smoothing.
- `smoothing=false`: Boolean flag to apply data smoothing. Set to `true` to smooth the data; `false` to skip smoothing.
- `type_of_loss="RE"`: Type of loss function used for optimization. Options are `"RE"` (relative error), `"L2"` (L2 norm), `"L2_derivative"`, and `"blank_weighted_L2"`.
- `blank_array=zeros(100)`: Array containing data of blanks for correction.
- `pt_smoothing_derivative=7`: Number of points for evaluating the specific growth rate. Uses interpolation if less than 2; otherwise, a sliding window approach is applied.
- `multiple_scattering_correction=false`: Boolean flag to perform multiple scattering correction. Set to `true` to apply correction, requiring a calibration curve.
- `calibration_OD_curve="NA"`: Path to the CSV file containing calibration data, used if `multiple_scattering_correction=true`.
- `method_multiple_scattering_correction="interpolation"`: Method for performing multiple scattering correction. Options are `"interpolation"` or `"exp_fit"` (adapted from Meyers et al., 2018).
- `thr_lowess=0.05`: Threshold parameter for lowess smoothing.
- `auto_diff_method=nothing`: Differentiation method for the optimizer, if required.
- `cons=nothing`: Constraints for optimization.
- `multistart=false`: Flag to enable or disable multistart optimization. Set to `true` to use multiple starting points.
- `n_restart=50`: Number of restarts for multistart optimization, used if `multistart=true`.
- `opt_params...`: Optional parameters for the optimizer (e.g., `lb=[0.1, 0.3], ub=[9.0, 1.0], maxiters=2000000`).

# Output:

- A data structure containing:
  1. `method`: A string describing the method used.
  2. Parameters array: `["name of model", "well", "param_1", "param_2", ..., "param_n", "maximum specific GR using ODE", "maximum specific GR using data", "objective function value (i.e., loss of the solution)"]`, where `"param_1", "param_2", ..., "param_n"` are the ODE model fit parameters.
  3. The numerical solution of the fitted ODE.
  4. Time coordinates corresponding to the fitted ODE.



"""
function fitting_one_well_ODE_constrained(
    data::Matrix{Float64}, # dataset first row times second row OD
    name_well::String, # name of the well
    label_exp::String, #label of the experiment
    model::String, # ode model to use
    param;
    integrator=Tsit5(), # selection of sciml integrator
    pt_avg=1, # numebr of the point to generate intial condition
    pt_smooth_derivative=7,
    smoothing=false, # the smoothing is done or not?
    type_of_smoothing="rolling_avg",
    type_of_loss="RE", # type of used loss
    blank_array=zeros(100), # data of all blanks
    multiple_scattering_correction=false, # if true uses the given calibration curve to fix the data
    method_multiple_scattering_correction="interpolation",
    calibration_OD_curve="NA",  #  the path to calibration curve to fix the data
    thr_lowess=0.05,
    multistart=false,
    n_restart=50,
    optimizer=BBO_adaptive_de_rand_1_bin_radiuslimited(),
    auto_diff_method=nothing,
    cons=nothing,
    opt_params...
)
    if multiple_scattering_correction == true

        data = correction_OD_multiple_scattering(data, calibration_OD_curve; method=method_multiple_scattering_correction)

    end
    #defining time interval

    max_t = data[1, end]
    min_t = data[1, 1]
    tspan = (min_t, max_t)
    tsteps = data[1, :]

    # smoothing data if required


    if smoothing == true
        data = smoothing_data(
            data;
            method=type_of_smoothing,
            pt_avg=pt_avg,
            thr_lowess=thr_lowess
        )
    end

    # setting initial conditions
    u0 = generating_IC(data, model, smoothing, pt_avg)

    #  definition  ode symbolic problem
    ODE_prob = model_selector(model, u0, tspan)

    ## defining loss function
    loss_function = select_loss_function(type_of_loss, data, ODE_prob, integrator, tsteps, blank_array)
    res = KinBiontSolve(loss_function,
        u0,
        param;
        opt=optimizer,
        auto_diff_method=auto_diff_method,
        multistart=multistart,
        n_restart=n_restart,
        cons=cons,
        opt_params...
    )


    #revalution of solution for plot an loss evaluation
    remade_solution = solve(remake(ODE_prob, p=res.u), integrator, saveat=tsteps)
    sol_time = reduce(hcat, remade_solution.t)
    sol_fin = reduce(hcat, remade_solution.u)
    sol_fin = sum(sol_fin, dims=1)

    # max_theoretical gr
    sol_fin, index_not_zero = remove_negative_value(sol_fin)

    data_th = transpose(hcat(sol_time[index_not_zero], sol_fin))


    max_th_gr = maximum(specific_gr_evaluation(Matrix(data_th), pt_smooth_derivative))

    # max empirical gr
    max_em_gr = maximum(specific_gr_evaluation(data, pt_smooth_derivative))

    res_temp = res.u
    loss_value = res.objective
    res_param = vectorize_df_results(label_exp,name_well, model, res_temp, max_th_gr, max_em_gr, loss_value)
    KinBiont_res_one_well = ("ODE", res_param, sol_fin, remade_solution.t)

    return KinBiont_res_one_well
end

#######################################################################
#######################################################################
"""
    fitting_one_well_custom_ODE(
    data::Matrix{Float64},
    name_well::String,
    label_exp::String,
    model::Any,
    param,
    n_equation::Int;
    integrator=Tsit5(),
    pt_avg=1,
    pt_smooth_derivative=0,
    smoothing=false,
    type_of_loss="RE",
    blank_array=zeros(100),
    multiple_scattering_correction=false,
    method_multiple_scattering_correction="interpolation",
    calibration_OD_curve="NA",
    thr_lowess=0.05,
    type_of_smoothing="lowess",
    multistart=false,
    n_restart=50,
    optimizer=BBO_adaptive_de_rand_1_bin_radiuslimited(),
    auto_diff_method=nothing,
    cons=nothing,
    opt_params...
    )

This function fits a user-defined ordinary differential equation (ODE) model to time-series data from a single well. It estimates the model parameters while applying specified constraints and optimizations.

# Arguments:

- `data::Matrix{Float64}`: The growth curve data. The first row contains time values, and the second row contains the observable values (e.g., OD).
- `model::Any`: The user-defined function representing the ODE model to be fitted. The function should define the ODE system to be solved.
- `name_well::String`: The name of the well.
- `label_exp::String`: The label for the experiment.
- `param`: Initial guess for the model parameters, provided as a vector of `Float64`.
- `n_equation::Int`: Number of ODE equations in the model.

# Key Arguments:

- `integrator=Tsit5()`: SciML integrator used for solving the ODE. For piecewise models, consider using `KenCarp4(autodiff=true)`.
- `optimizer=BBO_adaptive_de_rand_1_bin_radiuslimited()`: Optimizer for parameter estimation, from the BBO optimization library.
- `type_of_smoothing="lowess"`: Method for smoothing the data. Options include `"NO"`, `"rolling_avg"` (rolling average), and `"lowess"`.
- `pt_avg=1`: Size of the rolling average window for smoothing.
- `smoothing=false`: Boolean flag to apply data smoothing. Set to `true` to smooth the data; `false` to skip smoothing.
- `type_of_loss="RE"`: Type of loss function used for optimization. Options include `"RE"` (relative error), `"L2"` (L2 norm), `"L2_derivative"`, and `"blank_weighted_L2"`.
- `blank_array=zeros(100)`: Array containing data of blanks for correction.
- `pt_smoothing_derivative=0`: Number of points for evaluating the specific growth rate. If less than 2, uses interpolation; otherwise, uses a sliding window approach.
- `multiple_scattering_correction=false`: Boolean flag to apply multiple scattering correction. Set to `true` to apply correction, requiring a calibration curve.
- `calibration_OD_curve="NA"`: Path to the CSV file containing calibration data, used if `multiple_scattering_correction=true`.
- `method_multiple_scattering_correction="interpolation"`: Method for performing multiple scattering correction. Options are `"interpolation"` or `"exp_fit"` (adapted from Meyers et al., 2018).
- `thr_lowess=0.05`: Threshold parameter for lowess smoothing.
- `auto_diff_method=nothing`: Differentiation method for the optimizer, if required.
- `cons=nothing`: Constraints for optimization.
- `multistart=false`: Boolean flag to enable or disable multistart optimization. Set to `true` to use multiple starting points.
- `n_restart=50`: Number of restarts for multistart optimization, used if `multistart=true`.
- `opt_params...`: Optional parameters for the optimizer (e.g., `lb=[0.1, 0.3], ub=[9.0, 1.0], maxiters=2000000`).

# Output:

- A data structure containing:
  1. `method`: A string describing the method used.
  2. Parameters array: `["name of model", "well", "param_1", "param_2", ..., "param_n", "maximum specific GR using ODE", "maximum specific GR using data", "objective function value (i.e., loss of the solution)"]`, where `"param_1", "param_2", ..., "param_n"` are the ODE model fit parameters.
  3. The numerical solution of the fitted ODE.
  4. Time coordinates corresponding to the fitted ODE.

"""
function fitting_one_well_custom_ODE(
    data::Matrix{Float64}, # dataset first row times second row OD
    name_well::String, # name of the well
    label_exp::String, #label of the experiment
    model::Any, # ode model to use
    param,# initial guess param
    n_equation::Int; # number ode in the system
    integrator=Tsit5(), # selection of sciml integrator
    pt_avg=1, # numebr of the point to generate intial condition
    pt_smooth_derivative=0,
    smoothing=false, # the smoothing is done or not?
    type_of_loss="RE", # type of used loss
    blank_array=zeros(100), # data of all blanks
    multiple_scattering_correction=false, # if true uses the given calibration curve to fix the data
    method_multiple_scattering_correction="interpolation",
    calibration_OD_curve="NA",  #  the path to calibration curve to fix the data
    thr_lowess=0.05,
    type_of_smoothing="lowess",
    multistart=false,
    n_restart=50,
    optimizer=BBO_adaptive_de_rand_1_bin_radiuslimited(),
    auto_diff_method=nothing,
    cons=nothing,
    opt_params...
    )
    if multiple_scattering_correction == true

        data = correction_OD_multiple_scattering(data, calibration_OD_curve; method=method_multiple_scattering_correction)
    end

    #defining time interval
    max_t = data[1, end]
    min_t = data[1, 1]
    tspan = (min_t, max_t)
    tsteps = data[1, :]

    # smoothing data if required
    if smoothing == true
        data = smoothing_data(
            data;
            method=type_of_smoothing,
            pt_avg=pt_avg,
            thr_lowess=thr_lowess
        )

    end

    u0 = generating_IC_custom_ODE(data, n_equation, smoothing, pt_avg)

    #  definition  ode symbolic problem
    ODE_prob = ODEProblem(model, u0, tspan, nothing)

    loss_function =
        select_loss_function(type_of_loss, data, ODE_prob, integrator, tsteps, blank_array)

    res = KinBiontSolve(loss_function,
        u0,
        param;
        opt=optimizer,
        auto_diff_method=auto_diff_method,
        multistart=multistart,
        n_restart=n_restart,
        cons=cons,
        opt_params...
    )

    #revalution of solution for plot an loss evaluation
    remade_solution = solve(remake(ODE_prob, p=res.u), integrator, saveat=tsteps)
    sol_time = reduce(hcat, remade_solution.t)
    sol_fin = reduce(hcat, remade_solution.u)
    sol_fin = sum(sol_fin, dims=1)

    #max_theoretical gr
    sol_fin, index_not_zero = remove_negative_value(sol_fin)

    data_th = transpose(hcat(sol_time[index_not_zero], sol_fin))
    max_th_gr = maximum(specific_gr_evaluation(Matrix(data_th), pt_smooth_derivative))

    # max empirical gr
    max_em_gr = maximum(specific_gr_evaluation(data, pt_smooth_derivative))
    res_temp = res.u
    loss_value = res.objective

    res_param =
        [string(name_well), "custom_model", res_temp, max_th_gr, max_em_gr, loss_value]

    KinBiont_res_one_well = ("custom_ODE", res_param, data_th[2, :], data_th[1, :])


    return KinBiont_res_one_well
end

#######################################################################

"""
    ODE_Model_selection(
    data::Matrix{Float64},
    name_well::String,
    label_exp::String,
    models_list::Vector{String},
    param_array::Any;
    lb_param_array::Any=nothing,
    ub_param_array::Any=nothing,
    integrator=Tsit5(),
    pt_avg=3,
    beta_smoothing_ms=2.0,
    smoothing=false,
    type_of_smoothing="rolling_avg",
    thr_lowess=0.05,
    type_of_loss="L2",
    blank_array=zeros(100),
    pt_smooth_derivative=7,
    multiple_scattering_correction=false,
    method_multiple_scattering_correction="interpolation",
    calibration_OD_curve="NA",
    verbose=false,
    correction_AIC=true,
    multistart=false,
    n_restart=50,
    optimizer=BBO_adaptive_de_rand_1_bin_radiuslimited(),
    auto_diff_method=nothing,
    cons=nothing,
    opt_params...
    )

This function performs automatic model selection for multiple ODE models fitted to time-series data from a single well. It selects the best-fitting model based on the Akaike Information Criterion (AIC) or the corrected AIC (AICc).

# Arguments:

- `data::Matrix{Float64}`: The growth curve data. The first row contains time values, and the second row contains the observable values (e.g., OD).
- `name_well::String`: Name of the well.
- `label_exp::String`: Label of the experiment.
- `models_list::Vector{String}`: A vector of ODE model descriptions to evaluate.
- `param_array::Any`: Initial guess for the model parameters, provided as a vector or matrix.

# Key Arguments:

- `lb_param_array::Any=nothing`: Lower bounds for the parameters, compatible with the models. Use `nothing` for no bounds.
- `ub_param_array::Any=nothing`: Upper bounds for the parameters, compatible with the models. Use `nothing` for no bounds.
- `integrator=Tsit5()`: SciML integrator used for solving the ODEs. Consider `KenCarp4(autodiff=true)` for piecewise models.
- `optimizer=BBO_adaptive_de_rand_1_bin_radiuslimited()`: Optimizer from the BBO optimization library.
- `type_of_smoothing="rolling_avg"`: Method for smoothing the data. Options: `"NO"`, `"rolling_avg"` (rolling average), and `"lowess"`.
- `pt_avg=3`: Size of the rolling average window for smoothing.
- `beta_smoothing_ms=2.0`: Penalty parameter for evaluating AIC (or AICc).
- `smoothing=false`: Boolean flag to apply data smoothing. Set to `true` to smooth the data; `false` to skip smoothing.
- `type_of_loss="L2"`: Type of loss function used for optimization. Options include `"L2"` (L2 norm), `"RE"` (relative error), `"L2_derivative"`, and `"blank_weighted_L2"`.
- `blank_array=zeros(100)`: Array containing data of blanks for correction.
- `pt_smooth_derivative=7`: Number of points for evaluating the specific growth rate. Uses interpolation if less than 2; otherwise, uses a sliding window approach.
- `multiple_scattering_correction=false`: Boolean flag to apply multiple scattering correction. Set to `true` if correction is required, necessitating a calibration curve.
- `calibration_OD_curve="NA"`: Path to the CSV file with calibration data, used if `multiple_scattering_correction=true`.
- `method_multiple_scattering_correction="interpolation"`: Method for multiple scattering correction. Options are `"interpolation"` or `"exp_fit"` (adapted from Meyers et al., 2018).
- `thr_lowess=0.05`: Parameter for lowess smoothing.
- `correction_AIC=true`: Boolean flag to apply finite sample correction to AIC (AICc). Set to `true` to correct.
- `auto_diff_method=nothing`: Differentiation method for the optimizer, if required.
- `cons=nothing`: Constraints for optimization.
- `multistart=false`: Boolean flag to enable or disable multistart optimization. Set to `true` to use multiple starting points.
- `n_restart=50`: Number of restarts for multistart optimization, used if `multistart=true`.
- `opt_params...`: Optional parameters for the optimizer (e.g., `lb=[0.1, 0.3], ub=[9.0, 1.0], maxiters=2000000`).

# Output:

- A data structure containing:
  1. `method`: A string describing the method used for model selection.
  2. `best_model_params`: Matrix containing the parameters of the best model.
  3. `fit_numerical`: Numerical array of the best model fit.
  4. `fit_times`: Time coordinates corresponding to the best model fit.
  5. `model_stats`: Statistical summary for each model fitted, including AIC/AICc values.
  6. `best_model_score`: Score (AIC/AICc) of the best model.
  7. `best_model_params_values`: Parameters of the best model.
  8. `min_aic_or_aicc`: Minimum value of AIC or AICc for the best model.
  9. `best_model_string`: String representation of the best model.
  10. `all_params`: Parameters for all evaluated models.

"""
function ODE_Model_selection(
    data::Matrix{Float64}, # dataset first row times second row OD
    name_well::String, # name of the well
    label_exp::String, #label of the experiment
    models_list::Vector{String}, # ode model to use
    param_array::Any;
    lb_param_array::Any=nothing, # lower bound param
    ub_param_array::Any=nothing, # upper bound param
    integrator=Tsit5(), # selection of sciml integrator
    pt_avg=3, # number of the point to generate intial condition
    beta_smoothing_ms=2.0, # penality for AIC evaluation
    smoothing=false, # the smoothing is done or not?
    type_of_smoothing="rolling_avg",
    thr_lowess=0.05,
    type_of_loss="L2", # type of used loss
    blank_array=zeros(100), # data of all blanks
    pt_smooth_derivative=7,
    multiple_scattering_correction=false, # if true uses the given calibration curve to fix the data
    method_multiple_scattering_correction="interpolation",
    calibration_OD_curve="NA", #  the path to calibration curve to fix the data
    verbose=false,
    correction_AIC=true,
    multistart=false,
    n_restart=50,
    optimizer=BBO_adaptive_de_rand_1_bin_radiuslimited(),
    auto_diff_method=nothing,
    cons=nothing,
    opt_params...
)
    opt_param_temp = copy(opt_params)

    if multiple_scattering_correction == true

        data = correction_OD_multiple_scattering(data, calibration_OD_curve; method=method_multiple_scattering_correction)
    end

    # smooting if required
    if smoothing == true
        data = smoothing_data(
            data;
            method=type_of_smoothing,
            pt_avg=pt_avg,
            thr_lowess=thr_lowess
        )
    end

    # inizialization of array of results
    df_res_optimization = Array{Any}(nothing, length(models_list))
    rss_array = ["model", "Loss", "AIC_standard"]

    # generating time interval
    max_t = data[1, end]
    min_t = data[1, 1]
    tspan = (min_t, max_t)
    tsteps = data[1, :]
    data = Matrix(data)

    # common part loss definition
    ## defining loss function
    # loop over the models

    for mm in eachindex(models_list)
        if verbose == true
            println(string("fitting ", models_list[mm]))
        end
        # inizialization of some parameters

        temp_model = models_list[mm]

        temp_start_param = param_array[mm]

   

        # generating IC
        # setting initial conditions
        u0 = generating_IC(data, temp_model, smoothing, pt_avg)

        #  definition  ode symbolic problem
        ODE_prob = model_selector(temp_model, u0, tspan)

        loss_function = select_loss_function(
            type_of_loss,
            data,
            ODE_prob,
            integrator,
            tsteps,
            blank_array,
        )


        if !isnothing(lb_param_array)
            temp_param_lb = lb_param_array[mm]
            temp_param_ub = ub_param_array[mm]

            opt_param_temp = (opt_params...,
                lb=temp_param_lb,
                ub=temp_param_ub,
                )

            res = KinBiontSolve(loss_function,
                u0,
                temp_start_param;
                opt=optimizer,
                auto_diff_method=auto_diff_method,
                multistart=multistart,
                n_restart=n_restart,
                cons=cons,
                opt_param_temp...)

        else

            res = KinBiontSolve(loss_function,
                u0,
                temp_start_param;
                opt=optimizer,
                auto_diff_method=auto_diff_method,
                multistart=multistart,
                n_restart=n_restart,
                cons=cons,
                opt_params...)

        end



        #revalution of solution for plot an loss evaluation
        remade_solution = solve(remake(ODE_prob, p=res.u), integrator, saveat=tsteps)
        sol_time = reduce(hcat, remade_solution.t)
        sol_fin = reduce(hcat, remade_solution.u)
        sol_fin = sum(sol_fin, dims=1)

        param_number = length(temp_start_param)

        #  AICc = AICc_evaluation(param_number, beta_smoothing_ms, data[2, :], sol_fin, correction=correction_AIC)
        AICc = AICc_evaluation2(param_number, beta_smoothing_ms, data[2, :], res.objective, correction=correction_AIC)

        #max_theoretical gr
        sol_fin, index_not_zero = remove_negative_value(sol_fin)


        data_th = transpose(hcat(sol_time[index_not_zero], sol_fin))
        max_th_gr = maximum(specific_gr_evaluation(Matrix(data_th), pt_smooth_derivative))

        # max empirical gr
        max_em_gr = maximum(specific_gr_evaluation(data, pt_smooth_derivative))


        res_param = vectorize_df_results(
            label_exp,
            name_well,
            temp_model,
            res.u,
            max_th_gr,
            max_em_gr,
            res.objective,
        )


        results_to_be_pushed = [
            models_list[mm],
            res_param[end],
            AICc,
        ]
        rss_array = hcat(rss_array, results_to_be_pushed)
        df_res_optimization[mm] = res_param

    end
    AIC_array = rss_array[3, 2:end]
    min_AIC = minimum(AIC_array)
    index_minimal_AIC_model = findfirst(item -> item == min_AIC, AIC_array) + 1

    # string of the model choosen
    model = rss_array[1, index_minimal_AIC_model]

    # param of the best model
    param_min = df_res_optimization[index_minimal_AIC_model-1]
    param_out_full = copy(param_min)
    param_min = param_min[4:(end-3)]

    tsteps = data[1, :]
    tspan = (data[1, 1], data[1, end])
    u0 = generating_IC(data, model, smoothing, pt_avg)
    ODE_prob = model_selector(model, u0, tspan, param_min)
    sim = solve(ODE_prob, integrator, saveat=tsteps)
    sol_t = reduce(hcat, sim.u)
    sol_time = reduce(hcat, sim.t)
    sol_t = sum(sol_t, dims=1)
    sol_fin, index_not_zero = remove_negative_value(sol_t)

    data_th = transpose(hcat(sol_time[index_not_zero], sol_fin))

    sol_fin, index_not_zero = remove_negative_value(sol_fin)

    data_th = transpose(hcat(sol_time[index_not_zero], sol_fin))

    KinBiont_res_model_selection = ("ODE_model_selection",
        df_res_optimization,
        sol_fin,
        sol_time[index_not_zero],
        rss_array,
        minimum(rss_array[2, 2:end]),
        param_min,
        min_AIC,
        model,
        param_out_full
    )



    return KinBiont_res_model_selection

end



#######################################################################

"""
    one_well_morris_sensitivity(
    data::Matrix{Float64},
    name_well::String,
    label_exp::String,
    model::String,
    lb_param::Vector{Float64},
    ub_param::Vector{Float64};
    N_step_morris=7,
    integrator=Tsit5(),
    pt_avg=1,
    pt_smooth_derivative=7,
    write_res=false,
    smoothing=false,
    type_of_smoothing="rolling_avg",
    type_of_loss="RE",
    blank_array=zeros(100),
    multiple_scattering_correction=false,
    method_multiple_scattering_correction="interpolation",
    calibration_OD_curve="NA",
    optimizer=BBO_adaptive_de_rand_1_bin_radiuslimited(),
    auto_diff_method=nothing,
    cons=nothing,
    thr_lowess=0.05,
    multistart=false,
    n_restart=50,
    opt_params...
    )

This function performs Morris sensitivity analysis to evaluate the sensitivity of model parameters to variations in their initial guesses. This analysis is useful for assessing the robustness of nonlinear model fits.

# Arguments:

- `data::Matrix{Float64}`: The growth curve data. Time values are in the first row, and the fit observable (e.g., OD) is in the second row.
- `name_well::String`: Name of the well.
- `label_exp::String`: Label of the experiment.
- `model::String`: The ODE model to be used for fitting.
- `lb_param::Vector{Float64}`: Vector containing the lower bounds for the model parameters.
- `ub_param::Vector{Float64}`: Vector containing the upper bounds for the model parameters.

# Key Arguments:

- `N_step_morris=7`: Number of steps for the Morris sensitivity analysis.
- `param=lb_param .+ (ub_param.-lb_param)./2`: Initial guess for the model parameters, calculated as the midpoint of the lower and upper bounds.
- `integrator=Tsit5()`: SciML integrator used for solving the ODEs. Consider `KenCarp4(autodiff=true)` for piecewise models.
- `optimizer=BBO_adaptive_de_rand_1_bin_radiuslimited()`: Optimizer used for fitting the model.
- `type_of_smoothing="rolling_avg"`: Method for smoothing the data. Options include `"NO"`, `"rolling_avg"` (rolling average), and `"lowess"`.
- `pt_avg=7`: Size of the rolling average window for smoothing.
- `pt_smooth_derivative=7`: Number of points used for evaluating the specific growth rate. If less than 2, interpolation is used; otherwise, a sliding window approach is used.
- `smoothing=false`: Boolean flag indicating whether to apply smoothing to the data (`true`) or not (`false`).
- `type_of_loss="RE"`: Type of loss function used for optimization. Options include `"RE"` (relative error), `"L2"` (L2 norm), `"L2_derivative"`, and `"blank_weighted_L2"`.
- `blank_array=zeros(100)`: Array containing data of blanks for correction.
- `calibration_OD_curve="NA"`: Path to the CSV file with calibration data, used if `multiple_scattering_correction=true`.
- `multiple_scattering_correction=false`: Boolean flag to apply multiple scattering correction. Set to `true` if correction is required.
- `method_multiple_scattering_correction="interpolation"`: Method for multiple scattering correction. Options include `"interpolation"` or `"exp_fit"` (adapted from Meyers et al., 2018).
- `thr_lowess=0.05`: Parameter for lowess smoothing.
- `auto_diff_method=nothing`: Differentiation method for the optimizer, if required.
- `cons=nothing`: Constraints for optimization.
- `multistart=false`: Boolean flag to enable or disable multistart optimization. Set to `true` to use multiple starting points.
- `n_restart=50`: Number of restarts for multistart optimization, used if `multistart=true`.
- `opt_params...`: Optional parameters for the optimizer (e.g., `lb=[0.1, 0.3], ub=[9.0, 1.0], maxiters=2000000`).

# Output:

- A data structure containing:
  1. `method`: A string describing the Morris sensitivity analysis method used.
  2. `results_fit`: The results of the model fits for all starting parameter sets used in the sensitivity analysis.
  3. `initial_guesses`: Initial guess of parameters for each run of the sensitivity analysis.

"""
function one_well_morris_sensitivity(
    data::Matrix{Float64}, # dataset first row times second row OD
    name_well::String, # name of the well
    label_exp::String, #label of the experiment
    model::String, # ode model to use
    lb_param::Vector{Float64}, # lower bound param
    ub_param::Vector{Float64}; # upper bound param
    N_step_morris=7,
    integrator=Tsit5(), # selection of sciml integrator
    pt_avg=1, # numebr of the point to generate intial condition
    pt_smooth_derivative=7,
    write_res=false,
    smoothing=false, # the smoothing is done or not?
    type_of_smoothing="rolling_avg",
    type_of_loss="RE", # type of used loss
    blank_array=zeros(100), # data of all blanks
    multiple_scattering_correction=false, # if true uses the given calibration curve to fix the data
    method_multiple_scattering_correction="interpolation",
    calibration_OD_curve="NA",  #  the path to calibration curve to fix the data
    optimizer=BBO_adaptive_de_rand_1_bin_radiuslimited(),
    auto_diff_method=nothing,
    cons=nothing,
    thr_lowess = 0.05,
    multistart=false,
    n_restart=50,
    opt_params...
)

    # inizializing the results of sensitivity
    results_sensitivity = initialize_df_results(model)

    if write_res == true
        mkpath(path_to_results)
    end

    if multiple_scattering_correction == true

        data = correction_OD_multiple_scattering(data, calibration_OD_curve; method=method_multiple_scattering_correction)

    end
    #defining time interval
    max_t = data[1, end]
    min_t = data[1, 1]
    tspan = (min_t, max_t)
    tsteps = data[1, :]

    # smoothing data if required
    if smoothing == true
        data = smoothing_data(
            data;
            method=type_of_smoothing,
            pt_avg=pt_avg,
            thr_lowess=thr_lowess
        )
    end

    # setting initial conditions
    u0 = generating_IC(data, model, smoothing, pt_avg)

    #  definition  ode symbolic problem
    ODE_prob = model_selector(model, u0, tspan)
    loss_function =
        select_loss_function(type_of_loss, data, ODE_prob, integrator, tsteps, blank_array)

    param_combination =
        generation_of_combination_of_IC_morris(lb_param, ub_param, N_step_morris)

    for i = 1:size(param_combination)[2]

        param = param_combination[:, i]
        opt_params = (opt_params...,
        lb=lb_param,
        ub=ub_param,
        )
        res = KinBiontSolve(loss_function,
            u0,
            param;
            opt=optimizer,
            auto_diff_method=auto_diff_method,
            multistart=multistart,
            n_restart=n_restart,
            cons=cons,
            opt_params...)

        #revalution of solution for plot an loss evaluation
        remade_solution = solve(remake(ODE_prob, p=res.u), integrator, saveat=tsteps)
        sol_time = reduce(hcat, remade_solution.t)
        sol_fin = reduce(hcat, remade_solution.u)
        sol_fin = sum(sol_fin, dims=1)
        loss_value = res.objective

        #max_theoretical gr
        sol_fin, index_not_zero = remove_negative_value(sol_fin)

        data_th = transpose(hcat(sol_time[index_not_zero], sol_fin))
        max_th_gr = maximum(specific_gr_evaluation(Matrix(data_th), pt_smooth_derivative))

        # max empirical gr
        max_em_gr = maximum(specific_gr_evaluation(data, pt_smooth_derivative))
        res_temp = res.u
        res_param = vectorize_df_results(
            label_exp,
            name_well,
            model,
            res_temp,
            max_th_gr,
            max_em_gr,
            loss_value,
        )
        results_sensitivity = hcat(results_sensitivity, res_param)
    end

    if write_res == true
        mkpath(path_to_results)
        CSV.write(
            string(path_to_results, label_exp, "_results_sensitivity.csv"),
            Tables.table(Matrix(results_sensitivity)),
        )
        CSV.write(
            string(path_to_results, label_exp, "_configuration_tested.csv"),
            Tables.table(Matrix(param_combination)),
        )
    end

    KinBiont_res_sensitivity = ("ODE_Morris_sensitivity", results_sensitivity, param_combination)

    return KinBiont_res_sensitivity
end



"""
    selection_ODE_fixed_intervals(
    data::Matrix{Float64},
    name_well::String,
    label_exp::String,
    list_of_models::Vector{String},
    param_array,
    intervals_changepoints::Any;
    lb_param_array::Any=nothing,
    ub_param_array::Any=nothing,
    type_of_loss="L2",
    integrator=Tsit5(),
    smoothing=false,
    type_of_smoothing="lowess",
    thr_lowess=0.05,
    pt_avg=1,
    pt_smooth_derivative=0,
    multiple_scattering_correction=false,
    method_multiple_scattering_correction="interpolation",
    calibration_OD_curve="NA",
    beta_smoothing_ms=2.0,
    correction_AIC=true,
    optimizer=BBO_adaptive_de_rand_1_bin_radiuslimited(),
    multistart=false,
    n_restart=50,
    auto_diff_method=nothing,
    cons=nothing,
    opt_params...
    )

This function fits an Ordinary Differential Equation (ODE) model to segmented time-series data. Users provide fixed change points, and the function fits the models to each segment defined by these points. 

# Arguments:

- `data::Matrix{Float64}`: The growth curve data. Time values are in the first row, and the fit observable (e.g., OD) is in the second row.
- `name_well::String`: Name of the well.
- `label_exp::String`: Label of the experiment.
- `list_of_models::Vector{String}`: List of ODE models to be considered for fitting.
- `param_array`: Vector of initial guesses for model parameters.
- `intervals_changepoints::Any`: Array containing the list of change points, e.g., `[0.0, 10.0, 30.0]`. These define the segments for which models will be fitted.

# Key Arguments:

- `lb_param_array::Any`: Lower bounds for the parameters for each model (compatible with the models).
- `ub_param_array::Any`: Upper bounds for the parameters for each model (compatible with the models).
- `integrator=Tsit5()`: SciML integrator used for solving the ODEs. Use `KenCarp4(autodiff=true)` for piecewise models.
- `optimizer=BBO_adaptive_de_rand_1_bin_radiuslimited()`: Optimizer from optimizationBBO.
- `type_of_smoothing="lowess"`: Method for smoothing the data. Options include `"NO"`, `"rolling_avg"` (rolling average), and `"lowess"`.
- `pt_avg=1`: Size of the rolling average window for smoothing.
- `smoothing=false`: Boolean flag indicating whether to apply smoothing to the data (`true`) or not (`false`).
- `type_of_loss="L2"`: Type of loss function used for optimization. Options include `"RE"` (relative error), `"L2"` (L2 norm), `"L2_derivative"`, and `"blank_weighted_L2"`.
- `blank_array=zeros(100)`: Array containing data of blanks for correction.
- `pt_smooth_derivative=0`: Number of points for evaluation of specific growth rate. If less than 2, interpolation is used; otherwise, a sliding window approach is used.
- `calibration_OD_curve="NA"`: Path to the CSV file with calibration data, used if `multiple_scattering_correction=true`.
- `multiple_scattering_correction=false`: Boolean flag to apply multiple scattering correction. Set to `true` if correction is needed.
- `method_multiple_scattering_correction="interpolation"`: Method for multiple scattering correction. Options include `"interpolation"` or `"exp_fit"` (adapted from Meyers et al., 2018).
- `thr_lowess=0.05`: Parameter for lowess smoothing.
- `beta_smoothing_ms=2.0`: Penalty parameter for AIC (or AICc) evaluation.
- `auto_diff_method=nothing`: Differentiation method for the optimizer, if required.
- `cons=nothing`: Constraints for optimization.
- `multistart=false`: Boolean flag to enable or disable multistart optimization. Set to `true` to use multiple starting points.
- `n_restart=50`: Number of restarts for multistart optimization, used if `multistart=true`.
- `opt_params...`: Optional parameters for the optimizer (e.g., `lb=[0.1, 0.3], ub=[9.0, 1.0], maxiters=2000000`).

# Output:

If `res = selection_ODE_fixed_intervals(...)`:

- `res[1]`: Parameters of each segment.
- `res[2]`: Intervals corresponding to each ODE segment.
- `res[3]`: Time coordinates of the fitted solution.
- `res[4]`: Numerical values of the fitted solutions.
- `res[5]`: The fit loss score for each segment.


"""
function selection_ODE_fixed_intervals(
    data_testing::Matrix{Float64}, # dataset first row times second row OD
    name_well::String, # name of the well
    label_exp::String, #label of the experiment
    list_of_models::Vector{String}, # ode models to use
    param_array,
    intervals_changepoints::Any;
    lb_param_array::Any=nothing, # lower bound param
    ub_param_array::Any=nothing, # upper bound param
    type_of_loss="L2", # type of used loss
    integrator=Tsit5(), # selection of sciml integrator
    smoothing=false,
    type_of_smoothing="lowess",
    thr_lowess=0.05,
    pt_avg=1,
    pt_smooth_derivative=0,
    multiple_scattering_correction=false, # if true uses the given calibration curve to fix the data
    method_multiple_scattering_correction="interpolation",
    calibration_OD_curve="NA", #  the path to calibration curve to fix the data
    beta_smoothing_ms=2.0, #  parameter of the AIC penality
    correction_AIC=true,
    optimizer=BBO_adaptive_de_rand_1_bin_radiuslimited(),
    multistart=false,
    n_restart=50,
    auto_diff_method=nothing,
    cons=nothing,
    opt_params...)

    if multiple_scattering_correction == true
        data_testing = correction_OD_multiple_scattering(data_testing, calibration_OD_curve; method=method_multiple_scattering_correction)
    end
    total_loss = 0.0
    if smoothing == true
        data_testing = smoothing_data(
            data_testing;
            method=type_of_smoothing,
            pt_avg=pt_avg,
            thr_lowess=thr_lowess
        )
    end

    interval_changepoints = copy(intervals_changepoints)
    interval_changepoints = push!(interval_changepoints, data_testing[1, 1])
    interval_changepoints = push!(interval_changepoints, data_testing[1, end])
    interval_changepoints = sort(interval_changepoints)
    param_out = Vector{Vector{Any}}()
    composed_sol = Type{Any}
    composed_time = Type{Any}
    bc = [data_testing[1, 1], data_testing[2, 2]]

    for i = 2:(eachindex(interval_changepoints)[end])
        if i == 2
            tspan_array = findall((data_testing[1, :] .<= interval_changepoints[i]))
            data_temp = Matrix(
                transpose(hcat(data_testing[1, tspan_array], data_testing[2, tspan_array])),
            )
        else
            tspan_array_1 = findall((data_testing[1, :] .> interval_changepoints[i-1]))
            tspan_array_2 = findall((data_testing[1, :] .<= interval_changepoints[i]))
            tspan_array = intersect(tspan_array_1, tspan_array_2)
            data_temp = Matrix(
                transpose(hcat(data_testing[1, tspan_array], data_testing[2, tspan_array])),
            )
            # imposing_bounduary condition
            data_temp = hcat(bc, data_temp)
        end
        model_selection_results = ODE_Model_selection(
            data_temp, # dataset first row times second row OD
            name_well, # name of the well
            label_exp, #label of the experiment
            list_of_models, # ode model to use
            param_array;
            ub_param_array=ub_param_array, # lower bound param
            lb_param_array=lb_param_array, # upper bound param
            optimizer=optimizer, # selection of optimization method
            integrator=integrator, # selection of sciml integrator
            pt_avg=pt_avg, # number of the point to generate intial condition
            beta_smoothing_ms=beta_smoothing_ms, # penality for AIC evaluation
            smoothing=smoothing, # the smoothing is done or not?
            type_of_loss=type_of_loss, # type of used loss
            blank_array=zeros(100), # data of all blanks
            pt_smooth_derivative=pt_smooth_derivative,
            multiple_scattering_correction=false, # if true uses the given calibration curve to fix the data
            method_multiple_scattering_correction="interpolation",
            calibration_OD_curve="NA", #  the path to calibration curve to fix the data
            verbose=false,
            correction_AIC=correction_AIC,
            multistart=multistart,
            n_restart=n_restart,
            auto_diff_method=auto_diff_method,
            cons=cons,
            opt_params...
        )

        # selection of te model
        model = model_selection_results[9]
        total_loss = total_loss + model_selection_results[end][end]
        # param of the best model
        temp_res_win = model_selection_results[7]
        param_fitting = copy(temp_res_win)
        #temp_res_win = push!(temp_res_win,model)
        u0 = generating_IC(data_temp, model, smoothing, pt_avg)


        time_sol = model_selection_results[4]
        sol_fin = model_selection_results[3]


        time_bonduary = time_sol[end]
        value_bonduary = sol_fin[end]
        bc = [time_bonduary, value_bonduary]

        sol_fin, index_not_zero = remove_negative_value(sol_fin)

        if (
            pt_smooth_derivative > (length(data_temp[1, :]) - 2) &&
            length(data_temp[1, :]) > 3
        )
            emp_max_gr_of_segment =
                maximum(specific_gr_evaluation(data_temp, pt_smooth_derivative))
            data_th = transpose(hcat(time_sol[index_not_zero], sol_fin))
            th_max_gr_of_segment =
                maximum(specific_gr_evaluation(data_th, pt_smooth_derivative))
            temp_res_win = vcat(temp_res_win, th_max_gr_of_segment)
            temp_res_win = vcat(temp_res_win, emp_max_gr_of_segment)

        elseif length(data_temp[1, :]) <= 3

            emp_max_gr_of_segment = missing
            th_max_gr_of_segment = missing
            temp_res_win = vcat(temp_res_win, th_max_gr_of_segment)
            temp_res_win = vcat(temp_res_win, emp_max_gr_of_segment)
        else
            emp_max_gr_of_segment =
                maximum(specific_gr_evaluation(data_temp, pt_smooth_derivative))
            data_th = transpose(hcat(time_sol[index_not_zero], sol_fin))
            th_max_gr_of_segment =
                maximum(specific_gr_evaluation(data_th, pt_smooth_derivative))
            temp_res_win = vcat(temp_res_win, th_max_gr_of_segment)
            temp_res_win = vcat(temp_res_win, emp_max_gr_of_segment)
        end
        temp_res_win = vcat(model, temp_res_win)
        temp_res_win = vcat(temp_res_win, model_selection_results[6])
        temp_res_win = vcat(name_well, temp_res_win)
        temp_res_win = vcat(label_exp, temp_res_win)







        if i == 2
            composed_time = copy(time_sol[index_not_zero])
            composed_sol = copy(sol_fin)
            temp_res_win = push!(temp_res_win, i - 1)
            param_out = push!(param_out, temp_res_win)
        else


            composed_time = vcat(composed_time, time_sol[index_not_zero])
            composed_sol = vcat(composed_sol, sol_fin)
            temp_res_win = push!(temp_res_win, i - 1)
            param_out = push!(param_out, temp_res_win)
        end
    end

    composed_time, composed_sol = remove_replicate_data(composed_time, composed_sol)

    return param_out, interval_changepoints, composed_time, composed_sol, total_loss
end



"""
    segmentation_ODE(
    data_testing::Matrix{Float64},
    name_well::String,
    label_exp::String,
    list_of_models::Vector{String},
    param_array::Any,
    n_max_change_points::Int;
    lb_param_array::Any=nothing,
    ub_param_array::Any=nothing,
    detect_number_cpd=true,
    fixed_cpd=false,
    integrator=Tsit5(),
    type_of_loss="L2",
    type_of_detection="slinding_win",
    type_of_curve="original",
    pt_avg=1,
    smoothing=true,
    path_to_results="NA",
    win_size=14,
    pt_smooth_derivative=7,
    beta_smoothing_ms=2.0,
    multiple_scattering_correction=false,
    method_multiple_scattering_correction="interpolation",
    calibration_OD_curve="NA",
    save_all_model=false,
    method_peaks_detection="peaks_prominence",
    n_bins=40,
    type_of_smoothing="lowess",
    thr_lowess=0.05,
    correction_AIC=true,
    optimizer=BBO_adaptive_de_rand_1_bin_radiuslimited(),
    multistart=false,
    n_restart=50,
    auto_diff_method=nothing,
    cons=nothing,
    opt_params...
    )

This function performs model selection for ordinary differential equation (ODE) models across different segments of the input growth time series data. Segmentation is achieved using a change point detection algorithm, allowing the identification of multiple segments where different models can be fitted.

# Arguments:

- `data_testing::Matrix{Float64}`: The growth curve data. Time values are in the first row, and the fit observable (e.g., OD) is in the second row.
- `name_well::String`: Name of the well.
- `label_exp::String`: Label of the experiment.
- `list_of_models::Vector{String}`: List of ODE models to be considered for fitting.
- `param_array::Any`: Initial guesses for the model parameters.
- `n_max_change_points::Int`: Maximum number of change points to be considered. The function will evaluate models with varying numbers of change points up to this maximum.

# Key Arguments:

- `lb_param_array::Any`: Lower bounds for the parameters for each model (compatible with the models).
- `ub_param_array::Any`: Upper bounds for the parameters for each model (compatible with the models).
- `integrator=Tsit5()`: SciML integrator used for solving the ODEs. Use `KenCarp4(autodiff=true)` for piecewise models.
- `optimizer=BBO_adaptive_de_rand_1_bin_radiuslimited()`: Optimizer from optimizationBBO.
- `type_of_smoothing="lowess"`: Method of choice for smoothing the data. Options include `"NO"`, `"rolling_avg"` (rolling average), and `"lowess"`.
- `pt_avg=1`: Size of the rolling average window for smoothing.
- `pt_smooth_derivative=7`: Number of points for evaluation of the specific growth rate. If less than 2, interpolation is used; otherwise, a sliding window approach is used.
- `smoothing=true`: Boolean flag indicating whether to apply smoothing to the data (`true`) or not (`false`).
- `type_of_loss="L2"`: Type of loss function used for optimization. Options include `"RE"` (relative error), `"L2"` (L2 norm), `"L2_derivative"`, and `"blank_weighted_L2"`.
- `blank_array=zeros(100)`: Array containing data of blanks for correction.
- `calibration_OD_curve="NA"`: Path to the CSV file with calibration data, used if `multiple_scattering_correction=true`.
- `multiple_scattering_correction=false`: Boolean flag to apply multiple scattering correction. Set to `true` if correction is needed.
- `method_multiple_scattering_correction="interpolation"`: Method for multiple scattering correction. Options include `"interpolation"` or `"exp_fit"` (adapted from Meyers et al., 2018).
- `thr_lowess=0.05`: Parameter for lowess smoothing.
- `beta_smoothing_ms=2.0`: Penalty parameter for AIC (or AICc) evaluation.
- `type_of_detection="slinding_win"`: Method of change point detection. Options include `"slinding_win"` (sliding window approach) and `"lsdd"` (least square density difference).
- `type_of_curve="original"`: Defines the input curve for change point detection. Options include `"original"` for the raw time series and `"deriv"` for the specific growth rate time series.
- `method_peaks_detection="peaks_prominence"`: Method for peak detection in the dissimilarity curve. Options include `"peaks_prominence"` (orders peaks by prominence) and `"thr_scan"` (uses a threshold to identify peaks).
- `n_bins=40`: Number of bins used for threshold scanning if `method_peaks_detection="thr_scan"`.
- `detect_number_cpd=true`: Boolean flag to test all possible combinations of change points up to `n_max_change_points`. Set to `false` to only fit using `n_max_change_points`.
- `fixed_cpd=false`: Boolean flag to return the fit using exactly `n_max_change_points` change points if `true`.
- `win_size=14`: Size of the window used by change point detection algorithms.
- `path_to_results="NA"`: Path to save the results.
- `save_all_model=false`: Boolean flag to save all tested models if `true`.
- `auto_diff_method=nothing`: Differentiation method for the optimizer, if required.
- `cons=nothing`: Constraints for optimization.
- `multistart=false`: Boolean flag to enable or disable multistart optimization. Set to `true` to use multiple starting points.
- `n_restart=50`: Number of restarts for multistart optimization, used if `multistart=true`.
- `opt_params...`: Optional parameters for the optimizer (e.g., `lb=[0.1, 0.3], ub=[9.0, 1.0], maxiters=2000000`).

# Output:

If `res = segmentation_ODE(...)`:

- `res[1]`: String describing the method used for segmentation and model fitting.
- `res[2]`: Matrix of parameters for each ODE segment.
- `res[3]`: Numerical values of the fitted solutions.
- `res[4]`: Time coordinates of the fitted solutions.
- `res[5]`: Intervals corresponding to the identified change points.
- `res[6]`: AICc (or AIC) of the final model, indicating the goodness of fit.

"""
function segmentation_ODE(
    data_testing::Matrix{Float64}, # dataset x times y OD/fluorescence
    name_well::String, # name of the well
    label_exp::String, #label of the experiment
    list_of_models::Vector{String}, # ode model to use
    param_array::Any, #  param
    n_max_change_points::Int;
    lb_param_array::Any=nothing, # lower bound param
    ub_param_array::Any=nothing, # upper bound param
    detect_number_cpd=true,
    fixed_cpd=false,
    integrator=Tsit5(), # selection of sciml integrator
    type_of_loss="L2", # type of used loss
    type_of_detection="slinding_win",
    type_of_curve="original",
    pt_avg=1, # number of the point to generate intial condition
    smoothing=true, # the smoothing is done or not?
    path_to_results="NA",
    win_size=14, #  
    pt_smooth_derivative=7,
    beta_smoothing_ms=2.0,
    multiple_scattering_correction=false, # if true uses the given calibration curve to fix the data
    method_multiple_scattering_correction="interpolation",
    calibration_OD_curve="NA",  #  the path to calibration curve to fix the data
    save_all_model=false,
    method_peaks_detection="peaks_prominence",
    n_bins=40,
    type_of_smoothing="lowess",
    thr_lowess=0.05,
    correction_AIC=true,
    optimizer=BBO_adaptive_de_rand_1_bin_radiuslimited(),
    multistart=false,
    n_restart=50,
    auto_diff_method=nothing,
    cons=nothing,
    opt_params...
)

    # fitting single models
    change_point_list = Vector{Vector{Any}}()


    if n_max_change_points == 0 || detect_number_cpd == true
        res = ODE_Model_selection(
            data_testing, # dataset first row times second row OD
            name_well, # name of the well
            label_exp, #label of the experiment
            list_of_models, # ode model to use
            param_array;
            lb_param_array=lb_param_array, # lower bound param
            ub_param_array=ub_param_array, # upper bound param
            optimizer=optimizer, # selection of optimization method
            integrator=integrator, # selection of sciml integrator
            pt_avg=pt_avg, # number of the point to generate intial condition
            beta_smoothing_ms=beta_smoothing_ms, # penality for AIC evaluation
            smoothing=smoothing, # the smoothing is done or not?
            type_of_loss=type_of_loss, # type of used loss
            blank_array=zeros(100), # data of all blanks
            pt_smooth_derivative=pt_smooth_derivative,
            multiple_scattering_correction=multiple_scattering_correction, # if true uses the given calibration curve to fix the data
            method_multiple_scattering_correction=method_multiple_scattering_correction,
            calibration_OD_curve=calibration_OD_curve, #  the path to calibration curve to fix the data
            verbose=false,
            type_of_smoothing=type_of_smoothing,
            thr_lowess=thr_lowess,
            correction_AIC=correction_AIC,
            auto_diff_method=auto_diff_method,
            cons=cons,
            opt_params...)

        if save_all_model == true
            mkpath(path_to_results)
            CSV.write(
                string(
                    path_to_results,
                    label_exp,
                    "_segmented_ODE_scoring_",
                    name_well,
                    "_seg_0.csv",
                ),
                Tables.table(res[8]),
            )
            CSV.write(
                string(
                    path_to_results,
                    label_exp,
                    "_segmented_ODE_param_",
                    name_well,
                    "_seg_0.csv",
                ),
                Tables.table(Vector(res[end])),
            )

        end
        top_cps = [0.0]
        top_model = res[2]
        score_of_the_models = res[8]
        change_point_list = [0.0]
        change_point_to_plot = [0.0, 0.0, 0.0]
        time_points_to_plot = res[4]
        sol_to_plot = res[3]

    else

        top_model =[ "NO"]
        score_of_the_models = 10^9
        change_point_list = [0.0]
        change_point_to_plot = [0.0, 0.0, 0.0]
        time_points_to_plot = copy(data_testing[1, :])
        sol_to_plot = copy(data_testing[1, :])
        top_cps = [0.0]

    end
    #fitting all model with change points
    if n_max_change_points > 0

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
                n_max_change_points;
                type_of_detection=type_of_detection,
                type_of_curve=type_of_curve,
                pt_derivative=pt_smooth_derivative,
                size_win=win_size,
                method=method_peaks_detection,
                number_of_bin=n_bins,
            )

            combination_to_test = generation_of_combination_of_cpds(list_change_points_dev[2],
                n_fix=0)


        elseif fixed_cpd == true && detect_number_cpd == false

            list_change_points_dev = cpd_local_detection(
                data_testing_1,
                n_max_change_points;
                type_of_detection=type_of_detection,
                type_of_curve=type_of_curve,
                pt_derivative=pt_smooth_derivative,
                size_win=win_size,
                method=method_peaks_detection,
                number_of_bin=n_bins,
            )

            combination_to_test = generation_of_combination_of_cpds(list_change_points_dev[2],
                n_fix=n_max_change_points)



        elseif fixed_cpd == false && detect_number_cpd == false
            list_change_points_dev = cpd_local_detection(
                data_testing_1,
                n_max_change_points + 2;
                type_of_detection=type_of_detection,
                type_of_curve=type_of_curve,
                pt_derivative=pt_smooth_derivative,
                size_win=win_size,
                method=method_peaks_detection,
                number_of_bin=n_bins,
            )

            combination_to_test = generation_of_combination_of_cpds(list_change_points_dev[2],
                n_fix=n_max_change_points)


        end


        for i in 1:eachindex(combination_to_test)[end]

            cpd_temp = sort(combination_to_test[i])

            direct_search_results = selection_ODE_fixed_intervals(
                data_testing, # dataset first row times second row OD
                name_well, # name of the well
                label_exp, #label of the experiment
                list_of_models, # ode models to use
                param_array, # lower bound param
                cpd_temp;
                ub_param_array=ub_param_array, # upper bound param
                lb_param_array=lb_param_array, # lower bound param
                type_of_loss=type_of_loss, # type of used loss
                optimizer=optimizer, # selection of optimization method
                integrator=integrator, # selection of sciml integrator
                smoothing=smoothing,
                pt_avg=pt_avg,
                pt_smooth_derivative=pt_smooth_derivative,
                multiple_scattering_correction=multiple_scattering_correction, # if true uses the given calibration curve to fix the data
                method_multiple_scattering_correction=method_multiple_scattering_correction,
                calibration_OD_curve=calibration_OD_curve, #  the path to calibration curve to fix the data
                beta_smoothing_ms=beta_smoothing_ms, #  parameter of the AIC penality
                type_of_smoothing=type_of_smoothing,
                thr_lowess=thr_lowess,
                multistart=multistart,
                n_restart=n_restart,
                auto_diff_method=auto_diff_method,
                cons=cons,
                opt_params...
            )

            # composing piecewise penality

            n_param =
                sum([
                    length(direct_search_results[1][kk][4:(end-3)]) for
                    kk = 1:length(direct_search_results[1])
                ])

            n_param = n_param + length(cpd_temp)


            new_penality = AICc_evaluation2(n_param, beta_smoothing_ms, data_testing[2, :], direct_search_results[end], correction=correction_AIC)


            if new_penality <= score_of_the_models
                score_of_the_models = copy(new_penality)
                top_model = copy(direct_search_results[1])
                time_points_to_plot = copy(direct_search_results[3])
                sol_to_plot = copy(direct_search_results[4])
                change_point_to_plot = copy(direct_search_results[2])
                top_cps = copy(cpd_temp)
            end

            if save_all_model == true
                CSV.write(
                    string(
                        path_to_results,
                        label_exp,
                        "_segmented_ODE_",
                        name_well, "_conf_", i,
                        "_seg_",
                        length(cpd_temp) + 1,
                        ".csv",
                    ),
                    Tables.table(Vector(direct_search_results[1])),
                )
                CSV.write(
                    string(
                        path_to_results,
                        label_exp,
                        "_segmented_ODE_solution_",
                        name_well, "_conf_", i,
                        "_seg_",
                        length(cpd_temp) + 1,
                        ".csv",
                    ),
                    Tables.table(Vector(direct_search_results[4])),
                )
            end
        end
    end




    KinBiont_res_segmentation_ODE = ("ODE_segmentation", top_model, sol_to_plot, time_points_to_plot, top_cps, score_of_the_models)
    return KinBiont_res_segmentation_ODE
end

"""
    segment_gr_analysis(
    data::Matrix{Float64},
    name_well::String,
    label_exp::String;
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
    win_size=14,
    n_bins=40,
    method_peaks_detection="peaks_prominence"
    )

This function performs segmentation analysis on a time-series dataset by detecting change points. For each segment identified, it evaluates the minimum and maximum growth rate, the minimum and maximum derivative, the delta OD of the segment, and the times of change points.

# Arguments:

- `data::Matrix{Float64}`: The matrix containing the growth curve data. Time values should be in the first row, and the observable (e.g., OD) should be in the second row.
- `name_well::String`: Name of the well or sample under study.
- `label_exp::String`: Label or identifier for the experiment.

# Key Arguments:

- `n_max_change_points::Int`: Maximum number of change points to consider. If set to 0, the function will determine the number of change points based on the detection method and other parameters.
- `type_of_smoothing="rolling_avg"`: Method for smoothing the data. Options include `"NO"` (no smoothing), `"rolling_avg"` (rolling average), and `"lowess"` (locally weighted scatterplot smoothing).
- `pt_avg=7`: Size of the rolling average window for smoothing, applicable if `type_of_smoothing` is `"rolling_avg"`.
- `pt_smoothing_derivative=7`: Number of points used for the evaluation of specific growth rate. If less than 2, interpolation is used; otherwise, a sliding window approach is applied.
- `smoothing=false`: Boolean flag to enable or disable smoothing. Set to `true` to apply smoothing, or `false` to skip smoothing.
- `thr_lowess=0.05`: Parameter for lowess smoothing if `type_of_smoothing` is `"lowess"`.
- `type_of_detection="slinding_win"`: Method for detecting change points. Options include `"slinding_win"` (sliding window approach) and `"lsdd"` (least squares density difference).
- `type_of_curve="original"`: Specifies the input curve for change point detection. Options include `"original"` (for the raw time series) and `"deriv"` (for the specific growth rate time series).
- `method_peaks_detection="peaks_prominence"`: Method for detecting peaks in the dissimilarity curve. Options include `"peaks_prominence"` (orders peaks by prominence) and `"thr_scan"` (uses a threshold to select peaks).
- `n_bins=40`: Number of bins used to generate thresholds if `method_peaks_detection` is `"thr_scan"`.
- `win_size=14`: Size of the window used by the change point detection algorithms.

# Output:

If `res = segment_gr_analysis(...)`:

- `res[1]`: A string describing the method used for segmentation and analysis.
- `res[2]`: Array containing the parameters evaluated for each segment.
- `res[3]`: Intervals corresponding to the detected change points.
- `res[4]`: Preprocessed data, including smoothed values and calculated growth rates.

"""
function segment_gr_analysis(
    data::Matrix{Float64}, # dataset first row times second row OD
    name_well::String, # name of the well
    label_exp::String; #label of the experiment
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
    method_peaks_detection="peaks_prominence"
)


    if multiple_scattering_correction == true

        data = correction_OD_multiple_scattering(data, calibration_OD_curve; method=method_multiple_scattering_correction)
    end


    data = smoothing_data(
        data;
        method=type_of_smoothing,
        pt_avg=pt_avg,
        thr_lowess=thr_lowess
    )

    if n_max_change_points == 0



        res_temp = analyze_segment(label_exp, name_well, data, 0, pt_smoothing_derivative)
        res = res_temp[1]
        interval_changepoints = [0.0]
    else


        list_change_points_dev = cpd_local_detection(
            data,
            n_max_change_points;
            type_of_detection=type_of_detection,
            type_of_curve=type_of_curve,
            pt_derivative=pt_smoothing_derivative,
            size_win=win_size,
            method=method_peaks_detection,
            number_of_bin=n_bins,
        )

        combination_to_test = generation_of_combination_of_cpds(list_change_points_dev[2],
            n_fix=n_max_change_points)
        cpd_temp = sort(combination_to_test[1])


        interval_changepoints = copy(cpd_temp)
        interval_changepoints = push!(cpd_temp, data[1, 1])
        interval_changepoints = push!(cpd_temp, data[1, end])
        interval_changepoints = sort(cpd_temp)

        for i = 2:(eachindex(interval_changepoints)[end])
            if i == 2
                tspan_array = findall((data[1, :] .<= interval_changepoints[i]))
                data_temp = Matrix(
                    transpose(hcat(data[1, tspan_array], data[2, tspan_array])),)
                res_temp_1 = analyze_segment(label_exp, name_well, data_temp, i - 1, pt_smoothing_derivative)

                res = res_temp_1[1]

            else
                tspan_array_1 = findall((data[1, :] .> interval_changepoints[i-1]))
                tspan_array_2 = findall((data[1, :] .<= interval_changepoints[i]))
                tspan_array = intersect(tspan_array_1, tspan_array_2)
                data_temp = Matrix(
                    transpose(hcat(data[1, tspan_array], data[2, tspan_array])),
                )
                res_temp_1 = analyze_segment(label_exp, name_well, data_temp, i - 1, pt_smoothing_derivative)


                res = hcat(res, res_temp_1[1])


            end

        end

    end



    KinBiont_res_one_well = ("segment_analysis", res, interval_changepoints, data)



    return KinBiont_res_one_well
end


export fitting_one_well_Log_Lin
export fitting_one_well_ODE_constrained
export fitting_one_well_custom_ODE
export ODE_Model_selection
export one_well_morris_sensitivity
export selection_ODE_fixed_intervals
export segmentation_ODE
export segment_gr_analysis