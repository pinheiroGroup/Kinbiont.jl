using OptimizationBBO
using Distributions
using StatsBase
import Plots

#######################################################################
"""

    fitting_one_well_Log_Lin(
    data::Matrix{Float64},
    name_well::String,
    label_exp::String, 
    display_plots=false,
    save_plot=false,
    path_to_plot="NA", 
    type_of_smoothing="rolling_avg",
    pt_avg=7, 
    pt_smoothing_derivative=7, 
    pt_min_size_of_win=7, 
    type_of_win="maximum",
    threshold_of_exp=0.9, 
    multiple_scattering_correction=false, 
    method_multiple_scattering_correction="interpolation",
    calibration_OD_curve="NA", 
    thr_lowess=0.05 
    ) 
   
    
This function fits a logarithmic-linear model to the data of a single well to evaluate the specific growth rate ``\\lambda_s = \\frac{d}{dt}\\log(N(t))``. It uses a statistical 
threshold to define an exponetial window where the log-linear fit is performed.

# Arguments:
    
- `data::Matrix{Float64}`: The growth curve dataset. The first row represents the time, and the second row represents the fit observable (e.g., OD), see documentation.
- `name_well::String`: The name of the well.
- `label_exp::String`: The label of the experiment.

# Key Arguments:
- `path_to_plot="NA"`: String. Path to save the plots.
- `save_plot=false`: Bool. Options: "true" to save the plot, or "false" not to.
- `display_plots=true`: Bool. Options: "true" to display the plot in Julia, or "false" not to.
- `type_of_smoothing="rolling_avg"`: String. Method of choice to smooth the data. Options: "NO", "rolling_avg" (rolling average of the data), and "lowess".
- `pt_avg=7`: Int. Size of the rolling average window smoothing. 
- `pt_smoothing_derivative=7`: Int. Number of points for evaluation of the specific growth rate. If <2 it uses an interpolation algorithm. Otherwise, it uses a sliding window approach.
- `pt_min_size_of_win=7`: Int. The minimum size of the exponential windows in the number of smoothed points.
- `type_of_win="maximum"`: String. How the exponential phase window is selected ("maximum" or "global_thr").
- `threshold_of_exp=0.9`: Float. The quantile threshold in growth rate defininf the window where the log-linear fit is performed (a value between 0 and 1).
- `multiple_scattering_correction=false`: Bool. Options: "true" to perform the multiple scattering correction (requires a callibration curve) or "false" not to. 
- `calibration_OD_curve="NA"`: String. The path to the calibration curve (a .csv file). Used only if `multiple_scattering_correction=true`.
- `method_multiple_scattering_correction="interpolation"`: String. Method of choice to perform the multiple scattering curve inference. Options: '"interpolation"' or '"exp_fit"' (adapted from Meyers, A., Furtmann, C., & Jose, J., *Enzyme and microbial technology*, 118, 1-5., 2018). 
- `thr_lowess=0.05`: Float64 keyword argument of lowees smoothing.
- `start_exp_win_thr=0.05` minimum value (of OD) to consider the start of exp window

# Output: 

- `Matrix{Float64}`. An array with the following contents:

  `results_lin_log_fit = [label_exp, name_well, start of exp win, end of exp win, start of exp win, Maximum specific GR, specific GR, 2 sigma  CI of GR, doubling time,doubling time - 2 sigma, doubling time + 2 sigma, intercept log-lin fitting, 2 sigma intercept, R^2]`

- The plots of the log-linear fitting and the specific growth rate dynamics (if `save_plot=true` or `display_plots=true`).

"""
function fitting_one_well_Log_Lin(
    data::Matrix{Float64}, # dataset first row times second row OD
    name_well::String, # name of the well
    label_exp::String; #label of the experiment
    display_plots=false, # do plots or no
    save_plot=false,
    path_to_plot="NA", # where save plots
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


        index_of_max = argmax(specific_gr)[1]

        index_gr_max = findlast(x -> x > lb_of_distib, specific_gr[1:end])[1]
        
        index_gr_min = findfirst(x -> x > lb_of_distib, specific_gr[1:end])[1]

        while data_smooted[2,index_gr_min] <  start_exp_win_thr

            index_gr_min = index_gr_min +1
        end   



        t_start = specific_gr_times[index_gr_min]
        t_end = specific_gr_times[index_gr_max]
        
        index_of_t_start = findfirst(x -> x > t_start, data_smooted[1, :])[1]
        index_of_t_end = findall(x -> x > t_end, data_smooted[1, :])[1]



    end

    # checking the minimum size of the window before fitting
    if (index_of_t_end - index_of_t_start) < pt_min_size_of_win
        index_of_t_start = convert(Int, index_of_max - floor(pt_min_size_of_win / 2))
        index_of_t_end = convert(Int, index_of_max + floor(pt_min_size_of_win / 2))

        if index_of_t_start < 1
            index_of_t_start = 2
        end

        if index_of_t_end > length(data_smooted[1, :])
            index_of_t_end = length(data_smooted[1, :]) - 1
        end
    end

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

    if display_plots
        if_display = display
    else
        if_display = identity
    end

    if save_plot
        mkpath(path_to_plot)
    end

    # plotting if requested
    if_display(
        Plots.scatter(
            data_smooted[1, :],
            log.(data_smooted[2, :]),
            xlabel="Time",
            ylabel="Log(Arb. Units)",
            label=["Data " nothing],
            markersize=1,
            color=:black,
            title=string(label_exp, " ", name_well),
        ),
    )
    if_display(
        Plots.plot!(
            data_to_fit_times,
            Theoretical_fitting,
            ribbon=confidence_band,
            xlabel="Time ",
            ylabel="Log(Arb. Units)",
            label=[string("Fitting Log-Lin ") nothing],
            c=:red,
        ),
    )
    if_display(
        Plots.vline!(
            [data_to_fit_times[1], data_to_fit_times[end]],
            c=:black,
            label=[string("Window of exp. phase ") nothing],
        ),
    )
    if save_plot
        png(string(path_to_plot, label_exp, "_Log_Lin_Fit_", name_well, ".png"))
    end
    if_display(
        Plots.scatter(
            specific_gr_times,
            specific_gr,
            xlabel="Time ",
            ylabel="1 /time ",
            label=[string("Dynamics growth rate ") nothing],
            c=:red,
        ),
    )
    if_display(
        Plots.vline!(
            [data_to_fit_times[1], data_to_fit_times[end]],
            c=:black,
            label=[string("Window of exp. phase ") nothing],
        ),
    )
    if save_plot
        png(string(path_to_plot, label_exp, "_dynamics_gr_", name_well, ".png"))
    end

    return results_lin_log_fit
end


"""
    fitting_one_well_ODE_constrained(
    data::Matrix{Float64},
    name_well::String, 
    label_exp::String, 
    model::String, 
    lb_param::Vector{Float64},
    ub_param::Vector{Float64},
    param=lb_param .+ (ub_param .- lb_param) ./ 2,
    optmizator=BBO_adaptive_de_rand_1_bin_radiuslimited(), 
    integrator=Tsit5(), 
    display_plots=true, 
    save_plot=false,
    path_to_plot="NA", 
    pt_avg=1, 
    pt_smooth_derivative=7,
    smoothing=false,
    type_of_smoothing="rolling_avg",
    type_of_loss="RE",
    blank_array=zeros(100), 
    multiple_scattering_correction=false,
    method_multiple_scattering_correction="interpolation",
    calibration_OD_curve="NA",  
    PopulationSize=300,
    maxiters=2000000,
    abstol=0.00001,
    thr_lowess=0.05
    )

This function uses an ordinary differential equation (ODE) model to fit the data of a on a single well. It estimates the model parameters within specified lower and upper bounds.

# Arguments:

- `data::Matrix{Float64}`: The growth curve data. Time values are in the first row and the fit observable (e.g., OD) is in the second row, see documentation.
- `model::String`: ODE model of choice.
- `name_well::String`: Name of the well.
- `label_exp::String`: Label of the experiment.
- `lb_param::Vector{Float64}`: Lower bounds of the model parameters.
- `ub_param::Vector{Float64}`: Upper bounds of the model parameters.

# Key Arguments:

- `param=lb_param .+ (ub_param.-lb_param)./2`: Vector{Float64}. Used as the default initial guess for the model parameters.
- `integrator=Tsit5()`: sciML integrator. Use 'KenCarp4(autodiff=true)' to fit piecewise models.
- `optmizator=BBO_adaptive_de_rand_1_bin_radiuslimited()`: Optimizer from optimizationBBO.
- `save_plot=false`: Bool. Options: "true" to save the plot, or "false" not to.
- `display_plots=true`: Bool. Options: "true" to display the plot, or "false" not to.
- `type_of_smoothing="rolling_avg"`: String. Method of choice to smooth the data. Options: "NO", "rolling_avg" (rolling average of the data), and "lowess".
- `pt_avg=7`: Int. Size of the rolling average window smoothing. 
- `smoothing=false`: Bool. Options: "true" to smooth the data, or "false" not to.
- `type_of_loss:="RE" `: Type of loss function to be used. Options = "RE" (relative error), "L2" (L2 norm), "L2_derivative" (Xx) and "blank_weighted_L2" (Xx).
- `blank_array=zeros(100)`: Data of all blanks in single array.
- `pt_smoothing_derivative=7`: Int. Number of points for evaluation of specific growth rate. If <2 it uses interpolation algorithm otherwise a sliding window approach.
- `calibration_OD_curve="NA"`: String, The path where the .csv calibration data are located, used only if `multiple_scattering_correction=true`.
- `multiple_scattering_correction=false`: Bool. Options: "true" to perform the multiple scattering correction (requires a callibration curve) or "false" not to. 
- `calibration_OD_curve="NA"`: String. The path to the calibration curve (a .csv file). Used only if `multiple_scattering_correction=true`.
- `method_multiple_scattering_correction="interpolation"`: String. Method of choice to perform the multiple scattering curve inference. Options: '"interpolation"' or '"exp_fit"' (adapted from Meyers, A., Furtmann, C., & Jose, J., *Enzyme and microbial technology*, 118, 1-5., 2018). 
- `thr_lowess=0.05`: Float64. Argument of the lowess smoothing.
- `PopulationSize=100`: Size of the population of the optimization (Xx)
- `maxiters=2000000`: Stop criterion, the optimization stops when the number of iterations is bigger than `maxiters`.
- `abstol=0.00001`: Stop criterion, the optimization stops when the loss is smaller than `abstol`.


# Output (if `results_ODE_fit=fitting_one_well_ODE_constrained(...)`):

- `results_ODE_fit[1]`. An array with the following contents: 

  `["name of model", "well", "param_1", "param_2",..,"param_n", "maximum specific gr using ODE", "maximum specific gr using data", "objective function value (i.e. loss of the solution)"]`,
where `"param_1", "param_2", .., "param_n"` are the ODE model fit parameters as in the documentation.

- `results_ODE_fit[2]`. The time coordinates (Xx) of the fitted ODE. 

- `results_ODE_fit[3]`. The numerical solution of the fitted ODE.

- The fit plot if `do_plot=true`.

"""
function fitting_one_well_ODE_constrained(
    data::Matrix{Float64}, # dataset first row times second row OD
    name_well::String, # name of the well
    label_exp::String, #label of the experiment
    model::String, # ode model to use
    lb_param::Vector{Float64}, # lower bound param
    ub_param::Vector{Float64}; # upper bound param
    param=lb_param .+ (ub_param .- lb_param) ./ 2,# initial guess param
    optmizator=BBO_adaptive_de_rand_1_bin_radiuslimited(), # selection of optimization method
    integrator=Tsit5(), # selection of sciml integrator
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
    optf = Optimization.OptimizationFunction((x, p) -> loss_function(x))
    optprob_const =
        Optimization.OptimizationProblem(optf, param, u0, lb=lb_param, ub=ub_param)
    res = Optimization.solve(
        optprob_const,
        optmizator,
        PopulationSize=PopulationSize,
        maxiters=maxiters,
        abstol=abstol,
    )

    #revalution of solution for plot an loss evaluation
    remade_solution = solve(remake(ODE_prob, p=res.u), integrator, saveat=tsteps)
    sol_time = reduce(hcat, remade_solution.t)
    sol_fin = reduce(hcat, remade_solution.u)
    sol_fin = sum(sol_fin, dims=1)

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
        ),
    )
    if_display(
        Plots.plot!(
            remade_solution.t,
            sol_fin[1, 1:end],
            xlabel="Time",
            ylabel="Arb. Units",
            label=[string("Fitting ", model) nothing],
            c=:red,
        ),
    )
    if save_plot
        png(string(path_to_plot, label_exp, "_", model, "_", name_well, ".png"))
    end

    # max_theoretical gr
    sol_fin, index_not_zero = remove_negative_value(sol_fin)

    data_th = transpose(hcat(sol_time[index_not_zero], sol_fin))


    max_th_gr = maximum(specific_gr_evaluation(Matrix(data_th), pt_smooth_derivative))

    # max empirical gr
    max_em_gr = maximum(specific_gr_evaluation(data, pt_smooth_derivative))
    res_temp = res.u
    loss_value = res.objective

    res_param =
        vectorize_df_results(name_well, model, res_temp, max_th_gr, max_em_gr, loss_value)

    return res_param, remade_solution.t, sol_fin
end

#######################################################################
#######################################################################
"""
    fitting_one_well_custom_ODE(
    data::Matrix{Float64},
    name_well::String, 
    label_exp::String, 
    model::Any, 
    lb_param::Vector{Float64},
    ub_param::Vector{Float64},
    n_equation::Int,
    param=lb_param .+ (ub_param .- lb_param) ./ 2,
    optmizator=BBO_adaptive_de_rand_1_bin_radiuslimited(), 
    integrator=Tsit5(), 
    display_plots=false, 
    save_plot=false,
    path_to_plot="NA", 
    pt_avg=1, 
    pt_smooth_derivative=0,
    smoothing=false, 
    type_of_loss="RE",
    blank_array=zeros(100), 
    multiple_scattering_correction=false, 
    method_multiple_scattering_correction="interpolation",
    calibration_OD_curve="NA",  
    PopulationSize=300,
    maxiters=2000000,
    abstol=0.00001,
    thr_lowess=0.05,
    type_of_smoothing="lowess"
    )

This function is designed to fit a user-defined ordinary differential equation (ODE) model to time-series data of a single well.

# Arguments:

- `data::Matrix{Float64}`: The growth curve data. Time values are in the first row and the fit observable (e.g., OD) is in the second row, see documentation.
- `model::Any`: Function of the ODE model to be fitted. See the documentation for examples.
- `name_well::String`: Name of the well.
- `label_exp::String`: Label of the experiment.
- `lb_param::Vector{Float64}`: Lower bounds of the model parameters.
- `ub_param::Vector{Float64}`: Upper bounds of the model parameters.
- `n_equation::Int`: Number of ODEs in the model.

# Key Arguments:

- `param=lb_param .+ (ub_param.-lb_param)./2`: Vector{Float64}. Used as the default initial guess for the model parameters.
- `integrator=Tsit5()`: sciML integrator. Use 'KenCarp4(autodiff=true)' to fit piecewise models.
- `optmizator=BBO_adaptive_de_rand_1_bin_radiuslimited()`: Optimizer from optimizationBBO.
- `save_plot=false`: Bool. Options: "true" to save the plot, or "false" not to.
- `display_plots=true`: Bool. Options: "true" to display the plot, or "false" not to.
- `type_of_smoothing="rolling_avg"`: String. Method of choice to smooth the data. Options: "NO", "rolling_avg" (rolling average of the data), and "lowess".
- `pt_avg=7`: Int. Size of the rolling average window smoothing. 
- `smoothing=false`: Bool. Options: "true" to smooth the data, or "false" not to.
- `type_of_loss:="RE" `: Type of loss function to be used. Options = "RE" (relative error), "L2" (L2 norm), "L2_derivative" (Xx) and "blank_weighted_L2" (Xx).
- `blank_array=zeros(100)`: Data of all blanks in single array.
- `pt_smoothing_derivative=7`: Int. Number of points for evaluation of specific growth rate. If <2 it uses interpolation algorithm otherwise a sliding window approach.
- `calibration_OD_curve="NA"`: String, The path where the .csv calibration data are located, used only if `multiple_scattering_correction=true`.
- `multiple_scattering_correction=false`: Bool. Options: "true" to perform the multiple scattering correction (requires a callibration curve) or "false" not to. 
- `calibration_OD_curve="NA"`: String. The path to the calibration curve (a .csv file). Used only if `multiple_scattering_correction=true`.
- `method_multiple_scattering_correction="interpolation"`: String. Method of choice to perform the multiple scattering curve inference. Options: '"interpolation"' or '"exp_fit"' (adapted from Meyers, A., Furtmann, C., & Jose, J., *Enzyme and microbial technology*, 118, 1-5., 2018). 
- `thr_lowess=0.05`: Float64. Argument of the lowess smoothing.
- `PopulationSize=100`: Size of the population of the optimization (Xx)
- `maxiters=2000000`: Stop criterion, the optimization stops when the number of iterations is bigger than `maxiters`.
- `abstol=0.00001`: Stop criterion, the optimization stops when the loss is smaller than `abstol`.


# Output (if `results_ODE_fit =fitting_one_well_custom_ODE(...)`):
- `results_ODE_fit[1]`. An array with the following contents: 

  `["name of model", "well", "param_1", "param_2",..,"param_n", "maximum specific gr using ODE", "maximum specific gr using data", "objective function value (i.e. loss of the solution)"]`,
where `"param_1", "param_2", .., "param_n"` are the ODE model fit parameters as in the documentation.

- The fit plot if `do_plot=true`.

"""
function fitting_one_well_custom_ODE(
    data::Matrix{Float64}, # dataset first row times second row OD
    name_well::String, # name of the well
    label_exp::String, #label of the experiment
    model::Any, # ode model to use
    lb_param::Vector{Float64}, # lower bound param
    ub_param::Vector{Float64}, # upper bound param
    n_equation::Int; # number ode in the system
    param=lb_param .+ (ub_param .- lb_param) ./ 2,# initial guess param
    optmizator=BBO_adaptive_de_rand_1_bin_radiuslimited(), # selection of optimization method
    integrator=Tsit5(), # selection of sciml integrator
    display_plots=false, # do plots or no
    save_plot=false,
    path_to_plot="NA", # where save plots
    pt_avg=1, # numebr of the point to generate intial condition
    pt_smooth_derivative=0,
    smoothing=false, # the smoothing is done or not?
    type_of_loss="RE", # type of used loss
    blank_array=zeros(100), # data of all blanks
    multiple_scattering_correction=false, # if true uses the given calibration curve to fix the data
    method_multiple_scattering_correction="interpolation",
    calibration_OD_curve="NA",  #  the path to calibration curve to fix the data
    PopulationSize=300,
    maxiters=2000000,
    abstol=0.00001,
    thr_lowess=0.05,
    type_of_smoothing="lowess",
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
    optf = Optimization.OptimizationFunction((x, p) -> loss_function(x))
    optprob_const =
        Optimization.OptimizationProblem(optf, param, u0, lb=lb_param, ub=ub_param)
    res = Optimization.solve(
        optprob_const,
        optmizator,
        PopulationSize=PopulationSize,
        maxiters=maxiters,
        abstol=abstol,
    )

    #revalution of solution for plot an loss evaluation
    remade_solution = solve(remake(ODE_prob, p=res.u), integrator, saveat=tsteps)
    sol_time = reduce(hcat, remade_solution.t)
    sol_fin = reduce(hcat, remade_solution.u)
    sol_fin = sum(sol_fin, dims=1)

    if display_plots
        if_display = display
    else
        if_display = identity
    end

    if save_plot
        mkpath(path_to_plot)
    end

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
        ),
    )
    if_display(
        Plots.plot!(
            remade_solution.t,
            sol_fin[1, 1:end],
            xlabel="Time",
            ylabel="Arb. Units",
            label=[string("Fitting custom model") nothing],
            c=:red,
        ),
    )
    if save_plot
        png(string(path_to_plot, label_exp, "_custom_model_", name_well, ".png"))
    end

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

    return res_param
end

#######################################################################

"""

    ODE_Model_selection(
    data::Matrix{Float64},
    name_well::String,
    label_exp::String, 
    models_list::Vector{String}, 
    lb_param_array::Any, 
    ub_param_array::Any,
    optmizator=BBO_adaptive_de_rand_1_bin_radiuslimited(), 
    integrator=Tsit5(), 
    pt_avg=1,
    beta_penality=2.0,
    smoothing=false,
    type_of_smoothing="lowess",
    thr_lowess=0.05,
    type_of_loss="L2",
    blank_array=zeros(100),
    display_plot_best_model=false, 
    save_plot_best_model=false,
    path_to_plot="NA",
    pt_smooth_derivative=7,
    multiple_scattering_correction=false, 
    method_multiple_scattering_correction="interpolation",
    calibration_OD_curve="NA", 
    verbose=false,
    PopulationSize=300,
    maxiters=2000000,
    abstol=0.00001,
    correction_AIC=true
    )

Automatic model selection for multiple ODE model fits in the time series of a single well. The best-fitting model is chosen on the basis of the Akaike Information Criterion (AIC) or corrected AIC (AICc).

# Arguments:
- `data::Matrix{Float64}`: The growth curve data. Time values are in the first row and the fit observable (e.g., OD) is in the second row, see documentation.
- `name_well::String`: Name of the well.
- `label_exp::String`: Label of the experiment.
- `models_list::Vector{String}`: A vector of ODE models to evaluate.
- `lb_param_array::Any`: Lower bounds for the parameters (compatible with the models).
- `ub_param_array::Any`: Upper bounds for the parameters (compatible with the models).

# Key Arguments:

- `param=lb_param .+ (ub_param.-lb_param)./2`: Vector{Float64}. Used as the default initial guess for the model parameters.
- `integrator=Tsit5()`: sciML integrator. Use 'KenCarp4(autodiff=true)' to fit piecewise models.
- `optmizator=BBO_adaptive_de_rand_1_bin_radiuslimited()`: Optimizer from optimizationBBO.
- `save_plot=false`: Bool. Options: "true" to save the plot, or "false" not to.
- `display_plots=true`: Bool. Options: "true" to display the plot, or "false" not to.
- `type_of_smoothing="rolling_avg"`: String. Method of choice to smooth the data. Options: "NO", "rolling_avg" (rolling average of the data), and "lowess".
- `pt_avg=7`: Int. Size of the rolling average window smoothing. 
- `pt_smoothing_derivative=7`: Int. Number of points for evaluation of the specific growth rate. If <2 it uses interpolation algorithm otherwise a sliding window approach.
- `smoothing=false`: Bool. Options: "true" to smooth the data, or "false" not to.
- `type_of_loss:="RE" `: Type of loss function to be used. Options = "RE" (relative error), "L2" (L2 norm), "L2_derivative" (Xx) and "blank_weighted_L2" (Xx).
- `blank_array=zeros(100)`: Data of all blanks in single array.
- `calibration_OD_curve="NA"`: String, The path where the .csv calibration data are located, used only if `multiple_scattering_correction=true`.
- `multiple_scattering_correction=false`: Bool. Options: "true" to perform the multiple scattering correction (requires a callibration curve) or "false" not to. 
- `method_multiple_scattering_correction="interpolation"`: String. Method of choice to perform the multiple scattering curve inference. Options: '"interpolation"' or '"exp_fit"' (adapted from Meyers, A., Furtmann, C., & Jose, J., *Enzyme and microbial technology*, 118, 1-5., 2018). 
- `thr_lowess=0.05`: Float64. Argument of the lowess smoothing.
- `PopulationSize=100`: Size of the population of the optimization (Xx).
- `maxiters=2000000`: stop criterion, the optimization stops when the number of iterations is bigger than `maxiters`.
- `abstol=0.00001`: stop criterion, the optimization stops when the loss is smaller than `abstol`.
- `correction_AIC=true`: Bool. Options: "true" to perform the AIC finite samples correction or "false" not to.
- `beta_penality=2.0`: Penality parameters for the evaluation of AIC (or AICc).



# Output (if `Model_selection =ODE_Model_selection(...)`):

- `Model_selection[1]`: Matrix containing the loss and the AIC score for each model.
- `Model_selection[2]`: Tuple containing all the fitted models.
- `Model_selection[3]`: The best model's AIC score.
- `Model_selection[4]`: The best model's loss value. 
- `Model_selection[5]`: The best model's parameters. 
- `Model_selection[6]`: The best model's name. 
- `Model_selection[7]`: The fitted ODE numerical value. 
- The best model fit plot if `save_plot_best_model=true` or `display_plot_best_model=true` .
"""
function ODE_Model_selection(
    data::Matrix{Float64}, # dataset first row times second row OD
    name_well::String, # name of the well
    label_exp::String, #label of the experiment
    models_list::Vector{String}, # ode model to use
    lb_param_array::Any, # lower bound param
    ub_param_array::Any; # upper bound param
    optmizator=BBO_adaptive_de_rand_1_bin_radiuslimited(), # selection of optimization method
    integrator=Tsit5(), # selection of sciml integrator
    pt_avg=3, # number of the point to generate intial condition
    beta_penality=2.0, # penality for AIC evaluation
    smoothing=false, # the smoothing is done or not?
    type_of_smoothing="rolling_avg",
    thr_lowess=0.05,
    type_of_loss="L2", # type of used loss
    blank_array=zeros(100), # data of all blanks
    display_plot_best_model=false, # one wants the results of the best fit to be plotted
    save_plot_best_model=false,
    path_to_plot="NA",
    pt_smooth_derivative=7,
    multiple_scattering_correction=false, # if true uses the given calibration curve to fix the data
    method_multiple_scattering_correction="interpolation",
    calibration_OD_curve="NA", #  the path to calibration curve to fix the data
    verbose=false,
    PopulationSize=300,
    maxiters=2000000,
    abstol=0.00001,
    correction_AIC=true,)
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
        temp_param_lb = lb_param_array[mm]
        temp_param_ub = ub_param_array[mm]
        temp_start_param = temp_param_lb .+ (temp_param_ub .- temp_param_lb) ./ 2

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
        optf = Optimization.OptimizationFunction((x, p) -> loss_function(x))

        # solving optimization problem
        optprob_const = Optimization.OptimizationProblem(
            optf,
            temp_start_param,
            u0,
            lb=temp_param_lb,
            ub=temp_param_ub,
        )
        res = Optimization.solve(
            optprob_const,
            optmizator,
            PopulationSize=PopulationSize,
            maxiters=maxiters,
            abstol=abstol,
        )

        #revalution of solution for plot an loss evaluation
        remade_solution = solve(remake(ODE_prob, p=res.u), integrator, saveat=tsteps)
        sol_time = reduce(hcat, remade_solution.t)
        sol_fin = reduce(hcat, remade_solution.u)
        sol_fin = sum(sol_fin, dims=1)

        param_number = length(temp_start_param)

        #  AICc = AICc_evaluation(param_number, beta_penality, data[2, :], sol_fin, correction=correction_AIC)
        AICc = AICc_evaluation2(param_number, beta_penality, data[2, :], res.objective, correction=correction_AIC)

        #max_theoretical gr
        sol_fin, index_not_zero = remove_negative_value(sol_fin)


        data_th = transpose(hcat(sol_time[index_not_zero], sol_fin))
        max_th_gr = maximum(specific_gr_evaluation(Matrix(data_th), pt_smooth_derivative))

        # max empirical gr
        max_em_gr = maximum(specific_gr_evaluation(data, pt_smooth_derivative))
        res_temp = res.u
        loss_value = res.objective

        res_param = vectorize_df_results(
            name_well,
            temp_model,
            res_temp,
            max_th_gr,
            max_em_gr,
            loss_value,
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
    param_min = param_min[3:(end-3)]

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
    if display_plot_best_model
        if_display = display
    else
        if_display = identity
    end

    if save_plot_best_model
        mkpath(path_to_plot)
    end

    sol_fin, index_not_zero = remove_negative_value(sol_fin)

    data_th = transpose(hcat(sol_time[index_not_zero], sol_fin))
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
        ),
    )
    if_display(
        Plots.plot!(
            data_th[1, :],
            data_th[2, :],
            xlabel="Time",
            ylabel="Arb. Units",
            label=[string("Fitting ", model) nothing],
            c=:red,
        ),
    )
    if save_plot_best_model
        png(string(path_to_plot, label_exp, "_", model, "_", name_well, ".png"))
    end

    return rss_array,
    df_res_optimization,
    min_AIC,
    minimum(rss_array[2, 2:end]),
    param_min,
    model,
    data_th,
    param_out_full
end



#######################################################################

"""
    one_well_morris_sensitivity(
    data::Matrix{Float64},
    name_well::String, 
    label_exp::String, 
    model::String,
    lb_param::Vector{Float64}, 
    ub_param::Vector{Float64},
    N_step_morris=7,
    optmizator=BBO_adaptive_de_rand_1_bin_radiuslimited(), 
    integrator=Tsit5(),
    pt_avg=1, 
    pt_smooth_derivative=7,
    write_res=false,
    smoothing=false, 
    type_of_smoothing="lowess",
    type_of_loss="RE", 
    blank_array=zeros(100), 
    multiple_scattering_correction=false, 
    method_multiple_scattering_correction="interpolation",
    calibration_OD_curve="NA", 
    PopulationSize=300,
    maxiters=2000000,
    abstol=0.00001
    )

This function performs the Morris sensitivity analysis, which assesses the sensitivity of the fit parameters to variations of the initial guess (suitable for quality checks of nonlinear model fits. See https://docs.sciml.ai/GlobalSensitivity/stable/methods/morris/). 

# Arguments:

- `data::Matrix{Float64}`: The growth curve data. Time values are in the first row and the fit observable (e.g., OD) is in the second row, see documentation.
- `name_well::String`: Name of the well.
- `label_exp::String`: Label of the experiment.
- `model::String`: The ODE model of choice.
- `lb_param::Vector{Float64}`: Lower bounds for the parameters.
- `ub_param::Vector{Float64}`: Upper bounds for the parameters.


# Key Arguments:

- `N_step_morris=7`: Number of steps for the Morris sensitivity analysis.
- `param=lb_param .+ (ub_param.-lb_param)./2`: Vector{Float64}. Initial guess for the model parameters.
- `integrator=Tsit5()`: sciML integrator. Use 'KenCarp4(autodiff=true)' to fit piecewise models.
- `optmizator=BBO_adaptive_de_rand_1_bin_radiuslimited()`: Optimizer from optimizationBBO.
- `display_plots=true`: Bool. Options: "true" to display the plot, or "false" not to.
- `type_of_smoothing="rolling_avg"`: String. Method of choice to smooth the data. Options: "NO", "rolling_avg" (rolling average of the data), and "lowess".
- `pt_avg=7`: Int. Size of the rolling average window smoothing. 
- `pt_smoothing_derivative=7`: Int. Number of points for evaluation of specific growth rate. If <2 it uses interpolation algorithm otherwise a sliding window approach.
- `smoothing=false`: Bool. Options: "true" to smooth the data, or "false" not to.
- `type_of_loss:="RE" `: Type of loss function to be used. Options = "RE" (relative error), "L2" (L2 norm), "L2_derivative" (Xx) and "blank_weighted_L2" (Xx).
- `blank_array=zeros(100)`: Data of all blanks in a single array.
- `calibration_OD_curve="NA"`: String. The path where the .csv calibration data are located, used only if `multiple_scattering_correction=true`.
- `multiple_scattering_correction=false`: Bool. Options: "true" to perform the multiple scattering correction (requires a callibration curve) or "false" not to. 
- `method_multiple_scattering_correction="interpolation"`: String. Method of choice to perform the multiple scattering curve inference. Options: '"interpolation"' or '"exp_fit"' (adapted from Meyers, A., Furtmann, C., & Jose, J., *Enzyme and microbial technology*, 118, 1-5., 2018). 
- `thr_lowess=0.05`: Float64. Argument of the lowess smoothing.
- `PopulationSize=100`: Size of the population of the optimization (Xx).
- `maxiters=2000000`: Stop criterion, the optimization stops when the number of iterations is bigger than `maxiters`.
- `abstol=0.00001`: Stop criterion, the optimization stops when the loss is smaller than `abstol`.


# Output (if `results_ODE_morris_sensitivity =one_well_morris_sensitivity(...)`):

- `results_ODE_morris_sensitivity[1]`. A matrix with the optimzation parameters' initial guess in each column; same parameter order as in this [table](#ODE_list).
- `results_ODE_morris_sensitivity[2]`. A matrix with the following contents for each column: 

`["name of model", "well", "param_1", "param_2",.., "param_n", "maximum specific gr using ode", "maximum specific gr using data", "objective function value (i.e. loss of the solution)"]`,

where  `"param_1","param_2",..,"param_n"` are the parameters of the selected ODE as in this [table](#ODE_list). It can be saved into a .csv if `write_res=true`.

- The best fitting model plot if `save_plot_best_model=true` or `display_plot_best_model=true`.

"""
function one_well_morris_sensitivity(
    data::Matrix{Float64}, # dataset first row times second row OD
    name_well::String, # name of the well
    label_exp::String, #label of the experiment
    model::String, # ode model to use
    lb_param::Vector{Float64}, # lower bound param
    ub_param::Vector{Float64}; # upper bound param
    N_step_morris=7,
    optmizator=BBO_adaptive_de_rand_1_bin_radiuslimited(), # selection of optimization method
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
    PopulationSize=300,
    maxiters=2000000,
    abstol=0.00001,
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
    optf = Optimization.OptimizationFunction((x, p) -> loss_function(x))
    param_combination =
        generation_of_combination_of_IC_morris(lb_param, ub_param, N_step_morris)

    for i = 1:size(param_combination)[2]
        param = param_combination[:, i]
        optprob_const =
            Optimization.OptimizationProblem(optf, param, u0, lb=lb_param, ub=ub_param)
        res = Optimization.solve(
            optprob_const,
            optmizator,
            PopulationSize=PopulationSize,
            maxiters=maxiters,
            abstol=abstol,
        )

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

    return param_combination, results_sensitivity
end



"""
    selection_ODE_fixed_intervals(
    data_testing::Matrix{Float64}, 
    name_well::String, 
    label_exp::String, 
    list_of_models::Vector{String}, 
    list_lb_param::Any, 
    list_ub_param::Any, 
    intervals_changepoints::Any;
    type_of_loss="L2", 
    optmizator=BBO_adaptive_de_rand_1_bin_radiuslimited(), 
    integrator=Tsit5(), 
    smoothing=false,
    type_of_smoothing="lowess",
    thr_lowess=0.05,
    pt_avg=1,
    save_plot=false, 
    display_plots=false,
    path_to_plot="NA", 
    pt_smooth_derivative=7,
    multiple_scattering_correction=false, 
    method_multiple_scattering_correction="interpolation",
    calibration_OD_curve="NA", 
    beta_smoothing_ms=2.0, 
    PopulationSize=300,
    maxiters=2000000,
    abstol=0.0000000001,
    correction_AIC=true
    )
    

This function fits an ODE model at each segment of the time-series data. Change points are supplied by the user. 

# Arguments:

- `data_testing::Matrix{Float64}`:  The growth curve data. Time values are in the first row and the fit observable (e.g., OD) is in the second row, see documentation.
- `name_well::String`: Name of the well.
- `label_exp::String`: Label of the experiment.
- `list_of_models::Vector{String}`: List of the ODE models of choice.
- `list_lb_param::Any`: Lower bounds for the parameters (compatible with the models).
- `list_ub_param::Any`: Upper bounds for the parameters (compatible with the models).
- `intervals_changepoints::Any`: Array containing the list of change points, e.g., [0.0 10.0 30.0]. 


# Key Arguments:

- `param=lb_param .+ (ub_param.-lb_param)./2`: Vector{Float64}. Used as the default initial guess for the model parameters.
- `integrator=Tsit5()`: sciML integrator. Use 'KenCarp4(autodiff=true)' to fit piecewise models.
- `optmizator=BBO_adaptive_de_rand_1_bin_radiuslimited()`: Optimizer from optimizationBBO.
- `save_plot=false`: Bool. Options: "true" to save the plot, or "false" not to.
- `display_plots=true`: Bool. Options: "true" to display the plot, or "false" not to.
- `type_of_smoothing="rolling_avg"`: String. Method of choice to smooth the data. Options: "NO", "rolling_avg" (rolling average of the data), and "lowess".
- `pt_avg=7`: Int. Size of the rolling average window smoothing. 
- `smoothing=false`: Bool. Options: "true" to smooth the data, or "false" not to.
- `type_of_loss:="RE" `: Type of loss function to be used. Options = "RE" (relative error), "L2" (L2 norm), "L2_derivative" (Xx) and "blank_weighted_L2" (Xx).
- `blank_array=zeros(100)`: Data of all blanks in single array.
- `pt_smoothing_derivative=7`:Int. Number of points for evaluation of specific growth rate. If <2 it uses interpolation algorithm otherwise a sliding window approach.
- `calibration_OD_curve="NA"`: String, The path where the .csv calibration data are located, used only if `multiple_scattering_correction=true`.
- `multiple_scattering_correction=false`: Bool. Options: "true" to perform the multiple scattering correction (requires a callibration curve) or "false" not to. 
- `method_multiple_scattering_correction="interpolation"`: String. Method of choice to perform the multiple scattering curve inference. Options: '"interpolation"' or '"exp_fit"' (adapted from Meyers, A., Furtmann, C., & Jose, J., *Enzyme and microbial technology*, 118, 1-5., 2018). 
- `thr_lowess=0.05`: Float64. Argument of the lowess smoothing.
- `PopulationSize=100`: Size of the population of the optimization (Xx).
- `maxiters=2000000`: stop criterion, the optimization stops when the number of iterations is bigger than `maxiters`.
- `abstol=0.00001`: stop criterion, the optimization stops when the loss is smaller than `abstol`.
- `beta_penality=2.0` penality  parameters for AIC (or AICc) evaluation.

# Output (if `res =selection_ODE_fixed_intervals(...)`:

- `res[1]`. Parameters of each segment.
- `res[2]`. Interval of each ODE segment.
- `res[3]`. Time of the fitted solution.
- `res[4]`. Numerical value of the fitted solution.
- `res[5]`. The fit loss score. 
- The best fitting model plot if `save_plot_best_model=true` or `display_plot_best_model=true` .

"""
function selection_ODE_fixed_intervals(
    data_testing::Matrix{Float64}, # dataset first row times second row OD
    name_well::String, # name of the well
    label_exp::String, #label of the experiment
    list_of_models::Vector{String}, # ode models to use
    list_lb_param::Any, # lower bound param
    list_ub_param::Any, # upper bound param
    intervals_changepoints::Any;
    type_of_loss="L2", # type of used loss
    optmizator=BBO_adaptive_de_rand_1_bin_radiuslimited(), # selection of optimization method
    integrator=Tsit5(), # selection of sciml integrator
    smoothing=false,
    type_of_smoothing="lowess",
    thr_lowess=0.05,
    pt_avg=1,
    save_plot=false, # do plots or no
    display_plots=false,
    path_to_plot="NA", # where save plots
    pt_smooth_derivative=0,
    multiple_scattering_correction=false, # if true uses the given calibration curve to fix the data
    method_multiple_scattering_correction="interpolation",
    calibration_OD_curve="NA", #  the path to calibration curve to fix the data
    beta_smoothing_ms=2.0, #  parameter of the AIC penality
    PopulationSize=300,
    maxiters=2000000,
    abstol=0.0000000001,
    correction_AIC=true)

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
            list_lb_param, # lower bound param
            list_ub_param; # upper bound param
            optmizator=optmizator, # selection of optimization method
            integrator=integrator, # selection of sciml integrator
            pt_avg=pt_avg, # number of the point to generate intial condition
            beta_penality=beta_smoothing_ms, # penality for AIC evaluation
            smoothing=smoothing, # the smoothing is done or not?
            type_of_loss=type_of_loss, # type of used loss
            blank_array=zeros(100), # data of all blanks
            display_plot_best_model=false, # one wants the results of the best fit to be plotted
            path_to_plot="NA",
            pt_smooth_derivative=pt_smooth_derivative,
            multiple_scattering_correction=false, # if true uses the given calibration curve to fix the data
            method_multiple_scattering_correction="interpolation",
            calibration_OD_curve="NA", #  the path to calibration curve to fix the data
            verbose=false,
            PopulationSize=PopulationSize,
            maxiters=maxiters,
            abstol=abstol,
            correction_AIC=correction_AIC)

        # selection of te model
        model = model_selection_results[6]
        total_loss = total_loss + model_selection_results[end][end]
        # param of the best model
        temp_res_win = model_selection_results[5]
        param_fitting = copy(temp_res_win)
        #temp_res_win = push!(temp_res_win,model)
        u0 = generating_IC(data_temp, model, smoothing, pt_avg)

        # risimulate data

        remade_solution = ODE_sim_for_iterate(
            model, #string of the model
            u0, # starting condition
            data_temp[1, :], # start time of the sim
            integrator, # which sciml solver of ode
            param_fitting, # parameters of the ODE model
        )
        time_sol = reduce(hcat, remade_solution.t)
        sol_fin = reduce(hcat, remade_solution.u)
        sol_fin = sum(sol_fin, dims=1)
        time_bonduary = remade_solution.t[end]
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
        temp_res_win = vcat(temp_res_win, model_selection_results[4])
        temp_res_win = vcat(label_exp, temp_res_win)
        temp_res_win = vcat(name_well, temp_res_win)







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


    if display_plots
        if_display = display
    else
        if_display = identity
    end

    if save_plot
        mkpath(path_to_plot)
    end

    if_display(
        Plots.scatter(
            data_testing[1, :],
            data_testing[2, :],
            xlabel="Time",
            ylabel="Arb. Units",
            label=["Data " nothing],
            markersize=1,
            color=:black,
            title=string(label_exp, " ", name_well),
        ),
    )
    if_display(
        Plots.vline!(
            interval_changepoints[1:end],
            c=:black,
            label=["change points" nothing],
        ),
    )
    if_display(
        Plots.plot!(
            reduce(vcat, composed_time),
            reduce(vcat, composed_sol),
            xlabel="Time",
            ylabel="Arb. Units",
            label=[" fitting " nothing],
            color=:red,
            title=string(label_exp, " fitting ", name_well),
        ),
    )
    if save_plot
        png(
            string(
                path_to_plot,
                label_exp,
                "_model_selection_seg_",
                n_change_points,
                "_",
                name_well,
                ".png",
            ),
        )
    end

    return param_out, interval_changepoints, composed_time, composed_sol, total_loss
end



"""
    segmentation_ODE(
    data_testing::Matrix{Float64}, 
    name_well::String, 
    label_exp::String, 
    list_of_models::Vector{String}, 
    list_lb_param::Any,
    list_ub_param::Any,
    n_max_change_points::Int;
    detect_number_cpd=true,
    fixed_cpd=false,
    optmizator=BBO_adaptive_de_rand_1_bin_radiuslimited(),
    integrator=Tsit5(),
    type_of_loss="L2",
    type_of_detection="slinding_win",
    type_of_curve="original",
    pt_avg=1, 
    smoothing=true, 
    save_plot=false, 
    display_plot=false,
    path_to_plot="NA", 
    path_to_results="NA",
    win_size=14, 
    pt_smooth_derivative=7,
    penality_parameter=2.0,
    multiple_scattering_correction=false, 
    method_multiple_scattering_correction="interpolation",
    calibration_OD_curve="NA",  
    save_all_model=false,
    method_peaks_detection="peaks_prominence",
    n_bins=40,
    PopulationSize=300,
    maxiters=2000000,
    abstol=0.00001,
    type_of_smoothing="lowess",
    thr_lowess=0.05,
    correction_AIC=true)

This function performs model selection for ordinary differential equation (ODE) models in different segments of the input growth time series data. 
Segmentation is performed with a change points detection algorithm (see (Xx).)

# Arguments:

- `data_testing::Matrix{Float64}`:  The growth curve data. Time values are in the first row and the fit observable (e.g., OD) is in the second row, see documentation.
- `name_well::String`: Name of the well.
- `label_exp::String`: Label of the experiment.
- `list_of_models::Vector{String}`: List of the ODE models of choice.
- `list_lb_param::Any`: Lower bounds for the parameters (compatible with the models).
- `list_ub_param::Any`: Upper bounds for the parameters (compatible with the models).
- `n_max_change_points::Int`: Number of change points of choice, user defined. The results will have different a number of change points depending on the values of the key argument 'type_of_detection' and 'fixed_cpd'.


# Key Arguments:

- `integrator=Tsit5()`: sciML integrator. Use 'KenCarp4(autodiff=true)' to fit piecewise models.
- `optmizator=BBO_adaptive_de_rand_1_bin_radiuslimited()`: Optimizer from optimizationBBO.
- `save_plot=false`: Bool. Options: "true" to save the plot, or "false" not to.
- `display_plots=true`: Bool. Options: "true" to display the plot, or "false" not to.
- `type_of_smoothing="rolling_avg"`: String. Method of choice to smooth the data. Options: "NO", "rolling_avg" (rolling average of the data), and "lowess".
- `pt_avg=7`: Int. Size of the rolling average window smoothing. 
- `pt_smoothing_derivative=7`:Int. Number of points for evaluation of specific growth rate. If <2 it uses interpolation algorithm otherwise a sliding window approach.
- `smoothing=false`: Bool. Options: "true" to smooth the data, or "false" not to.
- `type_of_loss:="RE" `: Type of loss function to be used. Options = "RE" (relative error), "L2" (L2 norm), "L2_derivative" (Xx) and "blank_weighted_L2" (Xx).
- `blank_array=zeros(100)`: Data of all blanks in single array.
- `calibration_OD_curve="NA"`: String, The path where the .csv calibration data are located, used only if `multiple_scattering_correction=true`.
- `multiple_scattering_correction=false`: Bool. Options: "true" to perform the multiple scattering correction (requires a callibration curve) or "false" not to. 
- `method_multiple_scattering_correction="interpolation"`: String. Method of choice to perform the multiple scattering curve inference. Options: '"interpolation"' or '"exp_fit"' (adapted from Meyers, A., Furtmann, C., & Jose, J., *Enzyme and microbial technology*, 118, 1-5., 2018). 
- `thr_lowess=0.05`: Float64. Argument of the lowess smoothing.
- `PopulationSize=100`: Size of the population of the optimization (Xx).
- `maxiters=2000000`: stop criterion, the optimization stops when the number of iterations is bigger than `maxiters`.
- `abstol=0.00001`: stop criterion, the optimization stops when the loss is smaller than `abstol`.
- `beta_penality=2.0` penality  parameters for AIC (or AICc) evaluation.
- 'type_of_detection="slinding_win"': String. Change point detection method of choice. Options `"slinding_win"` (uses a slinding window approach), `"lsdd"` (uses least square density difference (LSDD) from ChangePointDetection.jl). 
- 'type_of_curve="original"': String. Defines the input curve for the change point detection. Options `"original"` for the original time series, and `"deriv"` for performing change point detection on the specific growth rate time series.
- `method_peaks_detection="peaks_prominence"`: How the peak detection is performed on the dissimilarity curve.  `"peaks_prominence"` orders the peaks by prominence. `thr_scan` uses a threshold to choose the peaks
- `n_bins=40`: Int. Used if `method_peaks_detection="thr_scan"`. Number of bins used to generate the threshold that has n_change_points peaks.
- 'detect_number_cpd=true': Bool. Options: true to test all possible combinations of 1, 2, .., n_change_points. The best model is defined on the basis of the AICc criteria. False not to test segment combinations. 
- 'fixed_cpd=false': Bool. Options: true to return the fit using top n_change_points. False not to.
- 'win_size=14': Int. Size of the window used by the cpd algorithms.
- 'path_to_results="NA"':String. Path to save the results. 
- 'save_all_model=false': Bool. Options: true to save all tested models. False not to.

Kimchi uses n_change_points but tests different combinations of the n_change_points+2 top change points if 'detect_number_cpd=false' and 'fixed_cpd=false'.


# Output (if `Model_selection =ODE_Model_selection(...)`:

- `res[1]`. Parameters of each segment.
- `res[2]`. Interval of each ODE segment.
- `res[3]`. Time of the fitted solution.
- `res[4]`. Numerical value of the fitted solution.
- The best fitting model plot if `save_plot_best_model=true` or `display_plot_best_model=true` .


"""
function segmentation_ODE(
    data_testing::Matrix{Float64}, # dataset x times y OD/fluorescence
    name_well::String, # name of the well
    label_exp::String, #label of the experiment
    list_of_models::Vector{String}, # ode model to use
    list_lb_param::Any, # lower bound param
    list_ub_param::Any, # upper bound param
    n_max_change_points::Int;
    detect_number_cpd=true,
    fixed_cpd=false,
    optmizator=BBO_adaptive_de_rand_1_bin_radiuslimited(), # selection of optimization method
    integrator=Tsit5(), # selection of sciml integrator
    type_of_loss="L2", # type of used loss
    type_of_detection="slinding_win",
    type_of_curve="original",
    pt_avg=1, # number of the point to generate intial condition
    smoothing=true, # the smoothing is done or not?
    save_plot=false, # do plots or no
    display_plot=false, # do plots or no
    path_to_plot="NA", # where save plots
    path_to_results="NA",
    win_size=14, #  
    pt_smooth_derivative=7,
    penality_parameter=2.0,
    multiple_scattering_correction=false, # if true uses the given calibration curve to fix the data
    method_multiple_scattering_correction="interpolation",
    calibration_OD_curve="NA",  #  the path to calibration curve to fix the data
    save_all_model=false,
    method_peaks_detection="peaks_prominence",
    n_bins=40,
    PopulationSize=300,
    maxiters=2000000,
    abstol=0.00001,
    type_of_smoothing="lowess",
    thr_lowess=0.05,
    correction_AIC=true)

    # fitting single models
    change_point_list = Vector{Vector{Any}}()
    if n_max_change_points == 0 || detect_number_cpd == true
        res = ODE_Model_selection(
            data_testing, # dataset first row times second row OD
            name_well, # name of the well
            label_exp, #label of the experiment
            list_of_models, # ode model to use
            list_lb_param, # lower bound param
            list_ub_param; # upper bound param
            optmizator=optmizator, # selection of optimization method
            integrator=integrator, # selection of sciml integrator
            pt_avg=pt_avg, # number of the point to generate intial condition
            beta_penality=penality_parameter, # penality for AIC evaluation
            smoothing=smoothing, # the smoothing is done or not?
            type_of_loss=type_of_loss, # type of used loss
            blank_array=zeros(100), # data of all blanks
            display_plot_best_model=false, # one wants the results of the best fit to be plotted
            save_plot_best_model=false,
            path_to_plot="NA",
            pt_smooth_derivative=pt_smooth_derivative,
            multiple_scattering_correction=multiple_scattering_correction, # if true uses the given calibration curve to fix the data
            method_multiple_scattering_correction=method_multiple_scattering_correction,
            calibration_OD_curve=calibration_OD_curve, #  the path to calibration curve to fix the data
            verbose=false,
            PopulationSize=PopulationSize,
            maxiters=maxiters,
            abstol=abstol,
            type_of_smoothing=type_of_smoothing,
            thr_lowess=thr_lowess,
            correction_AIC=correction_AIC)

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
                Tables.table(res[1]),
            )
            CSV.write(
                string(
                    path_to_results,
                    label_exp,
                    "_segmented_ODE_param_",
                    name_well,
                    "_seg_0.csv",
                ),
                Tables.table(Vector(res[2])),
            )

        end

        top_model = res[5]
        score_of_the_models = res[3]
        change_point_list = [0.0]
        change_point_to_plot = [0.0, 0.0, 0.0]
        time_points_to_plot = copy(data_testing[1, :])
        sol_to_plot = copy(reduce(vcat,res[7])[2:2:end] )
    else

        top_model = "NO"
        score_of_the_models = 10^9
        change_point_list = [0.0]
        change_point_to_plot = [0.0, 0.0, 0.0]
        time_points_to_plot = copy(data_testing[1, :])
        sol_to_plot = copy(data_testing[1, :])
    end
    #fitting all model with change points
    if n_max_change_points > 0

        # Fernanda modification to fix non smoothing (to be tested; start)
        # smooting if required
        if smoothing == true
            data_testing = smoothing_data(
                data_testing;
                method=type_of_smoothing,
                pt_avg=pt_avg,
                thr_lowess=thr_lowess
            )
        end
        # Fernanda modification (end)


        if detect_number_cpd == true

            list_change_points_dev = cpd_local_detection(
                data_testing,
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


        elseif fixed_cpd == true

            list_change_points_dev = cpd_local_detection(
                data_testing,
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



        else
            list_change_points_dev = cpd_local_detection(
                data_testing,
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
                list_lb_param, # lower bound param
                list_ub_param, # upper bound param
                cpd_temp;
                type_of_loss=type_of_loss, # type of used loss
                optmizator=optmizator, # selection of optimization method
                integrator=integrator, # selection of sciml integrator
                smoothing=smoothing,
                pt_avg=pt_avg,
                save_plot=false, # do plots or no 
                display_plots=false, # do plots or no
                path_to_plot="NA", # where save plots
                pt_smooth_derivative=pt_smooth_derivative,
                multiple_scattering_correction=multiple_scattering_correction, # if true uses the given calibration curve to fix the data
                method_multiple_scattering_correction=method_multiple_scattering_correction,
                calibration_OD_curve=calibration_OD_curve, #  the path to calibration curve to fix the data
                beta_smoothing_ms=penality_parameter, #  parameter of the AIC penality
                type_of_smoothing=type_of_smoothing,
                thr_lowess=thr_lowess,
            )

            # composing piecewise penality

            n_param =
                sum([
                    length(direct_search_results[1][kk][3:(end-4)]) for
                    kk = 1:length(direct_search_results[1])
                ])

            n_param = n_param + length(cpd_temp)


            new_penality = AICc_evaluation2(n_param, penality_parameter, data_testing[2, :], direct_search_results[end], correction=correction_AIC)


            if new_penality <= score_of_the_models
                score_of_the_models = copy(new_penality)
                top_model = copy(direct_search_results[1])
                time_points_to_plot = copy(direct_search_results[3])
                sol_to_plot = copy(direct_search_results[4])
                change_point_to_plot = copy(direct_search_results[2])
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


    println(size(time_points_to_plot))
    println(size(sol_to_plot))
    if display_plot
        if_display = display
    else
        if_display = identity
    end

    if save_plot == true
        mkpath(path_to_plot)
    end

    if_display(
        Plots.scatter(
            data_testing[1, :],
            data_testing[2, :],
            xlabel="Time",
            ylabel="Arb. Units",
            label=["Data " nothing],
            markersize=1,
            color=:black,
            title=string(label_exp, " ", name_well),
        ),
    )
    if_display(
        Plots.vline!(
            change_point_to_plot[2:end],
            c=:black,
            label=["change points" nothing],
        ),
    )
    if_display(
        Plots.plot!(
            reduce(vcat, time_points_to_plot),
            reduce(vcat, sol_to_plot),
            xlabel="Time",
            ylabel="Arb. Units",
            label=[" fitting " nothing],
            color=:red,
            title=string(label_exp, " fitting ", name_well),
        ),
    )
    if save_plot
        png(
            string(
                path_to_plot,
                label_exp,
                "_model_selection_seg_",
                length(change_point_to_plot[2:end]),
                "_",
                name_well,
                ".png",
            ),
        )
    end

    return top_model, time_points_to_plot, sol_to_plot,score_of_the_models
end

export fitting_one_well_Log_Lin
export fitting_one_well_ODE_constrained
export fitting_one_well_custom_ODE
export ODE_Model_selection
export one_well_morris_sensitivity
export selection_ODE_fixed_intervals
export segmentation_ODE
