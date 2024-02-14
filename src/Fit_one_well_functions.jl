
#######################################################################
"""
fitting single data functions log-lin
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
        r = 1:1:(length(data_smooted[2, :])-pt_smoothing_derivative)
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

        for yy = index_of_max:(length(specific_gr)-1)
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
        index_of_max = findfirst(x -> x == maximum(specific_gr), specific_gr)[1]
        index_gr_max =
            index_of_max +
            findfirst(x -> x < lb_of_distib, specific_gr[index_of_max:end])[1]
        index_gr_min = findlast(x -> x > lb_of_distib, specific_gr[1:index_of_max])[1]
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
    rho = cov(X, Y) / sqrt(var(X) * var(Y))
    d = TDist(N - 2)     # t-Student distribution with N-2 degrees of freedom
    cf = quantile(d, 0.975)  # correction factor for 95% confidence intervals (two-tailed distribution)
    confidence_band = cf * Cantrell_errors * sqrt.(1 / N .+ (data_to_fit_times .- mean(data_to_fit_times)) .^ 2 / var(data_to_fit_times) / (N - 1))
    println(confidence_band)


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

    return res_param, remade_solution.t, sol_fin[1, 1:end]
end

#######################################################################
#######################################################################
"""
fitting custom ODE
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
    pt_smooth_derivative=7,
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
model selection functions
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
    pt_avg=1, # number of the point to generate intial condition
    beta_penality=2.0, # penality for AIC evaluation
    smoothing=false, # the smoothing is done or not?
    type_of_smoothing="lowess",
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
)
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

    for mm = 1:length(models_list)
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

        AICc = AICc_evaluation(param_number, beta_penality, data[2, :], sol_fin)
        
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
    println(param_min)

    tsteps = data[1, :]
    tspan = (data[1, 1], data[1, end])
    u0 = generating_IC(data, model, smoothing, pt_avg)
    ODE_prob = model_selector(model, u0, tspan, param_min)
    sim = solve(ODE_prob, integrator, saveat=tsteps)
    sol_t = reduce(hcat, sim.u)
    sol_time = reduce(hcat, sim.t)
    sol_t = sum(sol_t, dims=1)
    sol_fin, index_not_zero = remove_negative_value(sol_fin)

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
morris sensitivity function
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
    type_of_smoothing="lowess",
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
ODE segementation fitting
"""

function selection_ODE_fixed_change_points(
    data_testing::Matrix{Float64}, # dataset first row times second row OD
    name_well::String, # name of the well
    label_exp::String, #label of the experiment
    list_of_models::Vector{String}, # ode models to use
    list_lb_param::Any, # lower bound param
    list_ub_param::Any, # upper bound param
    n_change_points::Int;
    type_of_loss="L2", # type of used loss
    optmizator=BBO_adaptive_de_rand_1_bin_radiuslimited(), # selection of optimization method
    integrator=Tsit5(), # selection of sciml integrator
    type_of_detection="lsdd",
    type_of_curve="original",
    smoothing=false,
    type_of_smoothing="lowess",
    thr_lowess=0.05,
    pt_avg=1,
    save_plot=false, # do plots or no
    display_plots=false,
    path_to_plot="NA", # where save plots
    win_size=2, # numebr of the point to generate intial condition
    pt_smooth_derivative=7,
    multiple_scattering_correction=false, # if true uses the given calibration curve to fix the data
    method_multiple_scattering_correction="interpolation",
    calibration_OD_curve="NA", #  the path to calibration curve to fix the data
    beta_smoothing_ms=2.0, #  parameter of the AIC penality
    method_peaks_detection="peaks_prominence",
    n_bins=40,
    PopulationSize=300,
    maxiters=2000000,
    abstol=0.000,
)

    if multiple_scattering_correction == true
        data_testing = correction_OD_multiple_scattering(data_testing, calibration_OD_curve; method=method_multiple_scattering_correction)
    end

    if smoothing == true
        data_testing = smoothing_data(
            data;
            method=type_of_smoothing,
            pt_avg=pt_avg,
            thr_lowess=thr_lowess
        )
    end

    list_change_points_dev = cpd_local_detection(
        data_testing,
        n_change_points;
        type_of_detection=type_of_detection,
        type_of_curve=type_of_curve,
        pt_derivative=pt_smooth_derivative,
        size_win=win_size,
        method=method_peaks_detection,
        number_of_bin=n_bins,
    )

    interval_changepoints = push!(list_change_points_dev[2], data_testing[1, 1])
    interval_changepoints = push!(list_change_points_dev[2], data_testing[1, end])
    interval_changepoints = sort(interval_changepoints)
    bc = [data_testing[1, 1], data_testing[2, 2]]
    param_out = Vector{Vector{Any}}()
    composed_sol = Type{Any}
    composed_time = Type{Any}

    for i = 2:(length(interval_changepoints))
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
        )

        # selection of te model
        model = model_selection_results[6]

        # param of the best model
        temp_res_win = model_selection_results[5]
        param_fitting = copy(temp_res_win)
        temp_res_win = vcat(model, temp_res_win)
        temp_res_win = vcat(temp_res_win, model_selection_results[4])

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
        value_bonduary = remade_solution.t[end]
        time_bonduary = sol_fin[end]
        bc = [value_bonduary, time_bonduary]

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
            interval_changepoints[2:end],
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

    return param_out, interval_changepoints, composed_time, composed_sol
end

function ODE_selection_NMAX_change_points(
    data_testing::Matrix{Float64}, # dataset x times y OD/fluorescence
    name_well::String, # name of the well
    label_exp::String, #label of the experiment
    list_lb_param::Any, # lower bound param
    list_ub_param::Any, # upper bound param
    list_of_models::Vector{String}, # ode model to use
    n_max_change_points::Int;
    optmizator=BBO_adaptive_de_rand_1_bin_radiuslimited(), # selection of optimization method
    integrator=Tsit5(), # selection of sciml integrator
    type_of_loss="L2", # type of used loss
    type_of_detection="lsdd",
    type_of_curve="original",
    pt_avg=1, # number of the point to generate intial condition
    smoothing=true, # the smoothing is done or not?
    save_plot=false, # do plots or no
    display_plot=false, # do plots or no
    path_to_plot="NA", # where save plots
    path_to_results="NA",
    win_size=2, # numebr of the point to generate intial condition
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
)

    # fitting single models
    change_point_list = Vector{Vector{Any}}()

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
    )

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
            Tables.table(Matrix(res[3])),
        )
        CSV.write(
            string(
                path_to_results,
                label_exp,
                "_segmented_ODE_param_",
                name_well,
                "_seg_0.csv",
            ),
            Tables.table(Vector(res[5])),
        )
        CSV.write(
            string(
                path_to_results,
                label_exp,
                "_segmented_ODE_solution_",
                name_well,
                "_seg_0.csv",
            ),
            Tables.table(Vector(res[7])),
        )
    end

    top_model = res[5]
    score_of_the_models = res[3]
    change_point_list = [0.0]
    change_point_to_plot = [0.0, 0.0, 0.0]
    time_points_to_plot = copy(data_testing[1, :])
    sol_to_plot = copy(res[7])

    #fitting all model with change points
    if n_max_change_points > 0
        for n = 1:n_max_change_points
            direct_search_results = selection_ODE_fixed_change_points(
                data_testing, # dataset first row times second row OD
                name_well, # name of the well
                label_exp, #label of the experiment
                list_of_models, # ode models to use
                list_lb_param, # lower bound param
                list_ub_param, # upper bound param
                n;
                type_of_loss=type_of_loss, # type of used loss
                optmizator=optmizator, # selection of optimization method
                integrator=integrator, # selection of sciml integrator
                type_of_detection=type_of_detection,
                type_of_curve=type_of_curve,
                smoothing=smoothing,
                pt_avg=pt_avg,
                save_plot=false, # do plots or no
                display_plots=false, # do plots or no
                path_to_plot="NA", # where save plots
                win_size=win_size, # numebr of the point to generate intial condition
                pt_smooth_derivative=pt_smooth_derivative,
                multiple_scattering_correction=multiple_scattering_correction, # if true uses the given calibration curve to fix the data
                method_multiple_scattering_correction=method_multiple_scattering_correction,
                calibration_OD_curve=calibration_OD_curve, #  the path to calibration curve to fix the data
                beta_smoothing_ms=penality_parameter, #  parameter of the AIC penality
                method_peaks_detection=method_peaks_detection,
                n_bins=n_bins,
                type_of_smoothing=type_of_smoothing,
                thr_lowess=thr_lowess,
            )

            # composing piecewise penality

            n_param =
                sum([
                    length(direct_search_results[1][kk][3:(end-5)]) for
                    kk = 1:length(direct_search_results[1])
                ]) + n_change_points

            new_penality = AICc_evaluation(n_param, penality_parameter, data, unique(direct_search_results[end]))


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
                        name_well,
                        "_seg_",
                        n,
                        ".csv",
                    ),
                    Tables.table(Vector(direct_search_results[1])),
                )
                CSV.write(
                    string(
                        path_to_results,
                        label_exp,
                        "_segmented_ODE_solution_",
                        name_well,
                        "_seg_",
                        n,
                        ".csv",
                    ),
                    Tables.table(Vector(direct_search_results[3])),
                )
            end
        end
    end

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

    return top_model, time_points_to_plot, sol_to_plot
end