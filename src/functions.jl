"""
 Declaring the ODE models
"""

#######################################################################
#

include("models.jl");
include("Fit_one_file_functions.jl");
include("pre_processing_functions.jl");

"""
Internal functions
"""

function model_selector(model::String, u0, tspan, param = nothing)
    """
      generate sciML ODE problem
    """

    ODE_prob = ODEProblem(models[model].func, u0, tspan, param)

    return ODE_prob
end

function specific_gr_evaluation(data_smooted::Any, pt_smoothing_derivative::Int)
    """
    specific gr evaluation with slinding window log-lin fitting
    data_smooted = matrix of data
    pt_smoothing_derivative = size of the win, if <2 the the numerical derivative of (log) data is evaluate with interpolation algorithm
    """

    if pt_smoothing_derivative > 1
        specific_gr = [
            CurveFit.curve_fit(
                LinearFit,
                data_smooted[1, r:(r+pt_smoothing_derivative)],
                log.(data_smooted[2, r:(r+pt_smoothing_derivative)]),
            ).coefs[2] for r = 1:1:(length(data_smooted[2, :])-pt_smoothing_derivative)
        ]
    else
        itp =
            interpolate((data_smooted[1, :],), log.(data_smooted[2, :]), Gridded(Linear()))
        specific_gr = only.(Interpolations.gradient.(Ref(itp), data_smooted[1, :]))
    end

    return specific_gr
end

function specific_gr_interpol_evaluation(data_testing)
    itp = interpolate((data_testing[1, :],), data_testing[2, :], Gridded(Linear()))
    specific_gr_interpol = only.(Interpolations.gradient.(Ref(itp), data_testing[1, :]))

    return specific_gr_interpol
end


function vectorize_df_results(
    well_name::String,
    model::String,
    res::Any,
    th_gr::Any,
    em_gr::Any,
    loss::Float64,
)
    """
        internal function to reder as vector the results of the hardcoded ODE
    """

    res_param = [
        [model, well_name]
        res[1:length(models[model].params)-5]
        [th_gr, em_gr, loss]
    ]

    return res_param
end

function initialize_df_results(model::String)
    """
        internal function to name the vectors the results of the hardcoded ODE
    """

    param_names = models[model].params

    return param_names
end

function guess_param(lb_param::Vector{Float64}, ub_param::Vector{Float64})
    """
    internal function to set the start of the optimization problem in the middle of the  box constrains
    """
    param = lb_param .+ (ub_param - lb_param) ./ 2

    return param
end





function generation_of_combination_of_IC_morris(
    lb_param::Vector{Float64},
    ub_param::Vector{Float64},
    N_step_morris::Int,
)

    starting_guess = copy(lb_param)
    delta_vec =
        [(ub_param[i] - lb_param[i]) / (N_step_morris + 1) for i = 1:length(ub_param)]

    # generating combinations of all possible parameters values
    combinations_par = copy(starting_guess)
    combinations_tot = copy(combinations_par)

    for k = 1:N_step_morris
        for ll = 1:(size(delta_vec)[1])
            combinations_par[ll] = combinations_par[ll] + delta_vec[ll]
            combinations_tot = hcat(combinations_tot, combinations_par)
        end
    end

    return combinations_tot
end

function generating_IC(data::Matrix{Float64}, model::String, smoothing::Bool, pt_avg::Int)
    # "starting condition  using data if smoothing average is used skip this part

    if smoothing == true
        u0 = [data[2, 1]]
    else
        u0 = [Statistics.mean(data[2, 1:pt_avg])]
    end

    if model == "HPM" ||
       model == "aHPM" ||
       model == "HPM_inhibition" ||
       model == "aHPM_inhibition" ||
       model == "ODEs_HPM_SR" ||
       model == "HPM_exp"
        # specific initial condition for system of ODEs
        # all the biomass starts as dormient
        if smoothing == true
            u0 = [data[2, 1], 0.0]
        else
            u0 = [Statistics.mean(data[2, 1:pt_avg]), 0.0]
        end
    end
    # "starting condition  using data if the model is an 3 state HPM the first equation is set equal to the starting pop and the other to zero

    if model == "HPM_3_death" ||
       model == "HPM_3_inhibition" ||
       model == "HPM_3_death_resistance" ||
       model == "aHPM_3_death_resistance"
        if smoothing == true
            u0 = [data[2, 1], 0.0, 0.0]
        else
            u0 = [Statistics.mean(data[2, 1:pt_avg]), 0.0, 0.0]
        end
    end

    return u0
end

function generating_IC_custom_ODE(
    data::Matrix{Float64},
    n_equation::Int,
    smoothing::Bool,
    pt_avg::Int,
)

    u0 = zeros(n_equation)
    # "starting condition  using data if smoothing average is used skip this part

    if smoothing == true
        u0[1] = data[2, 1]
    else
        u0[1] = Statistics.mean(data[2, 1:pt_avg])
    end

    return u0
end

function correction_OD_multiple_scattering(data::Matrix{Float64}, calibration_curve::String)

    od_calib = CSV.File(calibration_curve)
    names_of_cols = propertynames(od_calib)
    Od_real = od_calib[names_of_cols[1]]
    od_calib_array = Matrix(transpose(hcat(Od_real, od_calib[names_of_cols[2]])))
    soterd_calib = sort!(od_calib_array, rev = false, dims = 2)
    itp =
        interpolate(soterd_calib[1, :], soterd_calib[2, :], SteffenMonotonicInterpolation())
    extrap_spline = extrapolate(itp, 0)
    corrected_data = [extrap_spline(k) for k in data[2, :]]
    data_fin = Matrix(transpose(hcat(data[1, :], corrected_data)))

    return data_fin
end


#######################################################################
"""
sim ODE functions
"""

function ODE_sim(
    model::String, #string of the model
    n_start::Vector{Float64}, # starting condition
    tstart::Float64, # start time of the sim
    tmax::Float64, # final time of the sim
    delta_t::Float64, # delta t for poisson approx
    param_of_ode::Vector{Float64}; # parameters of the ODE model
    integrator = KenCarp4(), # which sciml solver of ode
)

    # defining time stepping
    t_steps = tstart:delta_t:tmax
    tspan = (tstart, tmax)
    u0 = n_start
    ODE_prob = model_selector(model, u0, tspan, param_of_ode)
    sim = solve(ODE_prob, integrator, saveat = t_steps)

    return sim
end

function ODE_sim_for_iterate(
    model::String, #string of the model
    n_start::Vector{Float64}, # starting condition
    array_time::Vector{Float64},
    integrator::Any, # which sciml solver of ode
    param_of_ode::Any, # parameters of the ODE model
)

    # defining time stepping
    t_steps = array_time
    tspan = (array_time[1], array_time[end])
    u0 = n_start
    ODE_prob = model_selector(model, u0, tspan, param_of_ode)
    sim = solve(ODE_prob, integrator, saveat = t_steps)

    return sim
end

#######################################################################
"""
sim stochastic function
"""

function stochastic_sim(
    model::String, #string of the model
    n_start::Int, # number of starting cells
    n_mol_start::Float64, # starting concentration of the limiting nutrients
    tstart::Float64, # start time of the sim
    tmax::Float64, # final time of the sim
    delta_t::Float64, # delta t for poisson approx
    k_1_val::Float64,
    k_2_val::Float64, # monod constant
    alpha_val::Float64, # massimum possible growth rate
    lambda::Float64, # lag time
    n_mol_per_birth::Float64,# nutrient consumed per division (conc)
    volume::Float64,
)

    #inizialization of times
    tot_pop = [copy(n_start)]
    times = [copy(tstart)]
    conc_of_nutriens = [copy(n_mol_start / volume)]
    n_times = floor((tmax - tstart) / delta_t)

    for i = 2:n_times
        # if for defining the lag phase

        if i * delta_t < lambda
            rate_per_cell = 0.0
        else

            if model == "Monod"
                factor = conc_of_nutriens[end] / (k_1_val + conc_of_nutriens[end])
            end

            if model == "Haldane"
                factor =
                    conc_of_nutriens[end] /
                    (k_1_val + conc_of_nutriens[end] + (conc_of_nutriens[end])^2 / k_2_val)
            end

            if model == "Blackman"
                factor = conc_of_nutriens[end] / (k_1_val)
            end

            if model == "Tesseir"
                factor = 1 - exp(conc_of_nutriens[end] * k_1_val)
            end

            if model == "Moser"
                factor =
                    conc_of_nutriens[end]^k_2_val /
                    (conc_of_nutriens[end]^k_2_val + k_1_val)
            end

            if model == "Aiba-Edwards"
                factor =
                    exp(-k_2_val / conc_of_nutriens[end]) * conc_of_nutriens[end] /
                    (conc_of_nutriens[end] + k_1_val)
            end

            if model == "Verhulst"
                factor = 1 - tot_pop[end] / k_1_val
            end

            rate_per_cell = alpha_val * factor
        end

        # evaluating the number of birth events with poisson approx
        total_birth_rate = rate_per_cell .* tot_pop[end] .* delta_t
        n_birth = rand(Poisson(total_birth_rate), 1)

        # updating states
        new_conc = max(0, conc_of_nutriens[end] - n_birth[1] * n_mol_per_birth / volume)
        net_pop_variation = n_birth[1]

        new_pop = max(0, tot_pop[end] + net_pop_variation)
        conc_of_nutriens = push!(conc_of_nutriens, new_conc)
        tot_pop = push!(tot_pop, new_pop)
        times = push!(times, times[end] + delta_t)

    end

    return tot_pop, conc_of_nutriens, times
end

#######################################################################

"""
plotting functions
"""
function plot_data(
    label_exp::String, #label of the experiment
    path_to_data::String, # path to the folder to analyze
    path_to_annotation::String;# path to the annotation of the wells
    path_to_plot = "NA", # path where to save Plots
    display_plots = true,# display plots in julia or not
    save_plot = false, # save the plot or not
    overlay_plots = true, # true a single plot for all dataset false one plot per well
    blank_subtraction = "NO", # string on how to use blank (NO,avg_subtraction,time_avg)
    average_replicate = false, # if true the average between replicates
    correct_negative = "thr_correction", # if "thr_correction" it put a thr on the minimum value of the data with blank subracted, if "blank_correction" uses blank distrib to impute negative values
    thr_negative = 0.01,  # used only if correct_negative == "thr_correction"
)

    """
    function that plot a dataset
    """
    annotation = CSV.File(string(path_to_annotation), header = false)
    names_of_annotated_df = [annotation[l][1] for l = 1:length(annotation)]
    # selcting blank wells
    properties_of_annotation = [annotation[l][2] for l = 1:length(annotation)]
    list_of_blank = names_of_annotated_df[findall(x -> x == "b", properties_of_annotation)]
    list_of_discarded =
        names_of_annotated_df[findall(x -> x == "X", properties_of_annotation)]
    list_of_blank = Symbol.(list_of_blank)

    # reading files
    dfs_data = CSV.File(path_to_data)

    # shaping df for the inference
    names_of_cols = propertynames(dfs_data)

    # excluding blank data and discarded wells
    names_of_cols = filter!(e -> !(e in list_of_blank), names_of_cols)
    if length(list_of_discarded) > 0
        names_of_cols = filter!(e -> !(e in list_of_discarded), names_of_cols)
    end

    times_data = dfs_data[names_of_cols[1]]
    if length(list_of_blank) > 0
        blank_array = reduce(vcat, [(dfs_data[k]) for k in list_of_blank])
        blank_array = convert(Vector{Float64}, blank_array)
    else
        blank_array = 0.0
    end

    ## BLANK ANALYSIS HERE
    if blank_subtraction == "avg_blank"
        blank_value = mean([mean(dfs_data[k]) for k in list_of_blank])
    elseif blank_subtraction == "time_blank"
        blank_value =
            [mean([dfs_data[k][j] for k in list_of_blank]) for j = 1:length(times_data)]
    else
        blank_value = zeros(length(times_data))
    end

    ## considering replicates
    list_replicate = unique(properties_of_annotation)
    list_replicate = filter!(e -> e != "b", list_replicate)

    if average_replicate == true
        new_data = times_data
        list_replicate = unique(properties_of_annotation)
        list_replicate = filter!(e -> e != "b", list_replicate)

        for replicate_temp in list_replicate
            names_of_replicate_temp =
                Symbol.(
                    names_of_annotated_df[findall(
                        x -> x == replicate_temp,
                        properties_of_annotation,
                    )]
                )
            replicate_mean = [
                mean([dfs_data[k][j] for k in names_of_replicate_temp]) for
                j = 1:length(times_data)
            ]
            new_data = hcat(new_data, replicate_mean)
        end

        new_data = DataFrame(new_data, :auto)
        rename!(new_data, vcat(:Time, reduce(vcat, Symbol.(list_replicate))))
        names_of_cols = propertynames(new_data)
        dfs_data = new_data
    end

    # creating the folder of data if one wants to save
    if save_plot == true
        mkpath(path_to_plot)
    end

    for well_name in names_of_cols[2:end]
        name_well = string(well_name)

        if average_replicate == true
            data_values = copy(dfs_data[!, well_name])
        else
            data_values = copy(dfs_data[well_name])
        end

        # blank subtraction
        data_values = data_values .- blank_value
        data = Matrix(transpose(hcat(times_data, data_values)))

        if correct_negative == "thr_correction"
            data = thr_negative_correction(data, thr_negative)
        end

        if correct_negative == "blank_correction"
            data = blank_distrib_negative_correction(data, blank_array)
        end

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
                        xlabel = "Time",
                        ylabel = "Arb. Units",
                        label = [name_well],
                        title = string(label_exp),
                        legend = :outertopright,
                    ),
                )
            else
                if_display(
                    Plots.plot!(
                        data[1, :],
                        data[2, :],
                        xlabel = "Time",
                        ylabel = "Arb. Units",
                        label = [name_well],
                        title = string(label_exp),
                        legend = :outertopright,
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
                    xlabel = "Time",
                    ylabel = "Arb. Units",
                    label = ["Data " nothing],
                    color = :black,
                    title = string(label_exp, " ", name_well),
                ),
            )
            if save_plots
                png(string(path_to_plot, label_exp, "_", name_well, ".png"))
            end
        end

    end
end

#######################################################################
"""
fitting single data functions log-lin
"""

function fitting_one_well_Log_Lin(
    data::Matrix{Float64}, # dataset first row times second row OD
    name_well::String, # name of the well
    label_exp::String; #label of the experiment
    display_plots = false, # do plots or no
    save_plot = false,
    path_to_plot = "NA", # where save plots
    type_of_smoothing = "rolling_avg", # option, NO, gaussian, rolling avg
    pt_avg = 7, # number of the point for rolling avg not used in the other cases
    pt_smoothing_derivative = 7, # number of poits to smooth the derivative
    pt_min_size_of_win = 7, # minimum size of the exp windows in number of smooted points
    type_of_win = "maximum", # how the exp. phase win is selected, "maximum" of "global_thr"
    threshold_of_exp = 0.9, # threshold of growth rate in quantile to define the exp windows
    multiple_scattering_correction = false, # if true uses the given calibration curve to fix the data
    calibration_OD_curve = "NA", #  the path to calibration curve to fix the data
    thr_lowess = 0.05, # keyword argument of lowees smoothing
)

    if multiple_scattering_correction == true
        data = correction_OD_multiple_scattering(data, calibration_OD_curve)
    end

    if type_of_smoothing == "rolling_avg"
        data_smooted = smoothing_data(data, pt_avg)
    elseif type_of_smoothing == "lowess"
        # lowess call here
        model_fit = lowess_model(data[1, :], data[2, :], thr_lowess)
        data_smooted = Matrix(transpose(hcat(data[1, :], model_fit)))
    else
        data_smooted = data
    end

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
    fitting_results = CurveFit.curve_fit(LinearFit, data_to_fit_times, data_to_fit_values)

    # residual calculation
    residual = [
        (
            (
                data_to_fit_values[ll] - fitting_results.coefs[2] * data_to_fit_times[ll] -
                fitting_results.coefs[1]
            )^2
        ) for ll = 1:length(data_to_fit_values)
    ]
    # sum_of_squares calculation
    sum_of_squares = [
        ((data_to_fit_values[ll] - mean(data_to_fit_values))^2) for
        ll = 1:length(data_to_fit_values)
    ]
    # coeff of determination
    rsquared = 1 - sum(residual) / sum(sum_of_squares)
    # confidence interval growth rate
    a = ((1 / (length(residual) - 2)) * sum(residual))
    std_error = sqrt(a / sum(sum_of_squares))
    p95 = ccdf(TDist(length(sum_of_squares) - 2), (1 - 0.05 / 2))
    confidence_coeff_2 = std_error * p95
    # confidence interval intercept
    a_1 = sqrt((1 / (length(residual))) * sum(data_to_fit_values .^ 2))
    std_error_intercept = std_error * a_1
    confidence_coeff_1 = std_error_intercept * p95

    ###
    # EVALUATING CONFIDENCE BANDS
    term_1 = ((1 / (length(residual) - 2)) * sum(residual))
    term_2 =
        (1 / (length(residual))) .+
        (data_to_fit_times .- mean(data_to_fit_values)) ./ (sum(sum_of_squares))
    confidence_band = p95 .* sqrt.(term_1 .* term_2)
    fitted_line = [
        fitting_results.coefs[2] * data_to_fit_times[ll] + fitting_results.coefs[1] for
        ll = 1:length(data_to_fit_times)
    ]
    # storing results

    results_lin_log_fit = [
        label_exp,
        name_well,
        data_to_fit_times[1],
        data_to_fit_times[end],
        specific_gr_times[index_of_max],
        gr_max,
        fitting_results.coefs[2],
        confidence_coeff_2,
        log(2) / (fitting_results.coefs[2]),
        log(2) / (fitting_results.coefs[2] - confidence_coeff_2),
        log(2) / (fitting_results.coefs[2] + confidence_coeff_2),
        fitting_results.coefs[1],
        confidence_coeff_1,
        rsquared,
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
            xlabel = "Time",
            ylabel = "Log(Arb. Units)",
            label = ["Data " nothing],
            markersize = 1,
            color = :black,
            title = string(label_exp, " ", name_well),
        ),
    )
    if_display(
        Plots.plot!(
            data_to_fit_times,
            fitted_line,
            ribbon = confidence_band,
            xlabel = "Time ",
            ylabel = "Log(Arb. Units)",
            label = [string("Fitting Log-Lin ") nothing],
            c = :red,
        ),
    )
    if_display(
        Plots.vline!(
            [data_to_fit_times[1], data_to_fit_times[end]],
            c = :black,
            label = [string("Window of exp. phase ") nothing],
        ),
    )
    if save_plot
        png(string(path_to_plot, label_exp, "_Log_Lin_Fit_", name_well, ".png"))
    end
    if_display(
        Plots.scatter(
            specific_gr_times,
            specific_gr,
            xlabel = "Time ",
            ylabel = "1 /time ",
            label = [string("Dynamics growth rate ") nothing],
            c = :red,
        ),
    )
    if_display(
        Plots.vline!(
            [data_to_fit_times[1], data_to_fit_times[end]],
            c = :black,
            label = [string("Window of exp. phase ") nothing],
        ),
    )
    if save_plot
        png(string(path_to_plot, label_exp, "_dynamics_gr_", name_well, ".png"))
    end

    return results_lin_log_fit
end

#######################################################################




function loss_L2_derivative(data, ODE_prob, integrator, p, tsteps)
    sol = solve(
        ODE_prob,
        integrator,
        p = p,
        saveat = tsteps,
        verbose = false,
        abstol = 1e-10,
        reltol = 1e-10,
    )
    sol_t = reduce(hcat, sol.u)
    sol_time = reduce(hcat, sol.t)
    sol_t = reduce(hcat, sol.u)
    sol_t = sum(sol_t, dims = 1)

    itp = interpolate((sol_time,), sol_t, Gridded(Linear()))
    derivative_theo = only.(Interpolations.gradient.(Ref(itp), sol_t))

    itp = interpolate((data[1, :],), data[2, :], Gridded(Linear()))
    derivative_data = only.(Interpolations.gradient.(Ref(itp), data[2, :]))

    if size(derivative_theo)[1] == size(derivative_data)[1]
        lossa =
            NaNMath.sum(abs2.((derivative_theo[2, :] .- derivative_data[1, 1:end]))) /
            length(derivative_data[2, :])
    else
        lossa = 10.0^9 * length(data[2, :])
    end

    return lossa, sol
end

function loss_blank_weighted_L2(data, ODE_prob, integrator, p, tsteps, blank_array)
    sol = solve(
        ODE_prob,
        integrator,
        p = p,
        saveat = tsteps,
        verbose = false,
        abstol = 1e-10,
        reltol = 1e-10,
    )
    sol_t = reduce(hcat, sol.u)
    sol_t = sum(sol_t, dims = 1)

    empirical_blank_distrib = (blank_array .- mean(blank_array))
    # generation of the empirica distrib respect the mean of the noise
    test = StatsBase.fit(Histogram, empirical_blank_distrib)
    binning_distrib = test.edges
    probabiltity_distrib = test.weights ./ sum(test.weights)

    if size(sol_t)[2] == size(data)[2]
        lossa = 0.0

        for ll = 1:size(sol_t)[2]
            dist = (data[2, ll] - sol_t[1, ll])
            index = findfirst(x -> (x > dist), binning_distrib[1])
            if (typeof(index) == Nothing || index > length(probabiltity_distrib))
                lossa = lossa + abs2.(dist) / length(data[2, :])
            else
                prob = probabiltity_distrib[index]
                lossa = lossa + abs2.((1 - prob) * dist) / length(data[2, :])
            end
        end
    else
        lossa = 10.0^9 * length(data[2, :])
    end

    return lossa, sol
end

function loss_L2(data, ODE_prob, integrator, p, tsteps)
    sol = solve(
        ODE_prob,
        integrator,
        p = p,
        saveat = tsteps,
        verbose = false,
        abstol = 1e-10,
        reltol = 1e-10,
    )
    sol_t = reduce(hcat, sol.u)
    sol_t = sum(sol_t, dims = 1)

    if size(sol_t)[2] == size(data)[2]
        lossa = NaNMath.sum(abs2.((data[2, :] .- sol_t[1, 1:end]))) / length(data[2, :])
    else
        lossa = 10.0^9 * length(data[2, :])
    end

    return lossa, sol
end

function loss_RE(data, ODE_prob, integrator, p, tsteps)
    sol = solve(
        ODE_prob,
        integrator,
        p = p,
        saveat = tsteps,
        verbose = false,
        abstol = 1e-10,
        reltol = 1e-10,
    )
    sol_t = reduce(hcat, sol.u)
    sol_t = sum(sol_t, dims = 1)

    if size(sol_t)[2] == size(data)[2]
        lossa =
            0.5 * NaNMath.sum(abs2.(1.0 .- (data[2, :] ./ sol_t[1, 1:end]))) /
            length(data[2, :])
    else
        lossa = 10.0^9 * length(data[2, :])
    end

    return lossa, sol
end

function select_loss_function(loss_name, data, ODE_prob, integrator, tsteps, blank_array)
    loss_functions = Dict(
        "L2" => loss_L2,
        "RE" => loss_RE,
        "blank_weighted_L2" => loss_blank_weighted_L2,
        "L2_derivative" => loss_L2_derivative,
    )

    if loss_name == "blank_weighted_L2"
        return (p) ->
            loss_functions[loss_name](data, ODE_prob, integrator, p, tsteps, blank_array)
    else
        return (p) -> loss_functions[loss_name](data, ODE_prob, integrator, p, tsteps)
    end
end

#######################################################################

"""
fitting single data function ODE
"""

function fitting_one_well_ODE_constrained(
    data::Matrix{Float64}, # dataset first row times second row OD
    name_well::String, # name of the well
    label_exp::String, #label of the experiment
    model::String, # ode model to use
    lb_param::Vector{Float64}, # lower bound param
    ub_param::Vector{Float64}; # upper bound param
    param = lb_param .+ (ub_param .- lb_param) ./ 2,# initial guess param
    optmizator = BBO_adaptive_de_rand_1_bin_radiuslimited(), # selection of optimization method
    integrator = KenCarp4(autodiff = true), # selection of sciml integrator
    display_plots = true, # display plots in julia or not
    save_plot = false,
    path_to_plot = "NA", # where save plots
    pt_avg = 1, # numebr of the point to generate intial condition
    pt_smooth_derivative = 7,
    smoothing = false, # the smoothing is done or not?
    type_of_smoothing = "rolling_avg",
    type_of_loss = "RE", # type of used loss
    blank_array = zeros(100), # data of all blanks
    multiple_scattering_correction = false, # if true uses the given calibration curve to fix the data
    calibration_OD_curve = "NA",  #  the path to calibration curve to fix the data
    PopulationSize = 300,
    maxiters = 2000000,
    abstol = 0.00001,
    thr_lowess = 0.05,
)

    if multiple_scattering_correction == true
        data = correction_OD_multiple_scattering(data, calibration_OD_curve)
    end

    #defining time interval

    max_t = data[1, end]
    min_t = data[1, 1]
    tspan = (min_t, max_t)
    tsteps = data[1, :]

    # smoothing data if required
    if type_of_smoothing == "rolling_avg"
        data = smoothing_data(data, pt_avg)
    elseif type_of_smoothing == "lowess"
        # lowess call here
        model_fit = lowess_model(data[1, :], data[2, :], thr_lowess)
        data = Matrix(transpose(hcat(data[1, :], model_fit)))
    else
        data = data
    end

    # setting initial conditions
    u0 = generating_IC(data, model, smoothing, pt_avg)

    #  definition  ode symbolic problem
    ODE_prob = model_selector(model, u0, tspan)

    ## defining loss function
    loss_function = select_loss_function(type_of_loss, data, ODE_prob, integrator, tsteps, blank_array)
    optf = Optimization.OptimizationFunction((x, p) -> loss_function(x))
    optprob_const =
        Optimization.OptimizationProblem(optf, param, u0, lb = lb_param, ub = ub_param)
    res = Optimization.solve(
        optprob_const,
        optmizator,
        PopulationSize = PopulationSize,
        maxiters = maxiters,
        abstol = abstol,
    )

    #revalution of solution for plot an loss evaluation
    remade_solution = solve(remake(ODE_prob, p = res.u), integrator, saveat = tsteps)
    sol_time = reduce(hcat, remade_solution.t)
    sol_fin = reduce(hcat, remade_solution.u)
    sol_fin = sum(sol_fin, dims = 1)

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
            xlabel = "Time",
            ylabel = "Arb. Units",
            label = ["Data " nothing],
            markersize = 2,
            color = :black,
            title = string(label_exp, " ", name_well),
        ),
    )
    if_display(
        Plots.plot!(
            remade_solution.t,
            sol_fin[1, 1:end],
            xlabel = "Time",
            ylabel = "Arb. Units",
            label = [string("Fitting ", model) nothing],
            c = :red,
        ),
    )
    if save_plot
        png(string(path_to_plot, label_exp, "_", model, "_", name_well, ".png"))
    end

    # here problem
    # max_theoretical gr
    data_th = vcat(sol_time, sol_fin)
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
    param = lb_param .+ (ub_param .- lb_param) ./ 2,# initial guess param
    optmizator = BBO_adaptive_de_rand_1_bin_radiuslimited(), # selection of optimization method
    integrator = KenCarp4(autodiff = true), # selection of sciml integrator
    display_plots = false, # do plots or no
    save_plot = false,
    path_to_plot = "NA", # where save plots
    pt_avg = 1, # numebr of the point to generate intial condition
    pt_smooth_derivative = 7,
    smoothing = false, # the smoothing is done or not?
    type_of_loss = "RE", # type of used loss
    blank_array = zeros(100), # data of all blanks
    multiple_scattering_correction = false, # if true uses the given calibration curve to fix the data
    calibration_OD_curve = "NA",  #  the path to calibration curve to fix the data
    PopulationSize = 300,
    maxiters = 2000000,
    abstol = 0.00001,
    thr_lowess = 0.05,
    type_of_smoothing = "lowess",
)

    if multiple_scattering_correction == true
        data = correction_OD_multiple_scattering(data, calibration_OD_curve)
    end

    #defining time interval
    max_t = data[1, end]
    min_t = data[1, 1]
    tspan = (min_t, max_t)
    tsteps = data[1, :]

    # smoothing data if required
    if type_of_smoothing == "rolling_avg" && smoothing == true
        data = smoothing_data(data, pt_avg)
    elseif type_of_smoothing == "lowess" && smoothing == true
        # lowess call here
        model_fit = lowess_model(data[1, :], data[2, :], thr_lowess)
        data = Matrix(transpose(hcat(data[1, :], model_fit)))
    else
        data = copy(data)
    end

    u0 = generating_IC_custom_ODE(data, n_equation, smoothing, pt_avg)

    #  definition  ode symbolic problem
    ODE_prob = ODEProblem(model, u0, tspan, nothing)

    loss_function =
        select_loss_function(type_of_loss, data, ODE_prob, integrator, tsteps, blank_array)
    optf = Optimization.OptimizationFunction((x, p) -> loss_function(x))
    optprob_const =
        Optimization.OptimizationProblem(optf, param, u0, lb = lb_param, ub = ub_param)
    res = Optimization.solve(
        optprob_const,
        optmizator,
        PopulationSize = PopulationSize,
        maxiters = maxiters,
        abstol = abstol,
    )

    #revalution of solution for plot an loss evaluation
    remade_solution = solve(remake(ODE_prob, p = res.u), integrator, saveat = tsteps)
    sol_time = reduce(hcat, remade_solution.t)
    sol_fin = reduce(hcat, remade_solution.u)
    sol_fin = sum(sol_fin, dims = 1)

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
            xlabel = "Time",
            ylabel = "Arb. Units",
            label = ["Data " nothing],
            markersize = 2,
            color = :black,
            title = string(label_exp, " ", name_well),
        ),
    )
    if_display(
        Plots.plot!(
            remade_solution.t,
            sol_fin[1, 1:end],
            xlabel = "Time",
            ylabel = "Arb. Units",
            label = [string("Fitting custom model") nothing],
            c = :red,
        ),
    )
    if save_plot
        png(string(path_to_plot, label_exp, "_custom_model_", name_well, ".png"))
    end

    #max_theoretical gr
    data_th = vcat(sol_time, sol_fin)
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
    optmizator = BBO_adaptive_de_rand_1_bin_radiuslimited(), # selection of optimization method
    integrator = KenCarp4(autodiff = true), # selection of sciml integrator
    pt_avg = 1, # number of the point to generate intial condition
    beta_penality = 2.0, # penality for AIC evaluation
    smoothing = false, # the smoothing is done or not?
    type_of_smoothing = "lowess",
    thr_lowess = 0.05,
    type_of_loss = "L2", # type of used loss
    blank_array = zeros(100), # data of all blanks
    display_plot_best_model = false, # one wants the results of the best fit to be plotted
    save_plot_best_model = false,
    path_to_plot = "NA",
    pt_smooth_derivative = 7,
    multiple_scattering_correction = false, # if true uses the given calibration curve to fix the data
    calibration_OD_curve = "NA", #  the path to calibration curve to fix the data
    verbose = false,
    PopulationSize = 300,
    maxiters = 2000000,
    abstol = 0.00001,
)

    if multiple_scattering_correction == true
        data = correction_OD_multiple_scattering(data, calibration_OD_curve)
    end

    # smooting if required
    if type_of_smoothing == "rolling_avg" && smoothing == true
        data = smoothing_data(data, pt_avg)
    elseif type_of_smoothing == "lowess" && smoothing == true
        # lowess call here
        model_fit = lowess_model(data[1, :], data[2, :], thr_lowess)
        data = Matrix(transpose(hcat(data[1, :], model_fit)))
    else
        data = copy(data)
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
            lb = temp_param_lb,
            ub = temp_param_ub,
        )
        res = Optimization.solve(
            optprob_const,
            optmizator,
            PopulationSize = PopulationSize,
            maxiters = maxiters,
            abstol = abstol,
        )

        #revalution of solution for plot an loss evaluation
        remade_solution = solve(remake(ODE_prob, p = res.u), integrator, saveat = tsteps)
        sol_time = reduce(hcat, remade_solution.t)
        sol_fin = reduce(hcat, remade_solution.u)
        sol_fin = sum(sol_fin, dims = 1)
        # here problem
        #max_theoretical gr
        data_th = vcat(sol_time, sol_fin)
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

        df_res_optimization[mm] = res_param
        data_size = size(data)[2]
        param_number = length(temp_start_param)
        results_to_be_pushed = [
            models_list[mm],
            res.objective,
            data_size * log(res.objective / data_size) + beta_penality * param_number,
        ]
        rss_array = hcat(rss_array, results_to_be_pushed)
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
    sim = solve(ODE_prob, integrator, saveat = tsteps)
    sol_t = reduce(hcat, sim.u)
    sol_time = reduce(hcat, sim.t)
    sol_t = sum(sol_t, dims = 1)
    data_th = vcat(sol_time, sol_t)

    if display_plot_best_model
        if_display = display
    else
        if_display = identity
    end

    if save_plot_best_model
        mkpath(path_to_plot)
    end

    data_th = vcat(sol_time, sol_t)
    if_display(
        Plots.scatter(
            data[1, :],
            data[2, :],
            xlabel = "Time",
            ylabel = "Arb. Units",
            label = ["Data " nothing],
            markersize = 2,
            color = :black,
            title = string(label_exp, " ", name_well),
        ),
    )
    if_display(
        Plots.plot!(
            data_th[1, :],
            data_th[2, :],
            xlabel = "Time",
            ylabel = "Arb. Units",
            label = [string("Fitting ", model) nothing],
            c = :red,
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
    N_step_morris = 7,
    optmizator = BBO_adaptive_de_rand_1_bin_radiuslimited(), # selection of optimization method
    integrator = KenCarp4(autodiff = true), # selection of sciml integrator
    pt_avg = 1, # numebr of the point to generate intial condition
    pt_smooth_derivative = 7,
    write_res = false,
    smoothing = false, # the smoothing is done or not?
    type_of_smoothing = "lowess",
    type_of_loss = "RE", # type of used loss
    blank_array = zeros(100), # data of all blanks
    multiple_scattering_correction = false, # if true uses the given calibration curve to fix the data
    calibration_OD_curve = "NA",  #  the path to calibration curve to fix the data
    PopulationSize = 300,
    maxiters = 2000000,
    abstol = 0.00001,
)

    # inizializing the results of sensitivity
    results_sensitivity = initialize_df_results(model)

    if write_res == true
        mkpath(path_to_results)
    end

    if multiple_scattering_correction == true
        data = correction_OD_multiple_scattering(data, calibration_OD_curve)
    end

    #defining time interval
    max_t = data[1, end]
    min_t = data[1, 1]
    tspan = (min_t, max_t)
    tsteps = data[1, :]

    # smoothing data if required
    if type_of_smoothing == "rolling_avg" && smoothing == true
        data = smoothing_data(data, pt_avg)
    elseif type_of_smoothing == "lowess" && smoothing == true
        # lowess call here
        model_fit = lowess_model(data[1, :], data[2, :], thr_lowess)
        data = Matrix(transpose(hcat(data[1, :], model_fit)))
    else
        data = copy(data)
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
            Optimization.OptimizationProblem(optf, param, u0, lb = lb_param, ub = ub_param)
        res = Optimization.solve(
            optprob_const,
            optmizator,
            PopulationSize = PopulationSize,
            maxiters = maxiters,
            abstol = abstol,
        )

        #revalution of solution for plot an loss evaluation
        remade_solution = solve(remake(ODE_prob, p = res.u), integrator, saveat = tsteps)
        sol_time = reduce(hcat, remade_solution.t)
        sol_fin = reduce(hcat, remade_solution.u)
        sol_fin = sum(sol_fin, dims = 1)
        loss_value = res.objective

        # here problem
        #max_theoretical gr
        data_th = vcat(sol_time, sol_fin)
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

#######################################################################

"""
change points functions
"""
# function from changepoints detections.jl with modification of thr usage


function getpoints_mod(profile; number_of_bin = 100)
    points_list = ()
    # list of points
    min_profile = minimum(profile)
    max_profile = maximum(profile)
    delta_profile = (max_profile - min_profile) / number_of_bin
    seq_thr = min_profile:delta_profile:max_profile

    for thr_temp in seq_thr[2:end]
        points = Int[]
        exceeded = false

        for (index, value) in enumerate(profile)
            if (value > thr_temp) && exceeded == false
                push!(points, index)
                exceeded = true
            elseif (value < thr_temp) && exceeded == true
                exceeded = false
            end
        end

        if thr_temp == seq_thr[2]
            points_list = [points]
        else
            points_list = push!(points_list, points)
        end
    end

    return reverse(points_list[1:(end-1)])
end


function cpd_local_detection(
    data::Matrix{Float64},
    n_max_cp::Int;
    type_of_detection = "lsdd",
    type_of_curve = "original",
    pt_derivative = 0,
    size_win = 2,
    method = "peaks_prominence",
    number_of_bin = 40,
)


    """
        perform the analyses with the required change point detection algorithm
        type_of_detection="lsdd" or piecewise linear fitting on the specific growth rate
        pt_derivative number of point to evaluate the derivative/specific gr (if 0 numerical derivative if >1 specific gr with that size of sliding window)
        size_win Int size of the used window in all of the methods
    """


    if type_of_detection == "lsdd" && type_of_curve == "deriv"
        list_of_cdps = cpd_lsdd_profile(
            data,
            n_max_cp;
            type_of_curve = "deriv",
            pt_deriv = pt_derivative,
            window_size = size_win,
            method = method,
            number_of_bin = number_of_bin,
        )
    elseif type_of_detection == "lsdd" && type_of_curve != "deriv"
        list_of_cdps = cpd_lsdd_profile(
            data,
            n_max_cp;
            pt_deriv = pt_derivative,
            window_size = size_win,
            type_of_curve = "original",
            method = method,
            number_of_bin = number_of_bin,
        )
    else
        type_of_detection != "lsdd"
        list_of_cdps = detect_list_change_points(
            data,
            n_max_cp;
            pt_deriv = pt_derivative,
            type_of_curve = type_of_curve,
            win_size = size_win,
            method = method,
            number_of_bin = number_of_bin,
        )
    end

    return list_of_cdps
end


function cpd_lsdd_profile(
    data::Matrix{Float64},
    n_max::Int;
    window_size = 2,
    type_of_curve = "original",
    pt_deriv = 0,
    method = "peaks_prominence",
    number_of_bin = 40,
)
    """
    evaluate change points using peak detection on a lsdd profile
    type_of_detection="lsdd" or piecewise linear fitting on the specific growth rate
    pt_derivative number of point to evaluate the derivative/specific gr (if 0 numerical derivative if >1 specific gr with that size of sliding window)
    size_win Int size of the used window in all of the methods
    """

    # evaluating the profile of lsdd on the data or on the derivative of the data
    if type_of_curve == "deriv"
        deriv = specific_gr_evaluation(data, pt_deriv)
        profile = ChangePointDetection.lsdd_profile(deriv; window = window_size)
    else
        profile = ChangePointDetection.lsdd_profile(data[2, :]; window = window_size)
    end

    # adding time to profile
    profile = convert.(Float64, profile)
    data_dissim = Matrix(transpose(hcat(data[1, 1:length(profile)], profile)))
    selected_change_point_index =
        peaks_detection(data_dissim, n_max; method = method, number_of_bin = number_of_bin)

    return selected_change_point_index
end



function detect_list_change_points(
    data::Matrix{Float64},
    n_max::Int;
    win_size = 2,
    method = "peaks_prominence",
    number_of_bin = 40,
    type_of_curve = "original",
    pt_deriv = 7,
)
    """
     evaluate change points using piecewise linear fitting
     n_max
     size_win Int size of the used window in all of the methods
    """
    if type_of_curve == "deriv"
        data = specific_gr_evaluation(data, pt_deriv)
    end

    curve_dissimilitary_deriv = curve_dissimilitary_lin_fitting(
        data, # dataset first row times second row OD
        1, # index of start
        win_size, # size sliding window
    )
    data_dissim = Matrix(
        transpose(
            hcat(
                data[1, convert.(Int, curve_dissimilitary_deriv[1, :])],
                curve_dissimilitary_deriv[2, :],
            ),
        ),
    )
    selected_change_point_index =
        peaks_detection(data_dissim, n_max; method = method, number_of_bin = number_of_bin)

    return selected_change_point_index
end

function peaks_detection(
    data::Matrix{Float64},
    n_max::Int;
    method = "peaks_prominence",
    number_of_bin = 40,
)

    """
    peaks detection
    n_max maximum number of peaks
    size_win Int size of the used window in all of the methods
    """
    if method == "peaks_prominence"
        index_of_peaks = findmaxima(data[2, :]; strict = true)
        array_prominence = peakproms(index_of_peaks[1], data[2, :])[2]
        index_prominence = peakproms(index_of_peaks[1], data[2, :])[1]

        if length(array_prominence) < n_max
            println("Warning: the max number of peaks is too much")
            top_prominence = sort(array_prominence)
        else
            top_prominence = sort(array_prominence)[((end-n_max)+1):end]
        end

        index_top_peaks = [findall(array_prominence .== i)[1] for i in top_prominence]
        selected_change_point_index = index_prominence[index_top_peaks]
        times_top_peaks = data[1, selected_change_point_index]
        values_top_peaks = data[2, selected_change_point_index]
    end

    if method == "thr_scan"
        if n_max == 1
            selected_change_point_index =
                findfirst(x -> x == maximum(data[2, :]), data[2, :])
        else
            selected_change_point_list =
                getpoints_mod(data[2, :], number_of_bin = number_of_bin)
            lenght_cdp_list = length.(selected_change_point_list)

            if n_max > maximum(lenght_cdp_list)
                println(
                    "Warning: this number of peaks is to much selecting the max number detected",
                )
                selected_change_point_index = selected_change_point_list[end]
            else
                selected_change_point_index =
                    selected_change_point_list[maximum(findlast(lenght_cdp_list .<= n_max))]

                if length(selected_change_point_index) != n_max
                    println(
                        "Warning: this number of peaks is not detected changing to nearest one smaller",
                    )
                end
            end
        end
        times_top_peaks = data[1, selected_change_point_index]
        values_top_peaks = data[2, selected_change_point_index]
    end

    return selected_change_point_index, times_top_peaks, values_top_peaks
end

function curve_dissimilitary_lin_fitting(
    data::Matrix{Float64}, # dataset first row times second row OD
    start_time_Index::Int,
    size_wind::Int, # size sliding window
)
    discrepancy_measure_curve = [start_time_Index, 0.0]
    ending = convert(Int, (length(data[2, :]) - floor(size_wind / 2) * 2))

    for index_t = start_time_Index:ending
        # defining the window
        middle_index = convert(Int, (index_t + floor(size_wind / 2)))
        end_index = convert(Int, (index_t + floor(size_wind / 2) * 2))
        win_1_data = data[2, index_t:middle_index]
        win_2_data = data[2, middle_index:end_index]
        win_tot_data = data[2, index_t:end_index]
        win_1_times = data[1, index_t:middle_index]
        win_2_times = data[1, middle_index:end_index]
        win_tot_times = data[1, index_t:end_index]

        #fitting total data
        data_total = Matrix(transpose(hcat(win_tot_times, win_tot_data)))
        fit_total = CurveFit.curve_fit(LinearFit, data_total[1, :], data_total[2, :])

        # residual calculation
        res_total = sum([
            abs((
                data_total[2, ll] - fit_total.coefs[2] * data_total[1, ll] -
                fit_total.coefs[1]
            )) for ll = 1:length(data_total[1, :])
        ])

        #fitting win 1
        data_1 = Matrix(transpose(hcat(win_1_times, win_1_data)))
        fit_1 = CurveFit.curve_fit(LinearFit, data_1[1, :], data_1[2, :])

        # residual calculation
        res_win_1 = sum([
            abs((data_1[2, ll] - fit_1.coefs[2] * data_1[1, ll] - fit_1.coefs[1])) for
            ll = 1:length(data_1[1, :])
        ])

        #fitting win 2
        data_2 = Matrix(transpose(hcat(win_2_times, win_2_data)))
        fit_2 = CurveFit.curve_fit(LinearFit, data_2[1, :], data_2[2, :])

        # residual calculation
        res_win_2 = sum([
            abs((data_2[2, ll] - fit_2.coefs[2] * data_2[1, ll] - fit_2.coefs[1])) for
            ll = 1:length(data_2[1, :])
        ])

        #evaluation of the cost
        cost = res_total - res_win_1 - res_win_2
        discrepancy_measure_curve =
            hcat(discrepancy_measure_curve, [index_t + floor(size_wind / 2), cost])
        # stop when first change point is fitted
    end

    return discrepancy_measure_curve
end

#######################################################################

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
    type_of_loss = "L2", # type of used loss
    optmizator = BBO_adaptive_de_rand_1_bin_radiuslimited(), # selection of optimization method
    integrator = KenCarp4(autodiff = true), # selection of sciml integrator
    type_of_detection = "lsdd",
    type_of_curve = "original",
    smoothing = false,
    type_of_smoothing = "lowess",
    thr_lowess = 0.05,
    pt_avg = 1,
    save_plot = false, # do plots or no
    display_plots = false,
    path_to_plot = "NA", # where save plots
    win_size = 2, # numebr of the point to generate intial condition
    pt_smooth_derivative = 7,
    multiple_scattering_correction = false, # if true uses the given calibration curve to fix the data
    calibration_OD_curve = "NA", #  the path to calibration curve to fix the data
    beta_smoothing_ms = 2.0, #  parameter of the AIC penality
    method_peaks_detection = "peaks_prominence",
    n_bins = 40,
    PopulationSize = 300,
    maxiters = 2000000,
    abstol = 0.000,
)

    if type_of_smoothing == "rolling_avg" && smoothing == true
        data_testing = smoothing_data(data_testing, pt_avg)
    elseif type_of_smoothing == "lowess" && smoothing == true
        # lowess call here
        model_fit = lowess_model(data_testing[1, :], data_testing[2, :], thr_lowess)
        data_testing = Matrix(transpose(hcat(data_testing[1, :], model_fit)))
    else
        data_testing = copy(data_testing)
    end

    if multiple_scattering_correction == true
        data_testing = correction_OD_multiple_scattering(data_testing, calibration_OD_curve)
    end

    list_change_points_dev = cpd_local_detection(
        data_testing,
        n_change_points;
        type_of_detection = type_of_detection,
        type_of_curve = type_of_curve,
        pt_derivative = pt_smooth_derivative,
        size_win = win_size,
        method = method_peaks_detection,
        number_of_bin = n_bins,
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
            optmizator = optmizator, # selection of optimization method
            integrator = integrator, # selection of sciml integrator
            pt_avg = pt_avg, # number of the point to generate intial condition
            beta_penality = beta_smoothing_ms, # penality for AIC evaluation
            smoothing = smoothing, # the smoothing is done or not?
            type_of_loss = type_of_loss, # type of used loss
            blank_array = zeros(100), # data of all blanks
            display_plot_best_model = false, # one wants the results of the best fit to be plotted
            path_to_plot = "NA",
            pt_smooth_derivative = pt_smooth_derivative,
            multiple_scattering_correction = false, # if true uses the given calibration curve to fix the data
            calibration_OD_curve = "NA", #  the path to calibration curve to fix the data
            verbose = false,
            PopulationSize = PopulationSize,
            maxiters = maxiters,
            abstol = abstol,
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
        sol_fin = sum(sol_fin, dims = 1)
        value_bonduary = remade_solution.t[end]
        time_bonduary = sol_fin[end]
        bc = [value_bonduary, time_bonduary]

        if (
            pt_smooth_derivative > (length(data_temp[1, :]) - 2) &&
            length(data_temp[1, :]) > 3
        )
            emp_max_gr_of_segment =
                maximum(specific_gr_evaluation(data_temp, pt_smooth_derivative))
            data_th = Matrix(vcat(time_sol, sol_fin))
            th_max_gr_of_segment =
                maximum(specific_gr_evaluation(data_th, pt_smooth_derivative))
            temp_res_win = vcat(temp_res_win, th_max_gr_of_segment)
            temp_res_win = vcat(temp_res_win, emp_max_gr_of_segment)
        elseif length(data_temp[1, :]) <= 3
            emp_max_gr_of_segment = "NA"
            th_max_gr_of_segment = "NA"
            temp_res_win = vcat(temp_res_win, th_max_gr_of_segment)
            temp_res_win = vcat(temp_res_win, emp_max_gr_of_segment)
        else
            emp_max_gr_of_segment =
                maximum(specific_gr_evaluation(data_temp, pt_smooth_derivative))
            data_th = Matrix(vcat(time_sol, sol_fin))
            th_max_gr_of_segment =
                maximum(specific_gr_evaluation(data_th, pt_smooth_derivative))
            temp_res_win = vcat(temp_res_win, th_max_gr_of_segment)
            temp_res_win = vcat(temp_res_win, emp_max_gr_of_segment)
        end

        if i == 2
            composed_time = copy(time_sol)
            composed_sol = copy(sol_fin)
            temp_res_win = push!(temp_res_win, i - 1)
            param_out = push!(param_out, temp_res_win)
        else
            composed_time = hcat(composed_time, reduce(hcat, time_sol))
            composed_sol = hcat(composed_sol, reduce(hcat, sol_fin))
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
            xlabel = "Time",
            ylabel = "Arb. Units",
            label = ["Data " nothing],
            markersize = 1,
            color = :black,
            title = string(label_exp, " ", name_well),
        ),
    )
    if_display(
        Plots.vline!(
            interval_changepoints[2:end],
            c = :black,
            label = ["change points" nothing],
        ),
    )
    if_display(
        Plots.plot!(
            reduce(vcat, composed_time),
            reduce(vcat, composed_sol),
            xlabel = "Time",
            ylabel = "Arb. Units",
            label = [" fitting " nothing],
            color = :red,
            title = string(label_exp, " fitting ", name_well),
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
    optmizator = BBO_adaptive_de_rand_1_bin_radiuslimited(), # selection of optimization method
    integrator = KenCarp4(autodiff = true), # selection of sciml integrator
    type_of_loss = "L2", # type of used loss
    type_of_detection = "lsdd",
    type_of_curve = "original",
    pt_avg = 1, # number of the point to generate intial condition
    smoothing = true, # the smoothing is done or not?
    save_plot = false, # do plots or no
    display_plot = false, # do plots or no
    path_to_plot = "NA", # where save plots
    path_to_results = "NA",
    win_size = 2, # numebr of the point to generate intial condition
    pt_smooth_derivative = 7,
    penality_parameter = 2.0,
    multiple_scattering_correction = "false", # if true uses the given calibration curve to fix the data
    calibration_OD_curve = "NA",  #  the path to calibration curve to fix the data
    save_all_model = false,
    method_peaks_detection = "peaks_prominence",
    n_bins = 40,
    PopulationSize = 300,
    maxiters = 2000000,
    abstol = 0.00001,
    type_of_smoothing = "lowess",
    thr_lowess = 0.05,
)

    # fitting single models
    length_dataset = length(data_testing[1, :])
    change_point_list = Vector{Vector{Any}}()

    res = ODE_Model_selection(
        data_testing, # dataset first row times second row OD
        name_well, # name of the well
        label_exp, #label of the experiment
        list_of_models, # ode model to use
        list_lb_param, # lower bound param
        list_ub_param; # upper bound param
        optmizator = optmizator, # selection of optimization method
        integrator = integrator, # selection of sciml integrator
        pt_avg = pt_avg, # number of the point to generate intial condition
        beta_penality = penality_parameter, # penality for AIC evaluation
        smoothing = smoothing, # the smoothing is done or not?
        type_of_loss = type_of_loss, # type of used loss
        blank_array = zeros(100), # data of all blanks
        display_plot_best_model = false, # one wants the results of the best fit to be plotted
        save_plot_best_model = false,
        path_to_plot = "NA",
        pt_smooth_derivative = pt_smooth_derivative,
        multiple_scattering_correction = multiple_scattering_correction, # if true uses the given calibration curve to fix the data
        calibration_OD_curve = calibration_OD_curve, #  the path to calibration curve to fix the data
        verbose = false,
        PopulationSize = PopulationSize,
        maxiters = maxiters,
        abstol = abstol,
        type_of_smoothing = type_of_smoothing,
        thr_lowess = thr_lowess,
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
                type_of_loss = type_of_loss, # type of used loss
                optmizator = optmizator, # selection of optimization method
                integrator = integrator, # selection of sciml integrator
                type_of_detection = type_of_detection,
                type_of_curve = type_of_curve,
                smoothing = smoothing,
                pt_avg = pt_avg,
                save_plot = false, # do plots or no
                display_plots = false, # do plots or no
                path_to_plot = "NA", # where save plots
                win_size = win_size, # numebr of the point to generate intial condition
                pt_smooth_derivative = pt_smooth_derivative,
                multiple_scattering_correction = multiple_scattering_correction, # if true uses the given calibration curve to fix the data
                calibration_OD_curve = calibration_OD_curve, #  the path to calibration curve to fix the data
                beta_smoothing_ms = penality_parameter, #  parameter of the AIC penality
                method_peaks_detection = method_peaks_detection,
                n_bins = n_bins,
                type_of_smoothing = type_of_smoothing,
                thr_lowess = thr_lowess,
            )

            # composing piecewise penality
            loss_penality = sum([
                direct_search_results[1][kk][(end-2)] for
                kk = 1:length(direct_search_results[1])
            ])
            n_param =
                sum([
                    length(direct_search_results[1][kk][3:(end-5)]) for
                    kk = 1:length(direct_search_results[1])
                ]) + n_change_points
            new_penality =
                length_dataset * log(loss_penality / length_dataset) +
                penality_parameter * n_param

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
            xlabel = "Time",
            ylabel = "Arb. Units",
            label = ["Data " nothing],
            markersize = 1,
            color = :black,
            title = string(label_exp, " ", name_well),
        ),
    )
    if_display(
        Plots.vline!(
            change_point_to_plot[2:end],
            c = :black,
            label = ["change points" nothing],
        ),
    )
    if_display(
        Plots.plot!(
            reduce(vcat, time_points_to_plot),
            reduce(vcat, sol_to_plot),
            xlabel = "Time",
            ylabel = "Arb. Units",
            label = [" fitting " nothing],
            color = :red,
            title = string(label_exp, " fitting ", name_well),
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

"
Testing part
"
function initialize_res_segmentation(
    data_df::Matrix{Float64},
    list_of_model_parameters::Any,
    number_of_segment::Any,
)

    # evaluation of the number of columns
    ncols = size(data_df, 2) * (number_of_segment + 1) + 1
    # evaluation of the number of rows
    nrow = maximum(length.(list_of_model_parameters)) + 7

    # inizialization of the matrix as full of missing
    matrix_result = missings(Any, nrow, ncols)
    # generation of the names of the rows
    matrix_result[1] = "well"
    matrix_result[2] = "label_exp"
    matrix_result[3] = "model"
    matrix_result[(end-3)] = "th_gr"
    matrix_result[(end-2)] = "em_gr"
    matrix_result[(end-1)] = "loss"

    for i = 4:(4+maximum(length.(list_of_model_parameters)))
        matrix_result[i, 1] = string("param_", i)
    end

    return matrix_result
end

# function to initialize the df for results of  model selection
function initialize_res_model_selection(list_of_model_parameters::Any)
  
    nmax_param = maximum(length.(list_of_model_parameters))
    # evaluation of the number of rows
    nrow = nmax_param + 6
    # inizialization of the matrix as full of missing
    matrix_result = missings(Any, nrow)
    # generation of the names of the rows
    matrix_result[1] = "well"
    matrix_result[2] = "label_exp"
    matrix_result[3] = "model"
    matrix_result[(end-2)] = "loss"
    matrix_result[(end-1)] = "th_gr"
    matrix_result[(end)] = "em_gr"

    for i = 4:(4+nmax_param-1)
        matrix_result[i] = string("param_", i - 3)
    end

    return matrix_result
end

# function to initialize the df for results of  model selection
function initialize_res_segmentation(list_of_model_parameters::Any)

    nmax_param = maximum(length.(list_of_model_parameters))

    # evaluation of the number of columns
    # evaluation of the number of rows
    nrow = nmax_param + 7

    # inizialization of the matrix as full of missing
    matrix_result = missings(Any, nrow)

    # generation of the names of the rows
    matrix_result[1] = "well"
    matrix_result[2] = "label_exp"
    matrix_result[3] = "model"
    matrix_result[(end-3)] = "loss"
    matrix_result[(end-2)] = "th_gr"
    matrix_result[(end-1)] = "em_gr"
    matrix_result[end] = "segment"

    for i = 4:(4+nmax_param-1)
        matrix_result[i] = string("param_", i - 3)
    end

    return matrix_result
end

# given optimization results for model selection add missing to mach the size  of the model of
function expand_res_model_selection(
    param_res::Any,
    list_of_model_parameters::Any,
    names_of_the_well::String,
)

    n_param = length(param_res) - 5
    nmax_param = maximum(length.(list_of_model_parameters))
    temp_output = missings(Any, nmax_param + 6)
    temp_output[1] = names_of_the_well
    temp_output[2] = param_res[1]
    temp_output[3] = param_res[2]
    temp_output[(end-2)] = param_res[(end)]
    temp_output[(end-1)] = param_res[(end-2)]
    temp_output[(end)] = param_res[(end)-1]

    for i = 3:(3+n_param)
        temp_output[i] = param_res[i-1]
    end

    return temp_output
end

# here the function to construct segmentation
function expand_res_segmentation(
    param_res::Any,
    list_of_model_parameters::Any,
    number_of_segment::Int,
    names_of_the_well::String,
    label_exp::String,
)

    nmax_param = maximum(length.(list_of_model_parameters))
    fin_output = Matrix{Any}

    for s = 1:number_of_segment
        n_param = length(param_res[s]) - 5
        temp_output = missings(Any, nmax_param + 7)
        temp_output[1] = names_of_the_well
        temp_output[2] = label_exp
        temp_output[3] = param_res[s][1]
        temp_output[(end-3)] = param_res[s][(end-3)]
        temp_output[(end-2)] = param_res[s][(end-2)]
        temp_output[(end-1)] = param_res[s][(end-1)]
        temp_output[(end)] = param_res[s][(end)]

        for i = 4:(4+n_param-1)
            temp_output[i] = param_res[s][i-2]
        end

        if s == 1
            fin_output = temp_output
        else
            fin_output = hcat(fin_output, temp_output)
        end
    end

    return fin_output
end

function initialize_df_results_ode_custom( list_of_model_parameters::Any)

    nmax_param =  length(list_of_model_parameters)

    # evaluation of the number of columns
    # evaluation of the number of rows
    nrow = length(list_of_model_parameters) + 5

    # inizialization of the matrix as full of missing   
    matrix_result = missings(Any, nrow)

    # generation of the names of the rows
    matrix_result[1] = "well"
    matrix_result[2] = "label_exp"
    matrix_result[(end - 2) ] = "th_gr"
    matrix_result[(end - 1) ] =  "em_gr"  
    matrix_result[(end)] =    "loss" 
    
    for i in 3:(3 + nmax_param-1)
                                            
        matrix_result[i] = string("param_",i-2)
    end
    return matrix_result
end    