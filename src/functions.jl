

#######################################################################
#

include("models.jl");
include("Fit_one_file_functions.jl");
include("pre_processing_functions.jl");
include("Fit_one_well_functions.jl");
include("loss_list.jl");
include("NL_models.jl");
include("NL_fit_one_well.jl");
include("NL_loss_list.jl");

"""
Internal functions
"""

function model_selector(model::String, u0, tspan, param=nothing)
    """
      generate sciML ODE problem
    """

    ODE_prob = ODEProblem(models[model].func, u0, tspan, param)

    return ODE_prob
end

function NL_model_selector(model::String, u0, tspan, param=nothing)
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
        for r = 1:1:(length(data_smooted[2, :])-pt_smoothing_derivative)
            M = [ones(pt_smoothing_derivative) data_smooted[1, r:(r+pt_smoothing_derivative-1)]]
            Y = log.(data_smooted[2, r:(r+pt_smoothing_derivative-1)])

            fit = M \ Y
            if r == 1
                specific_gr = fit[2]
            else
                specific_gr = vcat(specific_gr, fit[2])

            end
        end
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
    integrator=KenCarp4(), # which sciml solver of ode
)

    # defining time stepping
    t_steps = tstart:delta_t:tmax
    tspan = (tstart, tmax)
    u0 = n_start
    ODE_prob = model_selector(model, u0, tspan, param_of_ode)
    sim = solve(ODE_prob, integrator, saveat=t_steps)

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
    sim = solve(ODE_prob, integrator, saveat=t_steps)

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

#######################################################################

"""
change points functions
"""
# function from changepoints detections.jl with modification of thr usage


function getpoints_mod(profile; number_of_bin=100)
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
    type_of_detection="lsdd",
    type_of_curve="original",
    pt_derivative=0,
    size_win=2,
    method="peaks_prominence",
    number_of_bin=40,
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
            type_of_curve="deriv",
            pt_deriv=pt_derivative,
            window_size=size_win,
            method=method,
            number_of_bin=number_of_bin,
        )
    elseif type_of_detection == "lsdd" && type_of_curve != "deriv"
        list_of_cdps = cpd_lsdd_profile(
            data,
            n_max_cp;
            pt_deriv=pt_derivative,
            window_size=size_win,
            type_of_curve="original",
            method=method,
            number_of_bin=number_of_bin,
        )
    else
        type_of_detection != "lsdd"
        list_of_cdps = detect_list_change_points(
            data,
            n_max_cp;
            pt_deriv=pt_derivative,
            type_of_curve=type_of_curve,
            win_size=size_win,
            method=method,
            number_of_bin=number_of_bin,
        )
    end

    return list_of_cdps
end


function cpd_lsdd_profile(
    data::Matrix{Float64},
    n_max::Int;
    window_size=2,
    type_of_curve="original",
    pt_deriv=0,
    method="peaks_prominence",
    number_of_bin=40,
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
        profile = ChangePointDetection.lsdd_profile(deriv; window=window_size)
    else
        profile = ChangePointDetection.lsdd_profile(data[2, :]; window=window_size)
    end

    # adding time to profile
    profile = convert.(Float64, profile)
    data_dissim = Matrix(transpose(hcat(data[1, 1:length(profile)], profile)))
    selected_change_point_index =
        peaks_detection(data_dissim, n_max; method=method, number_of_bin=number_of_bin)

    return selected_change_point_index
end



function detect_list_change_points(
    data::Matrix{Float64},
    n_max::Int;
    win_size=2,
    method="peaks_prominence",
    number_of_bin=40,
    type_of_curve="original",
    pt_deriv=7,
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
        peaks_detection(data_dissim, n_max; method=method, number_of_bin=number_of_bin)

    return selected_change_point_index
end

function peaks_detection(
    data::Matrix{Float64},
    n_max::Int;
    method="peaks_prominence",
    number_of_bin=40,
)

    """
    peaks detection
    n_max maximum number of peaks
    size_win Int size of the used window in all of the methods
    """
    if method == "peaks_prominence"
        index_of_peaks = findmaxima(data[2, :]; strict=true)
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
                getpoints_mod(data[2, :], number_of_bin=number_of_bin)
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
        data_1 = Matrix(transpose(hcat(win_1_times, win_1_data)))
        data_2 = Matrix(transpose(hcat(win_2_times, win_2_data)))

        X_total = data_total[1, :]
        Y_total = data_total[2, :]
        X_1 = data_1[1, :]
        Y_1 = data_2[2, :]
        X_2 = data_1[1, :]
        Y_2 = data_2[2, :]
        N_1 = length(data_1[1, :])
        N_2 = length(data_2[1, :])
        N_tot = length(data_total[1, :])
        M_1 = [ones(N_1) X_1]
        M_2 = [ones(N_2) X_2]
        M_tot = [ones(N_tot) X_total]
        fit_1 = M_1 \ Y_1
        fit_2 = M_2 \ Y_2
        fit_total = M_tot \ Y_total




        # residual calculation
        res_total = sum([
            abs((
                data_total[2, ll] - fit_total[2] * data_total[1, ll] -
                fit_total[1]
            )) for ll = 1:length(data_total[1, :])
        ])


        # residual calculation
        res_win_1 = sum([
            abs((data_1[2, ll] - fit_1[2] * data_1[1, ll] - fit_1[1])) for
            ll = 1:length(data_1[1, :])
        ])

        #fitting win 2

        # residual calculation
        res_win_2 = sum([
            abs((data_2[2, ll] - fit_2[2] * data_2[1, ll] - fit_2[1])) for
            ll = 1:length(data_2[1, :])
        ])

        #evaluation of the cost
        cost = -res_total + res_win_1 + res_win_2
        discrepancy_measure_curve =
            hcat(discrepancy_measure_curve, [index_t + floor(size_wind / 2), cost])
        # stop when first change point is fitted
    end

    return discrepancy_measure_curve
end

#######################################################################



"
Testing part
"
function initialize_res_ms(
    list_of_model_parameters::Any;
    number_of_segment=0,
)
    if number_of_segment > 0

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



    elseif number_of_segment == 0

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
    end
    return matrix_result
end

function expand_res(
    param_res::Any,
    list_of_model_parameters::Any,
    names_of_the_well::String,
    label_exp::String;
    number_of_segment=0,)
    if number_of_segment == 0
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
        fin_output = copy(temp_output)
    else
        number_of_segment > 0
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


    end

    return fin_output
end



function initialize_df_results_ode_custom(list_of_model_parameters::Any)

    nmax_param = length(list_of_model_parameters)

    # evaluation of the number of columns
    # evaluation of the number of rows
    nrow = nmax_param + 5

    # inizialization of the matrix as full of missing   
    matrix_result = missings(Any, nrow)

    # generation of the names of the rows
    matrix_result[1] = "well"
    matrix_result[2] = "label_exp"
    matrix_result[(end-2)] = "th_gr"
    matrix_result[(end-1)] = "em_gr"
    matrix_result[(end)] = "loss"

    for i in 3:(3+nmax_param-1)

        matrix_result[i] = string("param_", i - 2)
    end
    return matrix_result
end


