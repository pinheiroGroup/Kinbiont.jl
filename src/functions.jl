using OrdinaryDiffEq
using Missings
using Statistics
using Interpolations
using DataFrames
using OptimizationMultistartOptimization
using SciMLBase
#######################################################################

include("ODE_models.jl");
include("Fit_one_file_functions.jl");
include("pre_processing_functions.jl");
include("Fit_one_well_functions.jl");
include("loss_list.jl");
include("NL_models.jl");
include("NL_fit_one_well.jl");
include("NL_fit_one_file.jl");
include("NL_loss_list.jl");
include("cpd_functions.jl");
include("ML_downstream.jl");
include("data_struct_Kimchi.jl");


function model_selector(model::String, u0, tspan, param=nothing)


    ODE_prob = ODEProblem(models[model].func, u0, tspan, param)

    return ODE_prob
end

"""
    specific_gr_evaluation(data_smooted::Any, 
    pt_smoothing_derivative::Int)


Function that evalauates specific gr evaluation with slinding window log-lin fitting

# Arguments:

- 'data_smooted':  matrix of data 2xn_time points, it is a single curve.
- 'pt_smoothing_derivative': Int size of the win, if <2 the the numerical derivative of (log) data is evaluate with interpolation algorithm

# Output:

- 'specific_gr' an array with the specific growth rate 
"""
function specific_gr_evaluation(data_smooted::Any, pt_smoothing_derivative::Int)

    specific_gr = 0.0
    if pt_smoothing_derivative > 1
        for r = 1:1:(eachindex(data_smooted[2, :])[end].-pt_smoothing_derivative)
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

function deriv_evaluation(data_testing)
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
    res_param = [
        [model, well_name]
        res[1:length(models[model].params)-5]
        [th_gr, em_gr, loss]
    ]

    return res_param
end


function initialize_df_results(model::String)
    param_names = models[model].params
    return param_names
end


function guess_param(lb_param::Vector{Float64}, ub_param::Vector{Float64})
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
        [(ub_param[i] - lb_param[i]) / (N_step_morris + 1) for i in eachindex(ub_param)]

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

function analyze_segment(label_exp, name_well, data, segment_number, pt_smoothing_derivative)
    if length(data[2, :]) >= pt_smoothing_derivative + 1
        specific_gr = specific_gr_evaluation(data, pt_smoothing_derivative)

        specific_gr_times = [
            (data[1, r] + data[1, (r+pt_smoothing_derivative)]) / 2 for
            r = 1:1:(eachindex(data[2, :])[end].-pt_smoothing_derivative)
        ]

        max_specific_gr = maximum(specific_gr)
        min_specific_gr = minimum(specific_gr)


        t_of_max = specific_gr_times[argmax(max_specific_gr)]
        index_of_max_od = findfirst(data[1, :] .> t_of_max)
        od_of_max = data[2, index_of_max_od]

        derivative = deriv_evaluation(data)
        min_deriv = maximum(derivative)
        max_deriv = minimum(derivative)




        res_deriv = [label_exp, name_well, max_specific_gr, min_specific_gr, t_of_max, od_of_max, max_deriv, min_deriv, data[1, end], data[2, end], segment_number]
    else
        res_deriv = [label_exp, name_well, missing, missing, missing, missing, missing, missing, data[1, end], data[2, end], segment_number]
        specific_gr = [missing]
        specific_gr_times = [missing]
        derivative = [missing]
    end
    return res_deriv, specific_gr, specific_gr_times, derivative


end

"""

    ODE_sim(
    model::String, 
    n_start::Vector{Float64}, 
    tstart::Float64, 
    tmax::Float64,
    delta_t::Float64, 
    param_of_ode::Vector{Float64};
    integrator=KenCarp4(),
    )


This function performs an ODE simulation of a model

# Arguments:
- `model::String`: The model to simulate. For the possible options please check the documentation.
- `n_start::Vector{Float64}`: The starting conditions.
- `tstart::Float64`: The start time of the simulation.
- `tmax::Float64`: The final time of the simulation.
- `delta_t::Float64`: The time step of the output.
- `param_of_ode::Vector{Float64}`: The parameters of the ODE model.
     
# Key argument:
- `integrator=KenCarp4() `: The chosen solver from the SciML ecosystem for ODE integration, default KenCarp4 algorithm. 

# Output:
    
- it returns a standard SciML output (i.e., if `sim =ODE_sim(...)`, then `sim.t` is the array of times and `sim.u` is the array of the simulation)
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

export ODE_sim

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


"""
    stochastic_sim(
    model::String,
    n_start::Int, 
    n_mol_start::Float64, 
    tstart::Float64, 
    tmax::Float64, 
    delta_t::Float64, 
    k_1_val::Float64,
    k_2_val::Float64, 
    alpha_val::Float64, 
    lambda::Float64, 
    n_mol_per_birth::Float64,
    volume::Float64,
    )


This function performs a stochastic simulation of a model, considering cell growth and nutrient consumption over time.

# Arguments:

- `model::String`: The model to simulate. Possible options "Monod","Haldane","Blackman","Tessier","Moser","Aiba-Edwards", and "Verhulst"
- `n_start::Int`: The number of starting cells.
- `n_mol_start::Float64`: The starting concentration of the limiting nutrient.
- `tstart::Float64`: The start time of the simulation.
- `tmax::Float64`: The final time of the simulation.
- `delta_t::Float64`: The time step for the Poisson approximation.
- `k_1_val::Float64`: The value of parameter k1.
- `k_2_val::Float64`: The value of the Monod constant.
- `alpha_val::Float64`: The maximum possible growth rate.
- `lambda::Float64`: The lag time, simulated as a zero growht time span at the start
- `n_mol_per_birth::Float64`: The nutrient consumed per division (mass).
- `volume::Float64`: The volume.


# Output (if `sim =stochastic_sim(...)`):

- `sim[1]`: array of the number of individuals in the population.
- `sim[2]`: array of the number of biomass equivalent mass of the limiting nutrient concentration.
- `sim[3]`: array of the times of the simulation. 
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
        n_birth = rand(Distributions.Poisson(total_birth_rate), 1)

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



function expand_res_seg(
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
    elseif number_of_segment > 0
        nmax_param = maximum(length.(list_of_model_parameters))
        fin_output = Matrix{Any}

        for s = 1:number_of_segment
            n_param = length(param_res[s]) - 7
            temp_output = missings(Any, nmax_param + 7)
            temp_output[1] = names_of_the_well
            temp_output[2] = label_exp
            temp_output[3] = param_res[s][3]
            temp_output[(end-3)] = param_res[s][(end-3)]
            temp_output[(end-2)] = param_res[s][(end-2)]
            temp_output[(end-1)] = param_res[s][(end-1)]
            temp_output[(end)] = param_res[s][(end)]

            for i = 4:(4+n_param-1)
                temp_output[i] = param_res[s][i]
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
    elseif number_of_segment > 0
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


function AICc_evaluation(n_param, beta_penality, data, data_th; correction=true)
    n_data = length(data)


    if n_data == length(data_th)
        if n_data > n_param - 2
            # RSS = sum(abs.(data_th .- data ) .^ 2)
            RSS = sum(abs.(data_th .- data) .^ 2)

            println(n_param)
            println(RSS)
            println(log(RSS / n_data))
            if correction == true
                correction_value = beta_penality * (((n_param + 1) * (n_param + 2)) / (n_data - n_param - 2))
            else
                correction_value = 0.0
            end
            AIC = +beta_penality * n_param + n_data * log(RSS / n_data)
            AICc = AIC + correction_value
        else
            AICc = 10^9

        end
    else
        AICc = 10^9


    end
    println(AICc)

    return AICc

end




function AICc_evaluation2(n_param, beta_penality, data, loss; correction=true)

    n_data = length(data)
    if n_data > n_param - 2
        # RSS = sum(abs.(data_th .- data ) .^ 2)
        RSS = loss

        if correction == true
            correction_value = beta_penality * (((n_param + 1) * (n_param + 2)) / (n_data - n_param - 2))
        else
            correction_value = 0.0
        end
        AIC = +beta_penality * n_param + n_data * log(RSS)
        AICc = AIC + correction_value
    else
        AICc = 10^9

    end




    return AICc

end



function remove_replicate_data(composed_time, composed_sol)
    duplicates_index = [0]
    for k = 2:length(composed_time)
        if composed_time[k-1] == composed_time[k]
            duplicates_index = vcat(k - 1, duplicates_index)

        end
    end
    if length(duplicates_index) > 1
        index_tot = 1:1:length(composed_time)
        index_to_use = setdiff(index_tot, duplicates_index)
        composed_time = composed_time[index_to_use]
        composed_sol = composed_sol[index_to_use]

    end

    return composed_time, composed_sol

end


function reading_annotation(path_to_annotation::Any)


    if typeof(path_to_annotation) == String

        annotation = CSV.File(string(path_to_annotation), header=false)
        names_of_annotated_df = [annotation[l][1] for l in eachindex(annotation)]
        # selecting blank wells
        properties_of_annotation = [annotation[l][2] for l in eachindex(annotation)]
        list_of_blank = names_of_annotated_df[findall(x -> x == "b", properties_of_annotation)]
        list_of_discarded =
            names_of_annotated_df[findall(x -> x == "X", properties_of_annotation)]
        list_of_blank = Symbol.(list_of_blank)
        list_of_discarded = Symbol.(list_of_discarded)
    else
        names_of_annotated_df = [""]
        properties_of_annotation = [""]
        list_of_blank = []
        list_of_discarded = []
    end


    return names_of_annotated_df, properties_of_annotation, list_of_blank, list_of_discarded
end


function KimchiSolve(loss_function,
    u0,
    p;
    opt,
    auto_diff_method=nothing,
    multistart=false,
    n_restart=50,
    cons=nothing,
    opt_params...)

    # generation of optimization fuction

    if isnothing(auto_diff_method) == true && isnothing(cons) == true

        optf = Optimization.OptimizationFunction((x, p) -> loss_function(x))

    elseif isnothing(auto_diff_method) == false && isnothing(cons) == true

        optf = Optimization.OptimizationFunction((x, p) -> loss_function(x), auto_diff_method)


    elseif isnothing(auto_diff_method) == true && isnothing(cons) == false

        optf = Optimization.OptimizationFunction((x, p) -> loss_function(x), cons=cons)


    else
        isnothing(auto_diff_method) == false && isnothing(cons) == false


        optf = Optimization.OptimizationFunction((x, p) -> loss_function(x), auto_diff_method, cons=cons)


    end


    prob = OptimizationProblem(optf, p, u0; opt_params...)


    if multistart == true

        sol = solve(prob, MultistartOptimization.TikTak(n_restart), opt)
        
    else

        sol = solve(prob, opt)

    end

    return sol
end


function KimchiSolve_NL(loss_function,
    u0,
    data;
    opt,
    auto_diff_method=nothing,
    multistart=false,
    n_restart=50,
    cons=nothing,
    opt_params...)

    # generation of optimization fuction

    if isnothing(auto_diff_method) == true && isnothing(cons) == true

        optf = Optimization.OptimizationFunction( loss_function)

    elseif isnothing(auto_diff_method) == false && isnothing(cons) == true

        optf = Optimization.OptimizationFunction( loss_function, auto_diff_method)


    elseif isnothing(auto_diff_method) == true && isnothing(cons) == false

        optf = Optimization.OptimizationFunction(loss_function, cons=cons)


    else
        isnothing(auto_diff_method) == false && isnothing(cons) == false


        optf = Optimization.OptimizationFunction(loss_function, auto_diff_method, cons=cons)


    end



    if multistart == true
        prob = OptimizationProblem(optf, u0 ; opt_params...)

        sol = solve(prob, MultistartOptimization.TikTak(n_restart), opt)
        
    else
        prob = OptimizationProblem(optf, u0; opt_params...)

        sol = solve(prob, opt)

    end

    return sol
end

export reading_annotation
export specific_gr_evaluation
export stochastic_sim
export ODE_sim
export initialize_df_results
export initialize_df_results_ode_custom
export expand_res
export expand_res_seg
export KimchiSolve
export generating_IC
export KimchiSolve_NL