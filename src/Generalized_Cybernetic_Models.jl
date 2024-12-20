# Define a generalized Cybernetic Model structure
struct Kinbiont_Cybernetic_Model_sim
    Bio_mass_conc::Float64             # Biomass concentration
    Substrate_concentrations::Vector{Float64}  # Array of substrate concentrations (S1, S2, ..., Sn)
    Protein_concentrations::Vector{Float64}    # Array of protein concentrations (P1, P2, ..., Pn)
    allocation_rule::Any                      # Function for the functional dependence of synthesis rate
    reaction::Union{Function, Nothing}  # Reaction function (optional, depends on the model)
    a::Vector{Union{Nothing,Float64}}               # Synthesis rate constants for proteins
    b::Vector{Union{Nothing,Float64}}                 # Degradation constants for proteins
    V_S::Vector{Union{Nothing,Float64}}              # Substrate utilization rate
    k_S::Vector{Union{Nothing,Float64}}                # Saturation constants for substrates
    Y_S::Vector{Union{Nothing,Float64}}                # Yield coefficients for cell mass per substrate
    cost::Any
    protein_thresholds::Any
end





function define_loss_function_cybernetic_ODEs(data, set_of_equation_to_fit, ODE_prob, integrator, tsteps)


    if !isnothing(set_of_equation_to_fit)
        index_of_eqs = set_of_equation_to_fit 
        index_of_data = set_of_equation_to_fit .+1 


    else

        index_of_eqs = 1:1:(size(data)[1]-1)
        index_of_data = index_of_eqs .+ 1

    end



    function loss_RE_Cybernentic(data, index_of_eqs, index_of_data, ODE_prob, integrator, p, tsteps)
        sol = solve(
            ODE_prob,
            integrator,
            p=p,
            saveat=tsteps,
            verbose=false,
            abstol=1e-10,
            reltol=1e-10,
        )


        sol_t = reduce(hcat, sol.u)



        if size(sol_t)[2] == size(data)[2]

            lossa =
                0.5 * NaNMath.sum(abs2.(1.0 .- (data[index_of_data, :] ./ sol_t[index_of_eqs, :]))) /
                length(data[2, :])

        else

            lossa = 10.0^9 * length(data[2, :])
        end


        return lossa, sol
    end



    return (p) ->
    loss_RE_Cybernentic(data, index_of_eqs, index_of_data, ODE_prob, integrator, p, tsteps)

end

# Custom constructor to allow initialization using keyword arguments
function Kinbiont_Cybernetic_Model_sim(; Bio_mass_conc, Substrate_concentrations, Protein_concentrations,
    allocation_rule, reaction, a, b, V_S, k_S, Y_S,cost,protein_thresholds)
    return Kinbiont_Cybernetic_Model_sim(Bio_mass_conc, Substrate_concentrations, Protein_concentrations,
    allocation_rule, reaction, a, b, V_S, k_S, Y_S,cost,protein_thresholds)
end





# Function for proportional allocation based on substrate concentrations
function proportional_allocation_rule(a, b, V_S, k_S, Y_S,P,S,cost,protein_thresholds)
    # Normalize substrate concentrations to create allocation vector
    alloc = S ./ sum(S)
    return alloc
end



function threshold_switching_rule(a, b, V_S, k_S, Y_S,P,S,cost,protein_thresholds)
    n = length(S)
    alloc = zeros(n)
    # Sort substrates by descending utilization rate (V_S)
    sorted_indices = sortperm(V_S, rev=true)

    for i in sorted_indices
        if S[i] > 0.000001
            alloc[i] = 1.0  # Allocate all resources to this substrate
            break
        end
    end

    # If no substrate is above the threshold, fallback to equal allocation
    if sum(alloc) == 0
        alloc .= 1.0 / n
    end

    return alloc
end



function cost_benefit_allocation_rule(a, b, V_S, k_S, Y_S,P,S,cost,protein_thresholds)
    # Calculate benefit-to-cost ratio
    benefits = V_S ./ costs
    # Normalize the allocation vector to sum to 1
    alloc = benefits ./ sum(benefits)
    return alloc
end



function feedback_controlled_allocation_rule(a, b, V_S, k_S, Y_S,P,S,cost,protein_thresholds)
    n = length(S)
    alloc = ones(n)

    for i in 1:n
        if P[i] > protein_thresholds[i]
            alloc[i] = 0.0  # Stop allocating to this substrate if protein exceeds threshold
        end
    end

    # Normalize the allocation vector to sum to 1
    alloc = alloc ./ sum(alloc)
    return alloc
end



function replace_nothing_with_param(Cybernetic_Model1)
    # Initialize a counter for cumulative param indices
    param_counter = 0

    # Define a helper function to replace `nothing` with cumulative `param[i]`
    function replace_in_array(arr)
        return [isnothing(val) ? "param[$(param_counter += 1)]" : val for val in arr]
    end

    # Replace `nothing` in each parameter array
    a_updated = replace_in_array(Cybernetic_Model1.a)
    b_updated = replace_in_array(Cybernetic_Model1.b)
    V_S_updated = replace_in_array(Cybernetic_Model1.V_S)
    k_S_updated = replace_in_array(Cybernetic_Model1.k_S)
    Y_S_updated = replace_in_array(Cybernetic_Model1.Y_S)

    return (a_updated, b_updated, V_S_updated, k_S_updated, Y_S_updated)
end




function Generate_Generalized_cybernetic_ODEs_problem(Cybernetic_Model1)
    # Extract the parameters
    n = length(Cybernetic_Model1.Substrate_concentrations)  # Number of substrates

    
    # Unpack parameters from the model
    a, b, V_S, k_S, Y_S = replace_nothing_with_param(Cybernetic_Model1)

    rule = Cybernetic_Model1.allocation_rule  # Function for synthesis rate dependency
    r = Cybernetic_Model1.reaction
   # Define allocation vector for resource allocation (normalized)
    if typeof(rule) == Vector{Float64} 
   
   
        alloc = rule./ sum(rule)  # Assuming proportional allocation based on substrate concentrations
    

    end









  function Cybernetic_Model_odes(du,u,t,param)

        # alias the state vector from the input (biomass, substrates, proteins)
        S = u[2:(n+1)]  # Substrate concentrations (S1, S2, ..., Sn)
        P = u[(n+2):end]  # Key protein concentrations (P1, P2, ..., Pn)
   
        if typeof(rule) != Vector{Float64} 
   
   
            alloc = rule(a, b, V_S, k_S, Y_S,P,S,cost,protein_thresholds)  # Assuming proportional allocation based on substrate concentrations
        
    
        end

    # Counter for params to handle dynamic replacement
    param_counter = 0

    # Function to get value or substitute param if needed
    get_val = (x) -> isa(x, String) ? param[param_counter += 1] : x


# Substrate consumption using Monod kinetics
    for i in 1:n
        V_S_val = get_val(V_S[i])
        k_S_val = get_val(k_S[i])

        du[i + 1]= -V_S_val * P[i] * S[i] / (k_S_val + S[i]) * u[1]  # Rate of change of each substrate
    end

    # Protein synthesis rate 
     for i in 1:n
        a_val = get_val(a[i])
        b_val = get_val(b[i])
        k_S_val = get_val(k_S[i])

           du[i + n + 1]  = a_val * alloc[i] * k_S_val - b_val * P[i] * u[1] # Rate of change of each key protein
      end

     
      Y_S_vec = [get_val(Y_S[i]) for i in 1:eachindex(Y_S)]
    # Cell growth dynamics (cumulative effect of substrate utilization)
            du[1] = -sum(Y_S_vec .* du[2:(n+1)] ) * u[1]   # Change in biomass concentration based on substrate utilization
    end


    return Cybernetic_Model_odes  # In-place modification of the rate of change vector

end


# Function to generate generalized cybernetic ODE system (corrected for the ODE solver)
function Generate_Generalized_cybernetic_ODEs_problem_simulation(Cybernetic_Model1)
    # Extract the parameters
    n = length(Cybernetic_Model1.Substrate_concentrations)  # Number of substrates

    
    # Unpack parameters from the model
    a = Cybernetic_Model1.a  # Synthesis rate constants for proteins
    b = Cybernetic_Model1.b  # Degradation constants for proteins
    rule = Cybernetic_Model1.allocation_rule  # Function for synthesis rate dependency
    V_S = Cybernetic_Model1.V_S  # Substrate utilization rate
    k_S = Cybernetic_Model1.k_S  # Saturation constants for substrates
    Y_S = Cybernetic_Model1.Y_S  # Yield coefficients for cell mass per substrate
    r = Cybernetic_Model1.reaction
    cost = Cybernetic_Model1.cost
    protein_thresholds = Cybernetic_Model1.protein_thresholds
    # Define allocation vector for resource allocation (normalized)
    if typeof(rule) == Vector{Float64} 
   
   
        alloc = rule./ sum(rule)  # Assuming proportional allocation based on substrate concentrations
    

    end









  function Cybernetic_Model_odes(du,u,t,param)

        # alias the state vector from the input (biomass, substrates, proteins)
        X_c = u[1]  # Biomass concentration (cells)
        S = u[2:(n+1)]  # Substrate concentrations (S1, S2, ..., Sn)
        P = u[(n+2):end]  # Key protein concentrations (P1, P2, ..., Pn)

        if typeof(rule) != Vector{Float64} 
   
   
            alloc = rule(a, b, V_S, k_S, Y_S,P,S,cost,protein_thresholds)  # Assuming proportional allocation based on substrate concentrations
        
    
        end

# Substrate consumption using Monod kinetics
    for i in 1:n
        du[i + 1]= - V_S[i] * P[i] * S[i] / (k_S[i] + S[i]) * u[1]  # Rate of change of each substrate
    end

    # Protein synthesis rate 
     for i in 1:n

           du[i + n + 1]  = a[i] * alloc[i] * k_S[i] - b[i] * P[i] * u[1] # Rate of change of each key protein
      end

    

    # Cell growth dynamics (cumulative effect of substrate utilization)
            du[1] = - sum(Y_S .* du[2:(n+1)]) *u[1]  # Change in biomass concentration based on substrate utilization
    end


    return Cybernetic_Model_odes  # In-place modification of the rate of change vector

end


function Kinbiont_Cybernetic_Model_simulation(model,tmin,tmax,deltaT; integrator = Tsit5())

    tspan =(tmin,tmax)
    tsteps = tmin:deltaT:tmax
    tsteps1 = [tsteps[i] for i in 1:length(tsteps)]
    u0 = [model.Bio_mass_conc, model.Substrate_concentrations..., model.Protein_concentrations...]
    Cybernetic_Model_odes = Generate_Generalized_cybernetic_ODEs_problem_simulation(model)
    oprob = ODEProblem(Cybernetic_Model_odes, u0, tspan)
    sol = solve(oprob, integrator,saveat=tsteps1)
   return sol

end

function fit_Cybernetic_models(data::Matrix{Float64}, # dataset first row times second row OD
    label_exp::String, #label of the experiment
    model::Any, # ode model to use
    param;
    set_of_equations_to_fit=nothing,
    integrator=Tsit5(), # selection of sciml integrator
    multistart=false,
    n_restart=50,
    optimizer=BBO_adaptive_de_rand_1_bin_radiuslimited(),
    auto_diff_method=nothing,
    cons=nothing,
    opt_params...)

    max_t = data[1, end]
    min_t = data[1, 1]
    tspan = (min_t, max_t)
    tsteps = data[1, :]


    # define ic 
    Start_IC  =  [model.Bio_mass_conc, model.Substrate_concentrations..., model.Protein_concentrations...]

    # Defining ODE problem

    ODE_function = Generate_Generalized_cybernetic_ODEs_problem(model)
    ODE_prob = ODEProblem(ODE_function, Start_IC, tspan)
    # to do
    loss_function = define_loss_function_cybernetic_ODEs(data, set_of_equations_to_fit, ODE_prob, integrator, tsteps)

    res = KinbiontSolve(loss_function,
        Start_IC,
        param;
        opt=optimizer,
        auto_diff_method=auto_diff_method,
        multistart=multistart,
        n_restart=n_restart,
        cons=cons,
        opt_params...
    )

    remade_solution = solve(remake(ODE_prob, p=res.u), integrator, saveat=tsteps)

    loss_value = res.objective
    # to do
    res_param = vectorize_df_results_Cybernetic_model(label_exp, model, res.u, loss_value)
    # to do
    Kinbiont_res_Cybernetic_model = ("Cybernetic_model", res_param, remade_solution)

    return Kinbiont_res_Cybernetic_model

end