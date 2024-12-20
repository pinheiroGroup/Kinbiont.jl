# Define a generalized Cybernetic Model structure
struct Kinbiont_Cybernetic_Model
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


# Custom constructor to allow initialization using keyword arguments
function Kinbiont_Cybernetic_Model(; Bio_mass_conc, Substrate_concentrations, Protein_concentrations,
    allocation_rule, reaction, a, b, V_S, k_S, Y_S,cost,protein_thresholds)
    return Kinbiont_Cybernetic_Model(Bio_mass_conc, Substrate_concentrations, Protein_concentrations,
    allocation_rule, reaction, a, b, V_S, k_S, Y_S,cost,protein_thresholds)
end

function loss_RE_Cybernetic(data, index_of_eqs, index_of_data, ODE_prob, integrator, p, tsteps)
    # Use the param 'p' as input
    sol = solve(
        ODE_prob,
        integrator,
        p=p,  # Pass param here
        saveat=tsteps,
        verbose=false,
        abstol=1e-10,
        reltol=1e-10,
    )
    # println(p)

    sol_t = reduce(hcat, sol.u)

  #  println(size(sol_t)[2])

    if size(sol_t)[2] == size(data)[2]
        lossa = 0.5 .* abs2.(1.0 .- (data[index_of_data, :] ./ sol_t[index_of_eqs, :]))
                # index 0/0
                lossa[ isnan.(lossa)] .= 0.0   
                lossa = NaNMath.sum(lossa)/ length(data[2, :])
    else
        lossa = 10.0^9 * length(data[2, :])
    end

    return lossa, sol
end

function define_loss_function_cybernetic_ODEs(data, set_of_equation_to_fit, ODE_prob, integrator, tsteps)


    if !isnothing(set_of_equation_to_fit)
        index_of_eqs = set_of_equation_to_fit 
        index_of_data = set_of_equation_to_fit .+1 


    else

        index_of_eqs = 1:1:(size(data)[1]-1)
        index_of_data = index_of_eqs .+ 1

    end



    function loss_Cybernetic(data,ODE_prob,p, tsteps, index_of_eqs, index_of_data,  integrator )
        # Use the param 'p' as input
        sol = solve(
            ODE_prob,
            integrator,
            p=p,  # Pass param here
            saveat=tsteps,
            verbose=false,
            abstol=1e-10,
            reltol=1e-10,
        )
      #   println(p)
        
        sol_1 = reduce(hcat, sol.u)

      #  println(size(sol_t)[2])

        if size(sol_1)[2] == size(data)[2]
            lossa = NaNMath.sum( abs2.(data[index_of_data, :] .- sol_1[index_of_eqs, :]))/ length(data[2, :])
                    # index 0/0
        else
            lossa = 10.0^9 * length(data[2, :])
        end
    
        return lossa, sol
    end
    



    return (p) -> loss_Cybernetic(data,ODE_prob,p, tsteps, index_of_eqs, index_of_data,  integrator )

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

function replace_nothing_with_param(Cybernetic_Model1, param)
    param_counter = 0

    # Helper function to replace `nothing` with "param[$counter]" and track counter
    function replace_in_array(arr)
        return [isnothing(val) ? param[(param_counter += 1)] : val for val in arr]
    end

    # Replace `nothing` in each parameter array
    a_updated = replace_in_array(Cybernetic_Model1.a)
    b_updated = replace_in_array(Cybernetic_Model1.b)
    V_S_updated = replace_in_array(Cybernetic_Model1.V_S)
    k_S_updated = replace_in_array(Cybernetic_Model1.k_S)
    Y_S_updated = replace_in_array(Cybernetic_Model1.Y_S)

    return (a_updated, b_updated, V_S_updated, k_S_updated, Y_S_updated)
end


    # Function to get value or substitute param if needed
function get_val(meta_param)
        if isa(meta_param, String)
            # Parse and evaluate the expression to get a numeric value
            return Meta.eval(Meta.parse(meta_param))
        end
        return meta_param  # If not a string, return the value as it is
end







function Cybernetic_Model_odes(du, u,param, t,  Cybernetic_Model1)
        n = length(Cybernetic_Model1.Substrate_concentrations)  # Number of substrates
        a, b, V_S, k_S, Y_S = replace_nothing_with_param(Cybernetic_Model1, param)
     #   println("Using params: a = ", a, ", b = ", b, ", V_S = ", V_S, ", k_S = ", k_S, ", Y_S = ", Y_S)

        # Extract substrates and proteins
        S = u[2:(n+1)]  # Substrate concentrations (S1, S2, ..., Sn)
        P = u[(n+2):end]  # Key protein concentrations (P1, P2, ..., Pn)
        
        rule = Cybernetic_Model1.allocation_rule  # Allocation rule for synthesis rate
        protein_thresholds = Cybernetic_Model1.protein_thresholds
        cost = Cybernetic_Model1.cost
        
        # If the allocation rule is a function, use it
        if typeof(rule) != Vector{Float64}
            alloc = rule(a, b, V_S, k_S, Y_S, P, S, cost, protein_thresholds)
        else
            alloc = rule ./ sum(rule)  # Normalize if it's a vector
        end
    
        # Rate of change for substrates and proteins
        for i in 1:n
            du[i + 1] = -V_S[i] * P[i] * S[i] / (k_S[i] + S[i]) * u[1]  # Substrate consumption
            du[i + n + 1] = a[i] * alloc[i] * k_S[i] - b[i] * P[i] * u[1]  # Protein synthesis rate
        end
    
        # Biomass growth dynamics
        du[1] = -sum(Y_S .* du[2:(n+1)]) * u[1]
end
    
# Function to generate generalized cybernetic ODE system with parameter replacement
function Generate_Generalized_cybernetic_ODEs_problem(Cybernetic_Model1)
    return (du, u,param, t) -> Cybernetic_Model_odes(du, u, param, t,  Cybernetic_Model1)
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









  function Cybernetic_Model_odes_sim(du,u,t,param)

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


    return Cybernetic_Model_odes_sim  # In-place modification of the rate of change vector

end


function Kinbiont_Cybernetic_Model_simulation(model,tmin,tmax,deltaT; integrator = Tsit5())

    tspan =(tmin,tmax)
    tsteps = tmin:deltaT:tmax
    tsteps1 = [tsteps[i] for i in 1:length(tsteps)]
    u0 = [model.Bio_mass_conc, model.Substrate_concentrations..., model.Protein_concentrations...]
    Cybernetic_Model_odes_sim = Generate_Generalized_cybernetic_ODEs_problem_simulation(model)
    oprob = ODEProblem(Cybernetic_Model_odes_sim, u0, tspan)
    sol = solve(oprob, integrator,saveat=tsteps1)
   return sol

end





function KinbiontSolve_test(loss_function, u0, param_g;
    opt,
    auto_diff_method=nothing,
    multistart=false,
    n_restart=50,
    cons=nothing,
    opt_params...)
    
    # Ensure param is passed to the optimization function properly
    if isnothing(auto_diff_method) && isnothing(cons)
        # Ensure that loss_function takes `p` (the parameters) and `param`
        optf = Optimization.OptimizationFunction((x,p) -> loss_function(x))
    elseif !isnothing(auto_diff_method) && isnothing(cons)
        optf = Optimization.OptimizationFunction((x,p) -> loss_function(x), auto_diff_method)
    elseif isnothing(auto_diff_method) && !isnothing(cons)
        optf = Optimization.OptimizationFunction((x,p) -> loss_function(x), cons=cons)
    else
        optf = Optimization.OptimizationFunction((x,p) -> loss_function(x), auto_diff_method, cons=cons)
    end

    # Set bounds for optimization if needed
    opt_params = check_bounds_opt(opt, param_g, opt_params...)

    # Define the optimization problem
    prob = OptimizationProblem(optf, param_g, u0; opt_params...)

    # Solve the optimization problem
    if multistart
        sol = solve(prob, MultistartOptimization.TikTak(n_restart), opt)
    else
        sol = solve(prob, opt)
    end

    return sol
end

function check_bounds_opt(opt,p_guess,
    opt_params...)
 
    if SciMLBase.requiresbounds(opt)

        opt_params = Dict(opt_params)




        if  !(:lb  in keys(opt_params)) && !(:ub  in keys(opt_params) )

            @warn "The used optimization method requires box bounds, Kinbiont.jl will use upper bounds that are 10 times the guess
             and lower bounds that are 10 times lower the guess.
             This choice can be suboptimal. Note that the Kinbiont.jl default optimizer requires box bounds to guaranteed a Real N(t) and positive parameters.
             To avoid to specify the bounds use you can use an optimizer that do not require it, e.g., `optimizer = NOMADOpt()`. 
             Note that numerical instabilities may occur.
             "
            opt_params = (opt_params..., lb= p_guess ./10 , ub = p_guess .*10 )

        end
    
    end
    
    return opt_params

end


function fit_Cybernetic_models(data::Matrix{Float64}, 
    label_exp::String, 
    model::Any, 
    param_g::Vector{Float64}; 
    set_of_equations_to_fit=nothing,
    integrator=Tsit5(), 
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

    # Define initial conditions (Start_IC)
    Start_IC = [model.Bio_mass_conc, model.Substrate_concentrations..., model.Protein_concentrations...]

    ODE_function = Generate_Generalized_cybernetic_ODEs_problem(model)
    ODE_prob = ODEProblem(ODE_function, Start_IC, tspan)
    # Define the loss function
    loss_function = define_loss_function_cybernetic_ODEs(data, set_of_equations_to_fit, ODE_prob, integrator, tsteps)

 
    res = KinbiontSolve_test(loss_function, Start_IC, param_g; 
    opt=optimizer, 
    auto_diff_method=auto_diff_method, 
    multistart=multistart, 
    n_restart=n_restart, 
    cons=cons, 
    opt_params...)

   return res

end    

