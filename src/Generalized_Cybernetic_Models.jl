# Define a generalized Cybernetic Model structure
struct Kinbiont_Cybernetic_Model_sim
    Bio_mass_conc::Float64             # Biomass concentration
    Substrate_concentrations::Vector{Float64}  # Array of substrate concentrations (S1, S2, ..., Sn)
    Protein_concentrations::Vector{Float64}    # Array of protein concentrations (P1, P2, ..., Pn)
    allocation_rule::Any                      # Function for the functional dependence of synthesis rate
    reaction::Union{Function, Nothing}  # Reaction function (optional, depends on the model)
    a::Vector{Float64}                 # Synthesis rate constants for proteins
    b::Vector{Float64}                 # Degradation constants for proteins
    V_S::Vector{Float64}               # Substrate utilization rate
    k_S::Vector{Float64}               # Saturation constants for substrates
    Y_S::Vector{Float64}               # Yield coefficients for cell mass per substrate
end

# Custom constructor to allow initialization using keyword arguments
function Kinbiont_Cybernetic_Model_sim(; Bio_mass_conc, Substrate_concentrations, Protein_concentrations,
    rule, reaction, a, b, V_S, k_S, Y_S)
    return Kinbiont_Cybernetic_Model_sim(Bio_mass_conc, Substrate_concentrations, Protein_concentrations,
     rule, reaction, a, b, V_S, k_S, Y_S)
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
    # Define allocation vector for resource allocation (normalized)
    if typeof(rule) == Vector{Float64} 
   
   
        alloc = rule./ sum(rule)  # Assuming proportional allocation based on substrate concentrations
    


    end

  function Cybernetic_Model_odes(du,u,t,param)

        # alias the state vector from the input (biomass, substrates, proteins)
        X_c = u[1]  # Biomass concentration (cells)
        S = u[2:(n+1)]  # Substrate concentrations (S1, S2, ..., Sn)
        P = u[(n+2):end]  # Key protein concentrations (P1, P2, ..., Pn)
        alloc = rule(S,k_S)  # Assuming proportional allocation based on substrate concentrations

    # Protein synthesis rate (using the functional form `rule`)
     for i in 1:n
           du[i + n + 1]  = a[i] * alloc[i] * k_S[i] - b[i] * P[i]  # Rate of change of each key protein
      end

    # Substrate consumption using Monod kinetics
       for i in 1:n
         du[i + 1]= -V_S[i] * P[i] * S[i] / (k_S[i] + S[i])  # Rate of change of each substrate
        end

    # Cell growth dynamics (cumulative effect of substrate utilization)
            du[1] = sum(Y_S .* du[(n+2):end])  # Change in biomass concentration based on substrate utilization
    end


    return Cybernetic_Model_odes  # In-place modification of the rate of change vector

end


function Kinbiont_Cybernetic_Model_simulation(model,tmin,tmax,deltaT; integrator = Tsit5())

    tspan =(tmin,tmax)
    tsteps = tmin:deltaT:tmin
    u0 = [model.Bio_mass_conc, model.Substrate_concentrations..., model.Protein_concentrations...]
    Cybernetic_Model_odes = Generate_Generalized_cybernetic_ODEs_problem_simulation(model)
    oprob = ODEProblem(Cybernetic_Model_odes, u0, tspan)
    prob = Generate_Generalized_cybernetic_ODEs_problem_simulation(model)

    sol = solve(oprob, integrator)
   return sol

end

