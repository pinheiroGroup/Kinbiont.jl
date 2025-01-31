
using Kinbiont
using DifferentialEquations
using OptimizationBBO
using NaNMath
using Plots
using Distributions
using Optimization




# Example of setting up the model and specifying known/unknown parameters
model = Kinbiont_Cybernetic_Model(
    Bio_mass_conc = 0.9, 
    Substrate_concentrations = [100.0, 5.0],  # Example concentrations of 3 substrates
    Protein_concentrations = [0.2, 0.1],  # Example concentrations for 3 proteins
    allocation_rule =  [1.0, 0.0],
    reaction = nothing,  # No specific reaction function in this case
    a = [0.5, 0.4],  # Synthesis rate constants
    b = [0.1, 0.1],  # Degradation constants for proteins
    V_S = [0.5, 0.4],  # Substrate utilization rate
    k_S = [0.1, 0.11],  # Saturation constants for substrates
    Y_S = [0.7, 0.7],  # Yield coefficients for cell mass per substrate
    cost = nothing,
    protein_thresholds =nothing
)

# Define the kimbiont data struct


model = Kinbiont_Cybernetic_Model(
    Bio_mass_conc = 1.01,
    Substrate_concentrations = [5.0, 5.0],
    Protein_concentrations = [0.0, 0.0],
 #   rule =  [1.0, 0.1],
 allocation_rule = threshold_switching_rule,    # Use the dynamic allocation rule
    reaction = nothing,
    cost = nothing,
    protein_thresholds = 0.01,
    a = [0.1, 0.4],
    b = [0.00001, 0.000001],
    V_S = [0.7, 0.1],
    k_S = [0.1, 0.11],
    Y_S = [0.07, 0.11]
)


simulation = Kinbiont_Cybernetic_Model_simulation(model,0.0,100.0,0.1; Integration_method = Tsit5())

plot(simulation)


plot(reduce(hcat,simulation.u)[1,:] )


data_to_fit = hcat(prob.t,reduce(hcat,prob.u)[1,:])
data_to_fit = hcat(data_to_fit,reduce(hcat,prob.u)[2,:])
data_to_fit = hcat(data_to_fit,reduce(hcat,prob.u)[3,:])
data_to_fit = hcat(data_to_fit,reduce(hcat,prob.u)[4,:])
data_to_fit = hcat(data_to_fit,reduce(hcat,prob.u)[5,:])
data_to_fit = permutedims(data_to_fit)

model_fit = Kinbiont_Cybernetic_Model(
    Bio_mass_conc = 1.01,
    Substrate_concentrations = [2.0, 5.0],
    Protein_concentrations = [0.0, 0.0],
 #   rule =  [1.0, 0.1],
 allocation_rule = proportional_allocation_rule,    # Use the dynamic allocation rule
    reaction = nothing,
    cost = nothing,
    protein_thresholds = 0.01,
    a = [nothing, 0.1],
    b = [0.00001, 0.000001],
    V_S = [nothing, 0.4],
    k_S = [0.1, 0.11],
    Y_S = [0.07, 0.11]
)

results = fit_Cybernetic_models(data_to_fit, 
"test", 
model_fit, 
[0.01,0.1]; 
set_of_equations_to_fit=nothing
)