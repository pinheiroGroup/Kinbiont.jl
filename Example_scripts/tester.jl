using Kinbiont
using DifferentialEquations
using OptimizationBBO
using Plots

using Catalyst
using DiffEqParamEstim



u0 = [:S => 301, :E => 100, :SE => 0, :P => 0]
ps = [:kB => 0.00166, :kD => 0.0001, :kP => 0.1]

model_Michaelis_Menten_enzyme_kinetics= @reaction_network begin
    kB, S + E --> SE
    kD, SE --> S + E
    kP, SE --> P + E
end






b=  Kinbiont_Reaction_network_sim("Michaelis_Menten",
    u0,
   0.0, # start time of the sim
   10.0, # final time of the sim
    0.1, # delta t for poisson approx
    ps; # parameters of the ODE model
    )

plot(b)




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

# Define initial conditions (biomass, substrates, proteins) for the model

# Time span for solving the ODE
tspan = (0.0, 10.0)
prob = Kinbiont_Cybernetic_Model_simulation(model, 0.0, 100.0, 0.1)
# Solve the ODE system using the DifferentialEquations package
plot(prob)
plot(reduce(hcat,prob.u)[1,:] )
plot(reduce(hcat,prob.u)[2,:] )
plot(reduce(hcat,prob.u)[3,:] )
plot(reduce(hcat,prob.u)[4,:] )
plot(reduce(hcat,prob.u)[5,:] )





model_1 = Kinbiont_Cybernetic_Model(
    Bio_mass_conc = 1.01,
    Substrate_concentrations = [2.0, 5.0],
    Protein_concentrations = [0.0, 0.0],
 #   rule =  [1.0, 0.1],
 allocation_rule = threshold_switching_rule,# sequential_allocation_rule,    # Use the dynamic allocation rule
    reaction = nothing,
    cost = nothing,
    protein_thresholds = nothing,
    a = [0.1, 0.4],
    b = [0.00001, 0.000001],
    V_S = [0.1, 0.4],
    k_S = [0.1, 10.1],
    Y_S = [0.07, 0.11]
)

probd = Kinbiont_Cybernetic_Model_simulation(model, 0.0, 80.0, 0.01)

plot(reduce(hcat,probd.u)[1,:] )
plot(log.(reduce(hcat,probd.u)[1,:] ))

plot(reduce(hcat,probd.u)[2,:] )
plot(reduce(hcat,probd.u)[3,:] )
plot(reduce(hcat,probd.u)[4,:] )
plot(reduce(hcat,probd.u)[5,:] )

data_cybernetic  = hcat(reduce(hcat,probd.t)[1,:],reduce(hcat,probd.u)[1,:])
data_cybernetic  = hcat(data_cybernetic,reduce(hcat,probd.u)[2,:] )
data_cybernetic  = hcat(data_cybernetic,reduce(hcat,probd.u)[3,:])
data_cybernetic  = hcat(data_cybernetic,reduce(hcat,probd.u)[4,:])
data_cybernetic  = hcat(data_cybernetic,reduce(hcat,probd.u)[5,:])

data_cybernetic = permutedims(data_cybernetic)


model_to_fit =  Kinbiont_Cybernetic_Model(
    Bio_mass_conc = 1.01,
    Substrate_concentrations = [2.0, 5.0],
    Protein_concentrations = [0.0, 0.0],
 #   rule =  [1.0, 0.1],
 allocation_rule = threshold_switching_rule,# sequential_allocation_rule,    # Use the dynamic allocation rule
    reaction = nothing,
    cost = nothing,
    protein_thresholds = nothing,
    a = [0.1, 0.4],
    b = [0.00001, 0.000001],
    V_S = [nothing, 0.4],
    k_S = [0.1, nothing],
    Y_S = [0.07, 0.11]
)


param_g = [0.2]
fit_CM = fit_Cybernetic_models(data_cybernetic, # dataset first row times second row OD
    "test", #label of the experiment
    model_to_fit, # ode model to use
    param_g;
    set_of_equations_to_fit=nothing,
    )

