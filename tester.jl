using Kinbiont
using DifferentialEquations
a =  ODEs_system_sim(
    "SIR", #string of the model
   [0.9,0.1,0.0], # starting condition
    0.0, # start time of the sim
    30.0, # final time of the sim
    0.10, # delta t for poisson approx
    [0.5,0.3]; # parameters of the ODE model
)





using Plots
scatter(a)
sol_time = reduce(hcat, a.t)

sol_t = reduce(hcat, a.u)

data = vcat(sol_time,sol_t)


Start_IC = [0.9,0.1,0.0]



aa = fit_ODEs_System(data,
                    "test",
                    "SIR", 
                    [0.1,0.5],
                    Start_IC;
                     set_of_equations_to_fit = [1]
)




plot!(aa[3])




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
model = Kinbiont_Cybernetic_Model_sim(
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

using DifferentialEquations
oprob = ODEProblem(prob, u0, tspan)

# Solve the problem
sol = solve(oprob, Tsit5())

# Visualize the solution
using Plots
plot(sol)



model = Kinbiont_Cybernetic_Model_sim(
    Bio_mass_conc = 1.01,
    Substrate_concentrations = [2.0, 5.0],
    Protein_concentrations = [0.0, 0.0],
 #   rule =  [1.0, 0.1],
   rule = dynamic_allocation_rule_2,# sequential_allocation_rule,    # Use the dynamic allocation rule
    reaction = nothing,
    a = [0.1, 0.4],
    b = [0.00001, 0.000001],
    V_S = [0.1, 0.4],
    k_S = [0.1, 0.11],
    Y_S = [0.07, 0.11]
)

prob = Kinbiont_Cybernetic_Model_simulation(model, 0.0, 80.0, 0.01)

plot(reduce(hcat,prob.u)[1,:] )
plot(log.(reduce(hcat,prob.u)[1,:] ))

plot(reduce(hcat,prob.u)[2,:] )
plot(reduce(hcat,prob.u)[3,:] )
plot(reduce(hcat,prob.u)[4,:] )
plot(reduce(hcat,prob.u)[5,:] )
