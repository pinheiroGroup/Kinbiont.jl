using Kinbiont
using DifferentialEquations
using OptimizationBBO
using Plots
using Catalyst
using DiffEqParamEstim
using Random
using Distributions

# --------------------------------------------------------------
# **Michaelis-Menten Reaction Network Simulation**
# --------------------------------------------------------------

# Initial conditions
u0 = [:S => 301, :E => 100, :SE => 0, :P => 0]
ps = [:kB => 0.00166, :kD => 0.0001, :kP => 0.1]

# Define Michaelis-Menten enzyme kinetics reaction network
model_Michaelis_Menten = @reaction_network begin
    kB, S + E --> SE
    kD, SE --> S + E
    kP, SE --> P + E
end 

# Run simulation
Simulation = Kinbiont.Kinbiont_Reaction_network_sim(
    "Michaelis_Menten",
    u0,
    0.0, 10.0, 0.1, # Start time, end time, step size
    ps
)

# Plot the simulation results
plot(Simulation)

# --------------------------------------------------------------
# **Generate Noisy Data for Fitting**
# --------------------------------------------------------------
sol_time = reduce(hcat, Simulation.t)
sol_t = reduce(hcat, Simulation.u)

# Add noise to the dataset
noise = rand(Uniform(-0.01, 0.05), size(sol_t))
sol_t_noise = sol_t .+ noise
data = vcat(sol_time, permutedims(sol_t_noise))

# Scatter plot of noisy data
scatter(data[1, :], data[2, :], label="S")
scatter!(data[1, :], data[3, :], label="E")
scatter!(data[1, :], data[4, :], label="SE")
scatter!(data[1, :], data[5, :], label="P")

# --------------------------------------------------------------
# **Fit Michaelis-Menten Model to Data**
# --------------------------------------------------------------
fit = RN_fit(data, model_Michaelis_Menten, u0, ps)
plot!(fit[4])  # Overlay fit result on the plot

# --------------------------------------------------------------
# **Glycolysis Reaction Network**
# --------------------------------------------------------------

model_Glycolysis = @reaction_network begin
    (kf1, kr1), Glucose + ATP <--> Glucose6P + ADP
    (kf2, kr2), Glucose6P <--> Fructose6P
    (kf3, kr3), Fructose6P + ATP <--> Fructose16BP + ADP
    (kf4, kr4), Fructose16BP <--> DHAP + GAP
    (kf5, kr5), DHAP <--> GAP
    (kf6, kr6), GAP + NADplus + Pi <--> BPG13 + NADH + Hplus
    (kf7, kr7), BPG13 + ADP <--> PG3 + ATP
    (kf8, kr8), PG3 <--> PG2
    (kf9, kr9), PG2 <--> PEP + H2O
    (kf10, kr10), PEP + ADP <--> Pyruvate + ATP
end

# Initial conditions
u0_glycolysis = [
    :Glucose => 100, :ATP => 200, :Glucose6P => 0, 
    :Fructose6P => 0, :Fructose16BP => 0, :DHAP => 0, 
    :GAP => 0, :NADplus => 100, :NADH => 0, :Pi => 100, 
    :BPG13 => 0, :PG3 => 0, :PG2 => 0, :PEP => 0, 
    :Pyruvate => 0, :ADP => 0, :H2O => 0, :Hplus => 0
]

# Parameters
ps_glycolysis = [
    :kf1 => 0.0011, :kr1 => 0.005, 
    :kf2 => 0.005, :kr2 => 0.002, 
    :kf3 => 0.02, :kr3 => 0.01, 
    :kf4 => 0.015, :kr4 => 0.007, 
    :kf5 => 0.01, :kr5 => 0.005, 
    :kf6 => 0.02, :kr6 => 0.01, 
    :kf7 => 0.03, :kr7 => 0.015, 
    :kf8 => 0.01, :kr8 => 0.005, 
    :kf9 => 0.008, :kr9 => 0.004, 
    :kf10 => 0.04, :kr10 => 0.02
]

# Simulate Glycolysis
Simulation = Kinbiont.Kinbiont_Reaction_network_sim(
    model_Glycolysis,
    u0_glycolysis,
    0.0, 30.0, 0.1, 
    ps_glycolysis
)

plot(Simulation)

# --------------------------------------------------------------
# **TCA Cycle Reaction Network**
# --------------------------------------------------------------
model_TCA_cycle = @reaction_network begin
    k1, Acetyl_CoA + Oxaloacetate + H2O --> Citrate
    k2, Citrate --> Isocitrate
    k3, Isocitrate + NAD_plus --> α_Ketoglutarate + NADH + CO2 + H_plus
    k4, α_Ketoglutarate + NAD_plus + H2O --> Succinyl_CoA + NADH + CO2 + H_plus
    k5, Succinyl_CoA --> Succinate
    k6, Succinate + NAD_plus --> Fumarate + NADH + H_plus
    k7, Fumarate + H2O --> Malate
    k8, Malate + NAD_plus --> Oxaloacetate + NADH + H_plus
end

# Initial conditions
u0_TCA = [
    :Acetyl_CoA => 1, :Citrate => 0, :Isocitrate => 0, :α_Ketoglutarate => 0,
    :Succinyl_CoA => 0, :Succinate => 0, :Fumarate => 0, :Malate => 0, :Oxaloacetate => 0,
    :H2O => 10, :NAD_plus => 10, :NADH => 0, :CO2 => 0, :H_plus => 0
]

# Parameters
ps_TCA = [
    :k1 => 1.0, :k2 => 1.0, :k3 => 1.0, :k4 => 1.0,
    :k5 => 1.0, :k6 => 1.0, :k7 => 1.0, :k8 => 1.0
]

# Simulate TCA Cycle
Simulation = Kinbiont.Kinbiont_Reaction_network_sim(
    model_TCA_cycle,
    u0_TCA,
    0.0, 30.0, 0.1,
    ps_TCA
)

plot(Simulation)
