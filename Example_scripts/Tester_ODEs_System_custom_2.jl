using Kinbiont
using DifferentialEquations
using CSV
using SymbolicRegression
using Plots
using StatsBase
using Distributions

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# The system described below represents an interaction between four variables:
# - **u1**: A reactant influenced by u4
# - **u2**: A product formed via an intermediate reaction
# - **u3**: Another intermediate compound
# - **u4**: Decreases as it drives the reaction forward
#
# Our goal is to examples how to fit the ODE system that is user defined 
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Step 1: Define the ODE System with an Unknown Functional Form
# ------------------------------------------------------------------------------
function model_1(du, u, param, t)
    # State variables:
    # u1 -> Reactant
    # u2 -> Intermediate product
    # u3 -> Final product
    # u4 -> Driving factor (e.g., enzyme, catalyst, or resource pool)

    # Parameter descriptions:
    # param[1] -> Rate of reactant conversion, dependent on u4
    # param[2] -> Degradation rate of intermediates
    # param[3] -> Transition rate from u2 to u3
    # param[4] -> Scaling factor for the first reaction

    # Differential equations:
    du[1] = param[1] * u[1] * u[4]                  # Reactant conversion
    du[2] = param[4] * du[1] - param[3] * u[2] - param[2] * u[2]  # Intermediate balance
    du[3] = param[3] * u[2] - param[2] * u[3]       # Final product formation
    du[4] = -du[1]                                  # Resource consumption
end

# ------------------------------------------------------------------------------
# Step 2: Set Initial Conditions and Parameter Ranges
# ------------------------------------------------------------------------------
u0 = [0.1, 0.0, 0.0, 1.0]  # Initial conditions for [u1, u2, u3, u4]

# True parameter values (used for simulation)
param = [0.1, 0.01, 0.5, 0.42]

# Define lower and upper bounds for parameter estimation
lb1 = [0.01, 0.0001, 0.0, 0.01]  # Lower bounds
ub1 = [0.2, 0.3, 1.1, 1.0]       # Upper bounds

# Initial parameter guess (midpoint between bounds)
param_guess = lb1 .+ (ub1 .- lb1) ./ 2  

# Noise level
noise_value = 0.01

# ------------------------------------------------------------------------------
# Step 3: Simulate the ODE System
# ------------------------------------------------------------------------------
Simulation = ODEs_system_sim(
    model_1, # ODE function
    u0,      # Initial conditions
    0.0,     # Start time
    50.0,    # End time
    1.0,     # Time step for Poisson approximation
    param    # True parameters
)

# ------------------------------------------------------------------------------
# Step 4: Add Noise to Simulated Data
# ------------------------------------------------------------------------------
sol_time = reduce(hcat, Simulation.t) # Extract time points
sol_t = reduce(hcat, Simulation.u)    # Extract simulation data

# Adding uniform random noise to mimic experimental uncertainty
noise_uniform = rand(Uniform(-0.05, 0.05), size(sol_t)[2])
sol_t_noise = [sol_t[i, :] .+ rand(Uniform(-0.05, 0.05), size(sol_t)[2]) for i in 1:size(sol_t)[1]]
sol_t_noise = permutedims(reduce(hcat, sol_t_noise))

# Combine noisy data with time points
data = vcat(sol_time, sol_t_noise)

# ------------------------------------------------------------------------------
# Step 5: Visualize Noisy Data
# ------------------------------------------------------------------------------
display(scatter(data[1, :], data[2, :], label="u1"))
display(scatter!(data[1, :], data[3, :], label="u2"))
display(scatter!(data[1, :], data[4, :], label="u3"))
display(scatter!(data[1, :], data[5, :], label="u4"))

# ------------------------------------------------------------------------------
# Step 6: Fit the Model Using Kinbiont
# ------------------------------------------------------------------------------
fit = fit_ODEs_System(
    data,
    "test",     # Label for dataset
    model_1,    # ODE model
    param_guess, # Initial parameter guess
    u0;         # Initial conditions
    lb=lb1,     # Lower bounds
    ub=ub1      # Upper bounds
)

# ------------------------------------------------------------------------------
# Step 7: Plot the Fitted Model
# ------------------------------------------------------------------------------
plot!(fit[3], label="Fitted Model")
