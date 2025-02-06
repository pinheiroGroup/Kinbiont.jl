using Kinbiont
using DifferentialEquations
using OptimizationBBO
using NaNMath
using Plots
using Distributions

# ----------------------------------------------
# Step 1: Generate Simulated Data with the SIR Model
# ----------------------------------------------

Simulation = ODEs_system_sim(
    "SIR",          # Name of the model (Susceptible-Infected-Recovered)
    [0.9, 0.1, 0.0], # Initial conditions: [S, I, R]
    0.0,            # Start time of the simulation
    30.0,           # End time of the simulation
    1.0,            # Time step for Poisson approximation
    [0.5, 0.3]      # Parameters of the ODE model: [Infection rate, Recovery rate]
)

# ----------------------------------------------
# Step 2: Plot the Simulated Data
# ----------------------------------------------

scatter(Simulation)  # Scatter plot of the simulation results

# ----------------------------------------------
# Step 3: Add Noise and Format Data for Kinbiont
# ----------------------------------------------

# Extracting time points and solution values
sol_time = reduce(hcat, Simulation.t) 
sol_t = reduce(hcat, Simulation.u)

# Adding uniform random noise to the simulation data
noise_uniform = rand(Uniform(-0.05, 0.05), size(sol_t)[2])
sol_t_noise = [sol_t[i, :] .+ rand(Uniform(-0.05, 0.05), size(sol_t)[2]) for i in 1:size(sol_t)[1]]
sol_t_noise = permutedims(reduce(hcat, sol_t_noise))

# Combining time and noisy solution data
data = vcat(sol_time, sol_t_noise)

# Plot the noisy data
scatter(data[1, :], data[2, :])  # Susceptible
scatter!(data[1, :], data[3, :]) # Infected
scatter!(data[1, :], data[4, :]) # Recovered

# ----------------------------------------------
# Step 4: Fit the Full Dataset to the SIR Model
# ----------------------------------------------

Start_IC = [0.9, 0.1, 0.0]  # Initial conditions

fit = fit_ODEs_System(
    data,
    "test",     # Label for the dataset
    "SIR",      # Model name
    [0.1, 0.5], # Initial guess for parameters
    Start_IC    # Initial conditions
)

plot!(fit[3])  # Plot the fitted results

# ----------------------------------------------
# Step 5: Remove "Recovered" (R) Data and Refit
# ----------------------------------------------

# Keep only time, S (Susceptible), and I (Infected)
data_reduced = hcat(data[1, :], data[2, :])
data_reduced = permutedims(hcat(data_reduced, data[3, :]))

# Fit using only the first two states (S and I)
fit = fit_ODEs_System(
    data_reduced,
    "test",
    "SIR",
    [0.1, 0.5], # Initial guess for parameters
    Start_IC;
    set_of_equation_to_fit = [1, 2]  # Only fit S and I equations
)

# ----------------------------------------------
# Step 6: Plot the Data and New Fit
# ----------------------------------------------

scatter(data[1, :], data[2, :])  # Susceptible
scatter!(data[1, :], data[3, :]) # Infected
scatter!(data[1, :], data[4, :]) # Recovered
plot!(fit[3])  # Plot the fitted model without R measurements
