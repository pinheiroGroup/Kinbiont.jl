using Kinbiont
using DifferentialEquations
using OptimizationBBO
using NaNMath
using Plots
using Distributions
using Optimization

# -------------------------------------------------------------------------------
# Cybernetic Models in Kinbiont
# -------------------------------------------------------------------------------
# This script demonstrates the following:
# 1. Setting up a **Kinbiont Cybernetic Model** with predefined parameters.
# 2. Running a simulation of biomass and substrate dynamics.
# 3. Fitting experimental data to the model using **Kinbiont's optimization tools**.
# -------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Step 1: Define a Cybernetic Model with Known Parameters
# ------------------------------------------------------------------------------

# Define the Kinbiont Cybernetic Model with specific parameters
model = Kinbiont_Cybernetic_Model(
    Bio_mass_conc = 1.01,  # Initial biomass concentration
    Substrate_concentrations = [5.0, 5.0],  # Concentrations of 2 substrates
    Protein_concentrations = [0.0, 0.0],  # Initial protein concentrations
    allocation_rule = threshold_switching_rule,  # Dynamic resource allocation rule
    reaction = nothing,  # No specific reaction function provided
    cost = nothing,  # No cost function
    protein_thresholds = 0.01,  # Protein activation threshold
    a = [0.1, 0.4],  # Synthesis rate constants for proteins
    b = [0.00001, 0.000001],  # Degradation constants for proteins
    V_S = [0.7, 0.1],  # Substrate utilization rates
    k_S = [0.1, 0.11],  # Saturation constants for substrates
    Y_S = [0.07, 0.11]  # Yield coefficients for biomass per substrate
)

# ------------------------------------------------------------------------------
# Step 2: Run a Simulation of the Cybernetic Model
# ------------------------------------------------------------------------------
# Simulate the model over time (from t=0 to t=100) with a time step of 0.1.
# We use the Tsit5 solver for integration.
simulation = Kinbiont_Cybernetic_Model_simulation(model, 0.0, 100.0, 0.1; Integration_method = Tsit5())

# Plot the results of the simulation
plot(simulation)

# ------------------------------------------------------------------------------
# Step 3: Prepare Data for Model Fitting
# ------------------------------------------------------------------------------
# Extract the data (biomass, substrate, and protein concentrations over time) from the simulation.
# Prepare the data in the format expected for model fitting.
data_to_fit = hcat(prob.t, reduce(hcat, prob.u)[1,:])
data_to_fit = hcat(data_to_fit, reduce(hcat, prob.u)[2,:])
data_to_fit = hcat(data_to_fit, reduce(hcat, prob.u)[3,:])
data_to_fit = hcat(data_to_fit, reduce(hcat, prob.u)[4,:])
data_to_fit = hcat(data_to_fit, reduce(hcat, prob.u)[5,:])
data_to_fit = permutedims(data_to_fit)  # Convert data to column-major order

# ------------------------------------------------------------------------------
# Step 4: Define a Cybernetic Model with Unknown Parameters for Fitting
# ------------------------------------------------------------------------------
# Set up a new Kinbiont model where certain parameters (a, V_S) are unknown and will be fitted.
model_fit = Kinbiont_Cybernetic_Model(
    Bio_mass_conc = 1.01,  # Initial biomass concentration
    Substrate_concentrations = [2.0, 5.0],  # Initial substrate concentrations
    Protein_concentrations = [0.0, 0.0],  # No initial protein concentrations
    allocation_rule = proportional_allocation_rule,  # Another dynamic allocation rule
    reaction = nothing,  # No specific reaction function
    cost = nothing,  # No cost function
    protein_thresholds = 0.01,  # Protein activation threshold
    a = [nothing, 0.1],  # One synthesis rate is unknown (to be fitted)
    b = [0.00001, 0.000001],  # Known degradation constants
    V_S = [nothing, 0.4],  # One substrate utilization rate is unknown (to be fitted)
    k_S = [0.1, 0.11],  # Saturation constants
    Y_S = [0.07, 0.11]  # Yield coefficients
)

# ------------------------------------------------------------------------------
# Step 5: Fit the Cybernetic Model to Experimental Data
# ------------------------------------------------------------------------------
# Use the `fit_Cybernetic_models` function to fit the model parameters to experimental data.
# The data_to_fit contains the time series of biomasses and other parameters for model fitting.
results = fit_Cybernetic_models(
    data_to_fit,  # Experimental data to fit
    "test",  # Name of the dataset for reference
    model_fit,  # Cybernetic model with unknown parameters to fit
    [0.01, 0.1];  # Initial guesses for the unknown parameters (a and V_S)
    set_of_equations_to_fit = nothing  # No additional equations to constrain fitting
)

# -------------------------------------------------------------------------------
# Notes:
# - **Kinbiont** models self-regulation, resource allocation, and optimization.
# - The **threshold_switching_rule** and **proportional_allocation_rule** allow the model to adjust resource allocation dynamically during simulation and fitting.
# - **Fitting** parameters to data can provide insights into how biological systems behave under different conditions.
# -------------------------------------------------------------------------------
