using Kinbiont
using DifferentialEquations
using OptimizationBBO
using NaNMath
using Plots
using Distributions

# ------------------------------------------------------------------------------
# Macromolecular Crowding and Enzyme Aggregation Model
# ------------------------------------------------------------------------------
# Background:
# The cytosol of cells is highly crowded with macromolecules, occupying 20-44% 
# of the total cellular volume. This crowding increases with the specific 
# growth rate (SGR) and affects metabolic rates due to:
# 1. Increased viscosity, restricting diffusion.
# 2. Shifts in biochemical equilibria, favoring macromolecular aggregation.
#
# One key unresolved aspect is how enzyme aggregation impacts metabolic dynamics. 
# To address this, we model enzyme aggregation as a reversible process:
# - Free enzymes (E) bind together to form inactive aggregates (En).
# - Aggregates can dissociate back into active enzymes with the help of chaperones.
# - This process directly affects the available enzymatic activity in the cell.
#
# Below, we implement a **simple kinetic mechanism** for this aggregation process.
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Step 1: Define the Enzyme Aggregation Model (ODE System)
# ------------------------------------------------------------------------------
function enzyme_aggregation(du, u, param, t)
    # Unpacking state variables:
    e, x, y, m = u  
    # e  -> Free active enzymes
    # x  -> Enzyme-substrate complex
    # y  -> Aggregates (inactive enzyme forms)
    # m  -> Substrate concentration

    # Unpacking parameters:
    k1, k2, k3, k4, k_cat, n = param  
    # k1, k2 -> Aggregation/dissociation rates
    # k3, k4 -> Substrate binding/release rates
    # k_cat  -> Catalytic rate
    # n      -> Hill coefficient (cooperativity in aggregation)

    # ODE System:
    du[1] = k4 * x - k3 * m * e + k2 * y^n - k1 * e + k_cat * x  # Free enzyme balance
    du[2] = k3 * m * e - k4 * x - k_cat * x                      # Enzyme-substrate complex
    du[3] = k1 * e - k2 * y^n                                    # Inactive enzyme aggregates
    du[4] = -du[1]                                              # Substrate degradation rate
end

# ------------------------------------------------------------------------------
# Step 2: Define Initial Conditions and Parameters
# ------------------------------------------------------------------------------
u0 = [1.0, 0.1, 0.1, 1.0]  # Initial conditions: [e, x, y, m]

param = [0.1, 0.1, 0.05, 0.05, 0.02, 2]  
# Parameter list: [k1, k2, k3, k4, k_cat, n]

# ------------------------------------------------------------------------------
# Step 3: Run the Simulation
# ------------------------------------------------------------------------------
Simulation = ODEs_system_sim(
    enzyme_aggregation, # Custom ODE function
    u0,  # Initial conditions
    0.0, # Start time
    30.0, # End time
    1.0, # Time step for Poisson approximation
    param # Model parameters
)

# ------------------------------------------------------------------------------
# Step 4: Plot the Simulated Data
# ------------------------------------------------------------------------------
scatter(Simulation)  # Scatter plot of the simulation results

# ------------------------------------------------------------------------------
# Step 5: Add Noise and Format Data for Kinbiont
# ------------------------------------------------------------------------------
sol_time = reduce(hcat, Simulation.t) # Extract time points
sol_t = reduce(hcat, Simulation.u)    # Extract solution values

# Adding uniform random noise to the simulated data
noise_uniform = rand(Uniform(-0.05, 0.05), size(sol_t)[2])
sol_t_noise = [sol_t[i, :] .+ rand(Uniform(-0.05, 0.05), size(sol_t)[2]) for i in 1:size(sol_t)[1]]
sol_t_noise = permutedims(reduce(hcat, sol_t_noise))

# Combine time and noisy solution data into a single matrix
data = vcat(sol_time, sol_t_noise)

# ------------------------------------------------------------------------------
# Step 6: Plot Noisy Data
# ------------------------------------------------------------------------------
display(scatter(data[1, :], data[2, :]))  # Free enzymes
display(scatter!(data[1, :], data[3, :])) # Enzyme-substrate complex
display(scatter!(data[1, :], data[4, :])) # Aggregates
display(scatter!(data[1, :], data[5, :])) # Substrate

# ------------------------------------------------------------------------------
# Step 7: Fit the Model Using Kinbiont
# ------------------------------------------------------------------------------
fit = fit_ODEs_System(
    data,
    "test",  # Label for dataset
    enzyme_aggregation,  # Custom ODE function
    param,  # Initial parameter guess
    u0  # Initial conditions
)

# ------------------------------------------------------------------------------
# Step 8: Plot the Fitted Model
# ------------------------------------------------------------------------------
plot!(fit[3])  # Overlay fitted model on the data
