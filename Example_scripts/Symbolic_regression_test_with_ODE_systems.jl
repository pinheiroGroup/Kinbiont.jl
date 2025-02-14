
using Kinbiont
using DifferentialEquations
using CSV
using SymbolicRegression
using Plots
using StatsBase
using Distributions
using OptimizationBBO
using Optimization
using NaNMath

# In this example, we simulate a system of ODEs with Kinbiont.
# One of the parameters in the system (param[4]) is affected by an experimental feature.
# The dependence of this parameter on the feature is **quadratic**: param[4] = param0 * (1 - feature)^2.
# The goal is to fit the ODE system to data generated under different feature values
# and then apply symbolic regression to recover the underlying relationship.


# Function defining the unknown quadratic response of the parameter to the feature
function unknown_response(feature)
    response = (1 - feature)^2
    return response
end

# Function defining time-dependent activity, used in the ODE system
function time_dep(time)
    activity = exp((- (time - 20)^2) / 20)
    return activity
end

# Defining the system of ODEs a tres state system. Where, the first state are bacteria reproduciing eating the nutrient (u[4), the production of the second state has a time dependent activity and the third state is produced by the second state with a constant maturation rate.
function model_1(du, u, param, t)
    # u[1]: Population 1
    # u[2]: Intermediate state 1
    # u[3]: Intermediate state 2
    # u[4]: Available resource
    
    du[1] = param[1] * u[1] * u[4]  # Growth equation
    du[2] = time_dep(t) * param[4] * u[1] - param[3] * u[2] - param[2] * u[2]  # Intermediate process
    du[3] = param[3] * u[2] - param[2] * u[3]  # Second intermediate process
    du[4] = -du[1]  # Resource depletion
end

# Initial conditions for the variables
u0 = [0.1, 0.0, 0.0, 1.0]

# Parameters: [growth rate, degradation rate, conversion rate, feature-dependent parameter]
param = [0.1, 0.001, 0.5, 0.42]
lb1 = [0.01, 0.0001, 0.1, 0.0]
ub1 = [0.2, 0.3, 1.1, 1.0]
param_guess = lb1 .+ (ub1 .- lb1) ./ 2

param0 = param[4]  # Store initial value of the feature-dependent parameter
noise_value = 0.01  # Noise level

# Defining the range of the perturbation on the feature
results_fit = Any
feature_range = 0.0:0.1:2.0

plot(0, 0)  # Empty plot for visualization

for f in feature_range
    # Quadratic dependence of param[4] on the feature
    param[4] = param0 * unknown_response(f)

    # Calling the simulation function
    Simulation = ODEs_system_sim(
        model_1,  # Model function
        u0,       # Initial conditions
        0.0,      # Start time
        100.0,    # End time
        2.0,      # Time step
        param     # Model parameters
    )

    # Extracting simulation results
    sol_time = reduce(hcat, Simulation.t)
    sol_t = reduce(hcat, Simulation.u)

    # Adding uniform random noise to simulation data
    sol_t_noise = [sol_t[i, :] .+ rand(Uniform(-0.05, 0.05), size(sol_t)[2]) for i in 1:size(sol_t)[1]]
    sol_t_noise = permutedims(reduce(hcat, sol_t_noise))
    
    data = vcat(sol_time, sol_t_noise)  # Final dataset

    # Plot data with noise for different system variables
    display(scatter(data[1, :], data[2, :]))
    display(scatter!(data[1, :], data[3, :]))
    display(scatter!(data[1, :], data[4, :]))
    display(scatter!(data[1, :], data[5, :]))

    # Fit ODE system to noisy data
    fit = fit_ODEs_System(
        data,
        string(f),
        model_1, 
        param_guess,
        u0;
        lb=lb1,
        ub=ub1
    )

    display(plot!(fit[3]))  # Plot fitted model results

    # Storing fitted results for symbolic regression
    if f == feature_range[1]
        results_fit = fit[2]
    else
        results_fit = vcat(results_fit, reduce(hcat, fit[2][2, :]))
    end
end

# Scatter plot of feature value vs. estimated parameter p4
scatter(results_fit[2:end, 1], results_fit[2:end, 6], xlabel="Feature value", ylabel="p4")

# Setting options for symbolic regression
options = SymbolicRegression.Options(
    binary_operators=[+, /, *, -],
    unary_operators=[square],
    constraints=nothing,
    elementwise_loss=nothing,
    loss_function=nothing,
    tournament_selection_n=12,
    tournament_selection_p=0.86,
    topn=12,
    complexity_of_operators=nothing,
    complexity_of_constants=nothing,
    complexity_of_variables=nothing,
    parsimony=0.05,
    dimensional_constraint_penalty=nothing,
    alpha=0.100000,
    maxsize=10,
    maxdepth=nothing
)

# Generating feature matrix for symbolic regression
feature_matrix = [[string(f), f] for f in feature_range]
feature_matrix = permutedims(reduce(hcat, feature_matrix))

results_fit[:, 2] = results_fit[:, 1]
results_fit = permutedims(results_fit)

# Performing symbolic regression on the fitted results
gr_sy_reg = Kinbiont.downstream_symbolic_regression(results_fit, feature_matrix, 6; options=options)

# Plot results of symbolic regression
scatter(results_fit[2, 2:end], results_fit[6, 2:end], xlabel="Feature value", ylabel="Growth rate")
hline!(unique(gr_sy_reg[3][:, 1]), label=["Eq. 1" nothing], line=(3, :green, :dash))
plot!(unique(results_fit[2, 2:end]), unique(gr_sy_reg[3][:, 2]), label=["Eq. 2" nothing], line=(3, :red))
plot!(unique(results_fit[2, 2:end]), unique(gr_sy_reg[3][:, 3]), label=["Eq. 3" nothing], line=(3, :blue, :dashdot))

