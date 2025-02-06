using DifferentialEquations
using CSV
using Plots
using StatsBase
using AbstractTrees
using MLJDecisionTreeInterface
using TreeRecipe
using DecisionTree
using Distributions
using Random

# In this example, we will use Kinbiont to generate data about a simple community of tree species (N_1, N_2, and N_3). 
# We suppose that the user can only measure the total biomass of the community (N_1 + N_2 + N_3)
# but they can vary the initial composition of the community. By fitting the biomass with a single empirical model, 
# we apply a decision tree algorithm using the initial composition of the community as features and 
# the parameters of the model as the quantity to predict. We show how we can decompose the effect 
# of the initial population size on the behavior of the full population.

# We initialize the array for the results
results_fit = Any

# Tree species ODE model: they compete for the same resource u[4], and u[3] is able to predate on u[1]
function model_1(du, u, param, t)
    # Define the ODEs
    du[1] = param[1] * u[1] * u[4] - param[4] * u[3] * u[1] 
    du[2] = param[2] * u[2] * u[4]
    du[3] = param[3] * u[2] * u[4] + param[4] * u[3] * u[1]
    du[4] = -param[1] * (u[1] + u[2] + u[3]) * u[4]
end

# Parameters [1/yield_1, 1/yield_2, 1/yield_3, predation (3->1)]
param = [0.05, 0.02, 0.03, 0.03]

# We define the noise to add to the simulated data
noise_value = 0.01

# We declare the ODE model, its upper/lower bounds, and initial guess
# Parameters to fit
ODE_models = "HPM"
ub_1 = [0.5, 5.1, 16.0]
lb_1 = [0.0001, 0.000001, 0.00]
p1_guess = lb_1 .+ (ub_1 .- lb_1) ./ 2

# We generate a random matrix of features (i.e., initial conditions)
# Define the dimensions of the matrix
cols = 3
n_experiment = 150

# Generate a random matrix with 0s and 0.1s
random_matrix = rand([0, 0.1], n_experiment, cols)
labels = string.(1:1:n_experiment)
random_matrix = hcat(labels, random_matrix)
random_matrix

# Defining the parameter values for the simulation 
t_min = 0.0
t_max = 200.0
delta_t = 8.0
noise_value = 0.05
blank_value = 0.08
plot(0, 0)

for f in 1:size(random_matrix)[1]

    # Defining the initial condition with different community compositions
    u0 = random_matrix[f, 2:4]
    u0 = push!(u0, 5.0)
    u0 = convert.(Float64, u0)

    # Calling the simulation function
    Simulation = ODEs_system_sim(
        model_1, # ODE system 
        u0, # Initial conditions
        t_min, # Start time of the simulation
        t_max, # Final time of the simulation
        delta_t, # Delta t for Poisson approximation
        param; # Parameters of the ODE model
    )

    # Plotting scatterplot of data without noise

    # Adding uniform random noise
    noise_uniform = rand(Uniform(-noise_value, noise_value), length(Simulation.t))

    data_t = reduce(hcat, Simulation.t)
    data_o = reduce(hcat, Simulation.u)
    data_o = data_o[1:3, :]

    data_o = sum(data_o, dims=1)

    data_OD = vcat(data_t, data_o)
    data_OD[2, :] = data_OD[2, :] .+ noise_uniform .+ blank_value
    # Plotting scatterplot of data with noise

    display(Plots.scatter!(data_OD[1, :], data_OD[2, :], xlabel="Time", ylabel="Arb. Units", label=nothing, color=:red, markersize=2, size=(300, 300)))

    # Fitting with Kinbiont
    results_ODE_fit = fitting_one_well_ODE_constrained(
        data_OD,
        string(random_matrix[f, 1]),
        "test_ODE",
        ODE_models,
        p1_guess;
        remove_negative_value=true,
        pt_avg=2,
        lb=lb_1,
        ub=ub_1
    )

    display(Plots.plot!(results_ODE_fit[4], results_ODE_fit[3], xlabel="Time", ylabel="Arb. Units", label=nothing, color=:red, markersize=2, size=(300, 300)))

    # Storing results
    if f == 1
        results_fit = results_ODE_fit[2]
    else
        results_fit = hcat(results_fit, results_ODE_fit[2])
    end
end

# Parameters and feature names for the decision tree

n_folds = 10
feature_names = ["CI_1", "CI_2", "CI_3"]
depth = -1 
seed = Random.seed!(1234)
feature_matrix = vcat(["label" "CI_1" "CI_2" "CI_3"], random_matrix)

# Calling the decision tree regression for the 6th row of results_fit. It represents the saturation value N_max
dt_gr = Kinbiont.downstream_decision_tree_regression(results_fit,
    feature_matrix,
    6; # Row to learn
    do_pruning=false,
    verbose=true,
    do_cross_validation=true,
    max_depth=depth,
    n_folds_cv=n_folds,
    seed=seed
)

# Wrap the decision tree model for visualization
wt = DecisionTree.wrap(dt_gr[1], (featurenames=feature_names,))

# Plot the decision tree
p2 = Plots.plot(wt, 0.9, 0.2; size=(1500, 700), connect_labels=["yes", "no"])

# Calling the decision tree regression for the 7th row of results_fit. It represents the maximum per capita growth rate
dt_gr = Kinbiont.downstream_decision_tree_regression(results_fit,
    feature_matrix,
    7; # Row to learn
    do_pruning=false,
    verbose=true,
    do_cross_validation=true,
    max_depth=depth,
    n_folds_cv=n_folds,
    seed=seed
)

# Wrap the decision tree model for visualization
wt = DecisionTree.wrap(dt_gr[1], (featurenames=feature_names,))

# Plot the decision tree
p2 = Plots.plot(wt, 0.9, 0.2; size=(1500, 700), connect_labels=["yes", "no"])
