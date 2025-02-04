using Kinbiont
using DifferentialEquations
using CSV
using Plots
using StatsBase
using AbstractTrees
using MLJDecisionTreeInterface
using TreeRecipe
using DecisionTree
using Distributions
using  Random
# Generate a dataset with an unknown dependence on a feature 


# defining the used ODE model 
results_fit =  Any

ODE_models = "baranyi_richards"

ub_1 = [0.2, 5.1, 500.0, 5.0]
lb_1 = [0.0001, 0.2, 0.00, 0.2]
p1_guess = lb_1 .+ (ub_1 .- lb_1) ./ 2


# defining the feature matrix of each simulation. suppose you have 3different antibiotics each with one  concentration


function transform_abx_vector(input_vector::Vector, mu::Float64)
    
    # Define concentration mapping rules
    concentration_map = Dict(
        (1, 0, 0) => 1.0 ,    # abx_1 -> μ
        (0, 1, 0) => 0.5 ,    # abx_2 -> 0.5μ
        (0, 0, 1) => 0.3 ,    # abx_3 -> 0.6μ
        (1, 1, 0) => 0.0 ,    # abx_1 + abx_2 -> 0μ
        (1, 0, 1) => 0.3 ,    # abx_1 + abx_3 -> 0.6μ
        (0, 1, 1) => 0.0 ,    # abx_2 + abx_3 -> 0μ
        (1, 1, 1) => 0.0,    # abx_1 + abx_2 + abx_3 -> 0.0μ
        (0, 0, 0) =>  1.0     # No antibiotics -> 0μ
    )
    
    mu_correct = concentration_map[Tuple(input_vector[2:end])] * mu   # Default to 0μ if not found
    return mu_correct
end

# generating random the matrix of the features
# Define the dimensions of the matrix
cols = 3
n_experiment = 100

# Generate a random matrix with 0s and 1s
random_matrix = rand(0:1, n_experiment, cols)
labels = string.(1:1:n_experiment)
random_matrix = hcat(labels,random_matrix)
random_matrix
# defining the parameters values for the simulation 
p_sim = [0.05, 1.0, 50.0, 1.0]
psim_1_0 =  p_sim[1]
p1_array = [transform_abx_vector(random_matrix[f,:], psim_1_0 ) for  f in 1:size(random_matrix)[1]]

t_min = 0.0
t_max = 800.0
n_start = [0.1]
delta_t = 10.0
noise_value = 0.03

plot(0, 0)
for f in 1:size(random_matrix)[1]

    # changing the parameters with unknown perturbation 

    p_sim[1] =transform_abx_vector(random_matrix[f,:], psim_1_0 )


    # Calling the simulation function
    sim = Kinbiont.ODE_sim("baranyi_richards", n_start, t_min, t_max, delta_t, p_sim)

    # Plotting scatterplot of data without noise

    #adding uniform random noise
    noise_unifom = rand(Uniform(-noise_value, noise_value), length(sim.t))


    data_t = reduce(hcat, sim.t)
    data_o = reduce(hcat, sim.u)
    data_OD = vcat(data_t, data_o)
    data_OD[2, :] = data_OD[2, :] .+ noise_unifom
    # ploting scatterplot of data with noise

    display(Plots.scatter!(data_OD[1, :], data_OD[2, :], xlabel="Time", ylabel="Arb. Units", label= nothing, color=:red, markersize=2, size=(300, 300)))

    results_ODE_fit = fitting_one_well_ODE_constrained(
        data_OD,
        string(random_matrix[f,1]),
        "test_ODE",
        "baranyi_richards",
        p1_guess;
        lb=lb_1,
        ub=ub_1
    )
    display(Plots.plot!(results_ODE_fit[4], results_ODE_fit[3], xlabel="Time", ylabel="Arb. Units", label= nothing, color=:red, markersize=2, size=(300, 300)))


    if f == 1
        results_fit = results_ODE_fit[2]
    else
        results_fit = hcat(results_fit, results_ODE_fit[2])
    end



end


# generating feature matrix
# the first column is the label as a string of the feature value we used for the fitting labeling


n_folds = 10
feature_names = ["abx_1","abx_2","abx_3"]

depth = -1 


# Set random seed for reproducibility
seed = Random.seed!(1234)



feature_matrix = vcat(["label" "abx_1" "abx_2" "abx_3"],random_matrix)

dt_gr = Kinbiont.downstream_decision_tree_regression(results_fit,
        feature_matrix,
        4;# row to learn
        do_pruning=false,
        pruning_accuracy=1.00,
        verbose=true,
        do_cross_validation=true,
        max_depth=depth,
        n_folds_cv=n_folds,
        seed=seed
    )


# Wrap the decision tree model for visualization
wt = DecisionTree.wrap(dt_gr[1], (featurenames = feature_names,))

# Plot the decision tree
p2 = Plots.plot(wt, 0.9, 0.2; size = (900, 400), connect_labels = ["yes", "no"])