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
 #tre species ODE model
function model_1(du, u, param, t)
    # Define the ODEs
    du[1] = param[1] * u[1] * u[4] - param[4]  * u[3] * u[1] 
    du[2] = param[2] * u[2] *  u[4]
    du[3] = param[3] * u[2]  * u[4] + param[4]  * u[3] * u[1]
    du[4]=   - param[1] * ( u[1] + u[2] + u[3] ) * u[4]
end
u0 = [0.1, 0.1, 0.1,5.0]  # Initial conditions for the variables

# Parameters
param = [0.05, 0.02, 0.03,  0.03]

param_guess = lb1 .+ (ub1 .- lb1) ./ 2

param0 = param[4]
noise_value = 0.01
# parameters to fit 
ODE_models = "HPM"

ub_1 = [0.5, 5.1, 16.0]
lb_1 = [0.0001, 0.000001, 0.00]
p1_guess = lb_1 .+ (ub_1 .- lb_1) ./ 2


# generating random the matrix of the features (i.e. initial conditions)
# Define the dimensions of the matrix
cols = 3
n_experiment = 150

# Generate a random matrix with 0s and 1s
random_matrix = rand([0,0.1], n_experiment, cols)
labels = string.(1:1:n_experiment)
random_matrix = hcat(labels,random_matrix)
random_matrix
# defining the parameters values for the simulation 

t_min = 0.0
t_max = 200.0
delta_t =8.0
noise_value = 0.05
blank_value = 0.08
plot(0, 0)


for f in 1:size(random_matrix)[1]

    # changing the parameters with unknown perturbation 
     println(f)
    u0 = random_matrix[f,2:4]
    u0 = push!(u0,5.0)
    u0 = convert.(Float64, u0)

    # Calling the simulation function
  # Calling the simulation function
   Simulation =  ODEs_system_sim(
    model_1, #string of the model
    u0, # starting condition
    t_min, # start time of the sim
    t_max, # final time of the sim
    delta_t, # delta t for poisson approx
    param; # parameters of the ODE model
    )
    # Plotting scatterplot of data without noise

    #adding uniform random noise
    noise_unifom = rand(Uniform(-noise_value, noise_value), length(Simulation.t))


    data_t = reduce(hcat, Simulation.t)
    data_o = reduce(hcat, Simulation.u)
    data_o = data_o[1:3,:]

    data_o = sum(data_o, dims=1)


    data_OD = vcat(data_t, data_o)
    data_OD[2, :] = data_OD[2, :] .+ noise_unifom .+ blank_value
    # ploting scatterplot of data with noise

    display(Plots.scatter!(data_OD[1, :], data_OD[2, :], xlabel="Time", ylabel="Arb. Units", label= nothing, color=:red, markersize=2, size=(300, 300)))
   # if the initial is  0.0,0.0,0.0,0.5, do not fit the model and set all numeric value of results_fit to 0






    results_ODE_fit = fitting_one_well_ODE_constrained(
        data_OD,
        string(random_matrix[f,1]),
        "test_ODE",
        ODE_models,
        p1_guess;
        remove_negative_value= true,
        pt_avg =2,
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
feature_names = ["CI_1","CI_2","CI_3"]

depth = -1 


# Set random seed for reproducibility
seed = Random.seed!(1234)



feature_matrix = vcat(["label" "CI_1" "CI_2" "CI_3"],random_matrix)

dt_gr = Kinbiont.downstream_decision_tree_regression(results_fit,
        feature_matrix,
        6;# row to learn
        do_pruning=false,
      #  pruning_accuracy=0.025,
        verbose=true,
        do_cross_validation=true,
        max_depth=depth,
        n_folds_cv=n_folds,
        seed=seed
    )


# Wrap the decision tree model for visualization
wt = DecisionTree.wrap(dt_gr[1], (featurenames = feature_names,))

# Plot the decision tree
p2 = Plots.plot(wt, 0.9, 0.2; size = (1500, 700), connect_labels = ["yes", "no"])







dt_gr = Kinbiont.downstream_decision_tree_regression(results_fit,
        feature_matrix,
        7;# row to learn
        do_pruning=false,
      #  pruning_accuracy=0.025,
        verbose=true,
        do_cross_validation=true,
        max_depth=depth,
        n_folds_cv=n_folds,
        seed=seed
    )


# Wrap the decision tree model for visualization
wt = DecisionTree.wrap(dt_gr[1], (featurenames = feature_names,))

# Plot the decision tree
p2 = Plots.plot(wt, 0.9, 0.2; size = (1500, 700), connect_labels = ["yes", "no"])