using Kinbiont
using CSV
using DataFrames
using Plots
using Statistics
# data from "Construction and Modeling of a Coculture Microplate for Real-Time Measurement of Microbial Interactions"
# strain annotation 
A0 = [1,11,21,31,41,51,13,33,53,7,17,37,47,57,38,48,58]
LP = [12,32,52,5,15,25,35,45,55,27,8,9,19,29,10,20,30]
LB = [4,14,24,34,44,54,16,36,56,18,28,39,49,59,40,50,60]

# load data
data = CSV.read("../data_examples/DATA_Gut_Microbiome_CoC.csv", DataFrame)
# convet time column from in elapsed time from the start in hours, the delta time in the data is 15 minutes

delta_time = 15/60
time_seq = 0:delta_time:(size(data,1)-1)*delta_time
time_seq = [ time_seq[i] for i in 1: length(time_seq)]
data[!,:Time] = time_seq 

index_tot = 2:1:size(data,2)
A0_index = A0.+1
LP_index = LP.+1
LB_index = LB.+1
blanS1_index = setdiff(index_tot, [A0_index; LP_index; LB_index])

### evaluating blank values by doing the mean of the blank columns
data_blanS1_values = Matrix(data[:,blanS1_index])
#blanS1_value = mean.(reduce(vcat,eachrow(data_blanS1_values)))
blanS1_value =  mean(data_blanS1_values,dims = (2))

# Intializing the model for fitting with Kinbiont
model ="aHPM"
results_matrix = Kinbiont.initialize_df_results(model)

p_guess = [0.2, 0.001, 1.00, 1.0]
ub_ahpm = [3.0, 1.0,3.00, 5.0]
lb_ahpm = [0.001, 0.0,0.2, 0.0]
# Fitting all AO data
# Not interacting A0
# AO, -
# [1,21,41];
A0 = [1,21,41].+1

data_A0 = data[:,A0] #.- blanS1_value

# mean of the replicates
data_A0_mean = mean(Matrix(data_A0),dims=2)

# Formatting data for Kinbiont
data_to_fit = permutedims( [time_seq data_A0_mean])

# Performing ODE fitting
results_ODE_fit = Kinbiont.fitting_one_well_ODE_constrained(
    data_to_fit, 
    "AO_NA",
    "CoC_DT",
    model,
    p_guess;
    lb = lb_ahpm,
    ub = ub_ahpm,
)
Plots.scatter(data_to_fit[1,:],data_to_fit[2,:], xlabel="Time [h]", ylabel="OD [Arb. Units]", label=["Data " nothing],color=:black,markersize =2 ,size = (300,300))
Plots.plot!(results_ODE_fit[4],results_ODE_fit[3], xlabel="Time [h]", ylabel="OD [Arb. Units]",label=["fit " nothing],color=:red,markersize =4 ,size = (300,300))
results_matrix = [results_matrix results_ODE_fit[2]]


#AO, LP
# [11,31,51,7]; AO


A0 =  [11,31,51,7].+1

data_A0 = data[:,A0] #.- blanS1_value

# mean of the replicates
data_A0_mean = mean(Matrix(data_A0),dims=2)

# Formatting data for Kinbiont
data_to_fit = permutedims( [time_seq data_A0_mean])

# Performing ODE fitting
results_ODE_fit = Kinbiont.fitting_one_well_ODE_constrained(
    data_to_fit, 
    "AO_LP",
    "CoC_DT",
    model,
    p_guess;
    lb = lb_ahpm,
    ub = ub_ahpm,
)

Plots.scatter(data_to_fit[1,:],data_to_fit[2,:], xlabel="Time [h]", ylabel="OD [Arb. Units]", label=["Data " nothing],color=:black,markersize =2 ,size = (300,300))
Plots.plot!(results_ODE_fit[4],results_ODE_fit[3], xlabel="Time [h]", ylabel="OD [Arb. Units]",label=["fit " nothing],color=:red,markersize =4 ,size = (300,300))
results_matrix = [results_matrix results_ODE_fit[2]]





# AO, LB
# [13,33,53,17]; 


A0 = [13,33,53,17].+1

data_A0 = data[:,A0] #.- blanS1_value

# mean of the replicates
data_A0_mean = mean(Matrix(data_A0),dims=2)

# Formatting data for Kinbiont
data_to_fit = permutedims( [time_seq data_A0_mean])

# Performing ODE fitting
results_ODE_fit = Kinbiont.fitting_one_well_ODE_constrained(
    data_to_fit, 
    "AO_LB",
    "CoC_DT",
    model,
    p_guess;
    lb = lb_ahpm,
    ub = ub_ahpm,
)
Plots.scatter(data_to_fit[1,:],data_to_fit[2,:], xlabel="Time [h]", ylabel="OD [Arb. Units]", label=["Data " nothing],color=:black,markersize =2 ,size = (300,300))
Plots.plot!(results_ODE_fit[4],results_ODE_fit[3], xlabel="Time [h]", ylabel="OD [Arb. Units]",label=["fit " nothing],color=:red,markersize =4 ,size = (300,300))
results_matrix = [results_matrix results_ODE_fit[2]]


# AO, AO
# [37,47,57];



A0 =[37,47,57].+1

data_A0 = data[:,A0] #.- blanS1_value

# mean of the replicates
data_A0_mean = mean(Matrix(data_A0),dims=2)

# Formatting data for Kinbiont
data_to_fit = permutedims( [time_seq data_A0_mean])

# Performing ODE fitting
results_ODE_fit = Kinbiont.fitting_one_well_ODE_constrained(
    data_to_fit, 
    "AO_AO",
    "CoC_DT",
    model,
    p_guess;
    lb = lb_ahpm,
    ub = ub_ahpm,
)
Plots.scatter(data_to_fit[1,:],data_to_fit[2,:], xlabel="Time [h]", ylabel="OD [Arb. Units]", label=["Data " nothing],color=:black,markersize =2 ,size = (300,300))
Plots.plot!(results_ODE_fit[4],results_ODE_fit[3], xlabel="Time [h]", ylabel="OD [Arb. Units]",label=["fit " nothing],color=:red,markersize =4 ,size = (300,300))
results_matrix = [results_matrix results_ODE_fit[2]]




# Creating feature Matrix


label_row = ["AO_NA","AO_LP","AO_LB","AO_AO"]
starting_0 = zeros(Int,4,3)
feature_matrix = hcat(label_row,starting_0)
feature_matrix[2,2] = 1 
feature_matrix[3,3] = 1 
feature_matrix[4,4] = 1 
feature_names = ["LP","LB","AO"]


# DT regression of Nmax against the feature matrix
using Random

# Parameters of the decision tree
depth = -1  # No depth limit

# Set random seed for reproducibility
seed = Random.seed!(1234)


# Decision tree regression ON Nmax
dt = Kinbiont.downstream_decision_tree_regression(results_matrix,
        feature_matrix,
        6;  # Row to learn
        do_pruning=false,
        verbose=true,
        do_cross_validation=false,
        max_depth=depth, 
        seed=seed,
        min_samples_leaf =1,
        min_purity_increase= 0.001, 
        min_samples_split=2
    )

# Visualizing the decision tree
using DecisionTree
using AbstractTrees
using MLJDecisionTreeInterface
using TreeRecipe
wt = DecisionTree.wrap(dt[1])
p2 = Plots.plot(wt, 0.9, 0.2; size=(1400, 700), connect_labels=["yes", "no"])


# Decision tree regression exponential max growht rate
dt = Kinbiont.downstream_decision_tree_regression(results_matrix,
        feature_matrix,
        8;  # Row to learn
        do_pruning=false,
        verbose=true,
        do_cross_validation=false,
        max_depth=depth, 
        seed=seed,
        min_samples_leaf =1,
        min_purity_increase= 0.0001, 
        min_samples_split=2
    )

# Visualizing the decision tree
using DecisionTree
using AbstractTrees
using MLJDecisionTreeInterface
using TreeRecipe
wt = DecisionTree.wrap(dt[1])
p2 = Plots.plot(wt, 0.9, 0.2; size=(1400, 700), connect_labels=["yes", "no"])
