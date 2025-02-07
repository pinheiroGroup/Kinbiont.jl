# [Examples and Tutorial](@id examples)

This section provides some copy-and-paste examples of Kinbiont.jl

```@contents
Pages = ["index.md"]
Depth = 3
```
To run all these example you need to call the following packages:


```julia
using Kinbiont
using CSV
using Plots
using Distributions
using Optimization
using OptimizationNLopt
using Tables
using SymbolicRegression
using DataFrames
using Statistics
using DelimitedFiles
using Random
using DecisionTree
using AbstractTrees
using MLJDecisionTreeInterface
using TreeRecipe
```
## Symbolic regression detection of laws

## Decision tree regression


## Symbolic regression on real data

This example demonstrates how to use the `Kinbiont` and `SymbolicRegression` packages to analyze kinetics data.

Set up paths to your data, annotation, calibration curve, and result directories (see examples):

```julia
path_to_data = "your_path/data_examples/plate_data.csv"
path_to_annotation = "your_path/data_examples/annotation.csv"
path_to_calib = "your_path/data_examples/cal_curve_example.csv"
path_to_results = "your_path//seg_res/"
```

We fit the data with  segmentation and 1 change point. First, we declare the models

```julia
model1 = "HPM_exp"
lb_param1 = [0.00001, 0.000001]
ub_param1 = [0.5, 1.5]
param_guess1 = [0.01, 0.01]

model2 = "logistic"
lb_param2 = [0.00001, 0.000001]
ub_param2 = [0.5, 2.5]
param_guess2 = [0.01, 1.01]

list_of_models = [model1, model2]
list_guess = [param_guess1, param_guess2]
list_lb = [lb_param1, lb_param2]
list_ub = [ub_param1, ub_param2]

n_change_points =1
```

We perform  the fitting:

```julia
fit_file = Kinbiont.segmentation_ODE_file(
    "seg_exp_1", # Label of the experiment
    path_to_data, # Path to the folder to analyze
    list_of_models, # ODE models to use
    list_guess, # Parameter guesses
    n_change_points;
    path_to_annotation=path_to_annotation, # Path to the annotation of the wells
    detect_number_cpd=false,
    fixed_cpd=false,
    path_to_calib=path_to_calib,
    multiple_scattering_correction=false, # If true, use calibration curve for data correction
    type_of_curve="deriv",
    pt_smooth_derivative=10,
    type_of_smoothing="lowess",
    verbose=true,
    write_res=false,
    path_to_results=path_to_results,
    win_size=12, # Number of points to generate initial condition
    smoothing=true,
    lb_param_array=list_lb, # Lower bounds for parameters
    ub_param_array=list_ub, # Upper bounds for parameters
    maxiters=200000
)
```
We load the annotation and  select results of a specific strain

```julia
annotation_test = CSV.File(path_to_annotation, header=false)
names_of_annotation = propertynames(annotation_test)
feature_matrix = hcat(annotation_test[:Column1], annotation_test[:Column3])

# Selecting strain S5
index_strain = findall(annotation_test[:Column2] .== "S5")
index_cc = findall(annotation_test[:Column3] .> 0.01)
to_keep = intersect(index_strain, index_cc)
feature_matrix = feature_matrix[to_keep, :]
wells_to_use = feature_matrix[:, 1]
index_res = Any

```

We give same order to the well to analyze and the features (can be skipped)
```julia
for i in wells_to_use
    if i == wells_to_use[1]
        index_res = findfirst(res_first_seg[:, 1:end] .== i)
    else
        index_res = hcat(index_res, findfirst(res_first_seg[:, 1:end] .== i))
    end
end

iii = [index_res[1, i][2] for i in eachindex(index_res[1, :])]
res_first_seg_ML = res_first_seg[:, iii]
res_first_seg_ML = hcat(res_first_seg[:, 1], res_first_seg_ML)
```

We add x = 0.0, y = 0.0 to data to take in consideration  not growing wells:

```julia

feature_matrix =vcat(feature_matrix,["zero" 0.0])
res_first_seg_ML=hcat(res_first_seg_ML , reduce(vcat,["zero" ,"zero", "zero", 0.0 ,  0.0 ,0.0 ,0.0 ,0.0 ,0.0 ]))


```

We Declare the options of symbolic regression:

```julia

options = SymbolicRegression.Options(
    binary_operators=[+, /, *, -],
    unary_operators=[],
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
```

We run the symbolic regression using dependent variable that is the 7th row of the Kinbiont results (i.e., the growth rate)

```julia

gr_sy_reg = Kinbiont.downstream_symbolic_regression(res_first_seg_ML, feature_matrix, 7; options=options)

scatter(feature_matrix[:, 2], res_first_seg_ML[7, 2:end], xlabel="Amino Acid concentration μM", ylabel="Growth rate [1/Min]", label=["Data" nothing])
hline!(unique(gr_sy_reg[3][:, 1]), label=["Eq. 1" nothing], line=(3, :green, :dash))
plot!(unique(convert.(Float64, feature_matrix[gr_sy_reg[4], 2])), unique(gr_sy_reg[3][:, 2]), label=["Eq. 2" nothing], line=(3, :red))
plot!(unique(convert.(Float64, feature_matrix[gr_sy_reg[4], 2])), unique(gr_sy_reg[3][:, 3]), label=["Eq. 3" nothing], line=(3, :blue, :dashdot))
plot!(unique(convert.(Float64, feature_matrix[gr_sy_reg[4], 2])), unique(gr_sy_reg[3][:, 4]), label=["Eq. 4" nothing], line=(2, :black))
plot!(unique(convert.(Float64, feature_matrix[gr_sy_reg[4], 2])), unique(gr_sy_reg[3][:, 5]), label=["Eq. 5" nothing], line=(2, :black))

```

### Decision Tree Regression Analysis on real data



In this example we analyze the already fitted data from:


We read the data from CSV files:

```julia
Kinbiont_res_test = readdlm("your_path/data_examples/Results_for_ML.csv", ',')
annotation_test = readdlm("your_path/data_examples/annotation_for_ML.csv", ',')
```


We define some variables for analysis:

```julia
ordered_strain = annotation_test[:, end]
n_folds = 10
feature_names = unique(annotation_test[1, 2:end])[2:(end-1)]

depth = -1 


# Set random seed for reproducibility
seed = Random.seed!(1234)
```
We perform decision tree regression only on "N. soli" strain and we  analyze the 9th row (i.e. growht rate) of the results:

```julia


index_strain = findall("N. soli".== ordered_strain)
feature_matrix = annotation_test[index_strain, 2:(end-1)]
Kinbiont_results = Kinbiont_res_test[:, index_strain]

dt_gr = Kinbiont.downstream_decision_tree_regression(Kinbiont_results,
        feature_matrix,
        9;# row to learn
        do_pruning=false,
        pruning_accuracy=1.00,
        verbose=true,
        do_cross_validation=true,
        max_depth=depth,
        n_folds_cv=n_folds,
        seed=seed
    )
```

We plot the tree
```julia
# Wrap the decision tree model for visualization
wt = DecisionTree.wrap(dt_gr[1], (featurenames = feature_names,))

# Plot the decision tree
p2 = Plots.plot(wt, 0.9, 0.2; size = (900, 400), connect_labels = ["yes", "no"])
```
