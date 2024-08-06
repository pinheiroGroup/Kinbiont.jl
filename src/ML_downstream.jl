using DecisionTree
using MLJDecisionTreeInterface
using DecisionTree
using AbstractTrees

"""
    downstream_decision_tree_regression(
    KinBiont_results::Matrix{Any},
    feature_matrix::Matrix{Any},
    row_to_learn::Int;
    max_depth = -1,
    verbose = true,
    pruning_purity = 1.0,
    min_samples_leaf = 5,
    min_samples_split = 2,
    min_purity_increase = 0.0,
    n_subfeatures = 0,
    do_pruning = true,
    pruning_accuracy = 1.0,
    seed = 3,
    do_cross_validation = false,
    n_folds_cv = 3,
    )

This function performs regression using a decision tree algorithm on results from KinBiont. It includes options for tree pruning, cross-validation, and feature selection.

# Arguments:

- `KinBiont_results::Matrix{Any}`: The matrix containing results from fitting one or more files using KinBiont. 

- `feature_matrix::Matrix{Any}`: Matrix of features used for machine learning analysis. The number of rows should match the number of columns (minus one) in the `KinBiont_results`, with the first column containing well names to align with the well names in the second row of `KinBiont_results`.

- `row_to_learn::Int`: The index of the row in the `KinBiont_results` matrix that will be used as the target for the regression model.

# Key Arguments:

- `max_depth::Int = -1`: Maximum depth of the decision tree. If -1, there is no maximum depth.

- `verbose::Bool = true`: If true, the function will output additional details during the training process.

- `pruning_purity::Float64 = 1.0`: Purity threshold for pruning. If set to 1.0, no pruning will be performed.

- `min_samples_leaf::Int = 5`: Minimum number of samples required to be at a leaf node.

- `min_samples_split::Int = 2`: Minimum number of samples required to split an internal node.

- `min_purity_increase::Float64 = 0.0`: Minimum increase in purity required for a split.

- `n_subfeatures::Int = 0`: Number of features to select at random for splitting. If 0, all features are considered.

- `do_pruning::Bool = true`: If true, post-inference impurity pruning will be performed.

- `pruning_accuracy::Float64 = 1.0`: Purity threshold used for post-pruning. A value of 1.0 means no pruning.

- `seed::Int = 3`: Random seed for reproducibility.

- `do_cross_validation::Bool = false`: If true, performs n-fold cross-validation.

- `n_folds_cv::Int = 3`: Number of folds for cross-validation. Ignored if `do_cross_validation` is false.

# Outputs:

1. `tree_model::Any`: The trained decision tree model.

2. `importance_score::Vector{Float64}`: Importance scores of each feature used in the model.

3. `importance_rank::Vector{Int}`: Ranking of features based on their importance scores.

4. `cross_validation_score::Union{Float64, Nothing}`: Cross-validation score if `do_cross_validation` is true; otherwise, `nothing`.

5. `samples_per_leaf::Matrix{Int}`: Matrix where each row represents a leaf node and the columns represent the samples associated with that leaf.


"""
function downstream_decision_tree_regression(KinBiont_results::Matrix{Any}, # output of KinBiont results
  feature_matrix::Matrix{Any},
  row_to_learn::Int;
  max_depth = -1,
  verbose = true,
  pruning_purity = 1.0,
  min_samples_leaf = 5,
  min_samples_split = 2,  
  min_purity_increase = 0.0, 
  n_subfeatures = 0,
  do_pruning = true,
  pruning_accuracy = 1.0,
  seed = 3,
  do_cross_validation = false,
  n_folds_cv = 3,
)
 max_depth = convert(Int, max_depth)


  feature_names = string.(feature_matrix[1,2:end])

  names_of_the_wells_res = KinBiont_results[2, 1:end]
  names_of_the_wells_annotation = feature_matrix[1:end, 1]
  wells_to_use = intersect(names_of_the_wells_res, names_of_the_wells_annotation)


  index_res = Int
  index_annotation = Int

  index_res =  findfirst.(isequal.(wells_to_use), (names_of_the_wells_res,))
  index_annotation =findfirst.(isequal.(wells_to_use), (names_of_the_wells_annotation,))

  index_res = index_res[index_res.!=nothing]
  index_annotation =index_annotation[index_annotation.!=nothing]

  # order results and annotation by well in the same order

  #index_res[1, :] = index_res[1, :] 

  output = convert.(Float64, KinBiont_results[row_to_learn, index_res])

  predictors = convert.(Float64, feature_matrix[index_annotation, 2:end])





  model = build_tree(output, predictors, 
    0,
    max_depth,
    min_samples_leaf,
    min_samples_split,
    min_purity_increase;
    rng = seed)

  if do_pruning ==  true
    model = prune_tree(model, pruning_accuracy)
  end 

  r2 = Any
  if do_cross_validation ==  true
    r2 =  nfoldCV_tree(output, predictors,
    n_folds_cv,
    pruning_purity,
    max_depth,
    min_samples_leaf,
    min_samples_split,
    min_purity_increase;
    verbose = verbose,
    rng = seed)

  end

  imp_1 = impurity_importance(model)
  imp_2 = split_importance(model)
  #imp_3 = permutation_importance(model)

  wt = DecisionTree.wrap(model, (featurenames = feature_names,))

  leaves = collect(AbstractTrees.Leaves(wt))
  wt = DecisionTree.wrap(model, (featurenames = feature_names,))
  leaves_values = ["values","cluster"]
  values = Any

  for i in eachindex(leaves)
      values = leaves[i].leaf.values
      index_leaf = zeros(length(values))
      index_leaf .= i
      temp = hcat(index_leaf,values)
      leaves_values = hcat(leaves_values,transpose(temp))

  end

  if verbose == true
    tree_out = DecisionTree.print_tree(model, max_depth)
  end
    return model, imp_1, imp_2, r2 , leaves_values
end






"""
    downstream_symbolic_regression(
    KinBiont_results::Matrix{Any},
    feature_matrix::Matrix{Any},
    row_to_learn::Int;
    options = SymbolicRegression.Options(),
    )

This function performs symbolic regression on the results obtained from fitting models using KinBiont. It uses a feature matrix to train a symbolic regression model to predict a specific row of the KinBiont results.

# Arguments:

- `KinBiont_results::Matrix{Any}`: The matrix containing results from fitting one or more files using KinBiont. 

- `feature_matrix::Matrix{Any}`: Matrix of features used for machine learning analysis. The number of rows in this matrix should match the number of columns (minus one) in the `KinBiont_results`, with the first column containing well names to match the features with the well names in the second row of `KinBiont_results`.

- `row_to_learn::Int`: The index of the row in the `KinBiont_results` matrix that will be the target for machine learning inference.

# Key Arguments:

- `options::SymbolicRegression.Options()`: Options for the symbolic regression process. This argument uses the `SymbolicRegression.Options` class, allowing customization of the symbolic regression parameters. See [SymbolicRegression.jl API documentation](https://astroautomata.com/SymbolicRegression.jl/stable/api/#SymbolicRegression.CoreModule.OptionsStructModule.Options) for details.

# Outputs:

- `trees::Vector{SymbolicRegression.Tree}`: A vector of trees representing the hall of fame results from the symbolic regression process.

- `res_output::Matrix{Any}`: A matrix containing the hall of fame results, where each row includes:
1. Complexity score of the equation.
2. Mean Squared Error (MSE) of the equation.
3. The symbolic equation itself.

- `predictions::Matrix{Float64}`: For each equation in the hall of fame, this matrix contains the predicted values for each sample. Columns represent the different equations, and rows correspond to the samples.

- `index_annotation::Vector{Int}`: An index vector used to order the rows of the `feature_matrix` to match the columns of the `KinBiont_results`.

"""
function downstream_symbolic_regression(kinBiont_results,
  feature_matrix,
  row_to_learn;
  options = SymbolicRegression.Options(),
)





  names_of_the_wells_res = KinBiont_results[2, 1:end]
  names_of_the_wells_annotation = feature_matrix[1:end, 1]
  wells_to_use = intersect(names_of_the_wells_res, names_of_the_wells_annotation)


  index_res = Int
  index_annotation = Int





  index_res =  findfirst.(isequal.(wells_to_use), (names_of_the_wells_res,))
  index_annotation =findfirst.(isequal.(wells_to_use), (names_of_the_wells_annotation,))

  index_res = index_res[index_res.!=nothing]
  index_annotation =index_annotation[index_annotation.!=nothing]

  # order results and annotation by well in the same order

  #index_res[1, :] = index_res[1, :] 

  output = convert.(Float64, KinBiont_results[row_to_learn, index_res])


  predictors =Matrix(transpose(  convert.(Float64, feature_matrix[index_annotation,2])))

  #hall_of_fame = SymbolicRegression.equation_search(predictors,output,options =options)
  hall_of_fame = SymbolicRegression.equation_search(predictors,output)

  dominating = calculate_pareto_frontier(hall_of_fame)
  trees = [member.tree for member in dominating]



  res_output = ["Complexity","MSE","Equation"]

  predictions = Any

  for t in trees


    output2, did_succeed = eval_tree_array(t, predictors, options)
    if t == trees[1]
      predictions =  output2
    else

      predictions =  hcat(predictions,output2)

    end


  end



  for member in dominating
      complexity = compute_complexity(member, options)
      loss = member.loss
      string = string_tree(member.tree, options)

      res_output = hcat(res_output,[complexity,loss,string])
  end




  return trees,res_output,predictions,index_annotation

end










export downstream_decision_tree_regression
export downstream_symbolic_regression