using SymbolicRegression
using MLJ
using Zygote
using DecisionTree
using MLJDecisionTreeInterface
using DecisionTree
using AbstractTrees

"""
    downstream_decision_tree_regression(
    kimchi_results::Matrix{Any},
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


Function that 

# Arguments:

- `kimchi_results::Matrix{Any}`: The matrix of results of fitting one file with Kimchi. Compatible functions `fit_file_ODE`,`fit_file_custom_ODE`, `ODE_model_selection_file`, `segmentation_ODE_file`, `fit_NL_model_file`, `fit_NL_model_selection_file`, and `fit_NL_segmentation_file`.
- `feature_matrix::Matrix{Any}`: Matrix of the features for the ML analysis. Important, the number of rows of this file should be te number of columns (minus one) of the kimchi_results, and the first column should have a name of the well in order to mach the feature with the names of the well in second row of kimchi_results. 
- `row_to_learn::Int`: which row of the matrix `kimchi_results` will be the target of the ML inference.

# Key Arguments:

- `max_depth = -1`, Int, maximum depth of the decision tree ( -1, no maximum)
- `verbose = true`
- `pruning_purity = 1.0`
- `min_samples_leaf = 5`: Int, the minimum number of samples each leaf needs to have
- `min_samples_split = 2`: Int, the minimum number of samples in needed for a split (
- `min_purity_increase = 0.0`: Float, minimum purity needed for a split 
- `n_subfeatures = 0`Int,  number of features to select at random ( 0, keep all)
- `do_pruning = true` perform or not a post inference impurity pruning
- `pruning_accuracy = 1.0`:Float, purity threshold used for post-pruning (1.0, no pruning)
- `seed = 3`: the random seed
- `do_cross_validation = false`:Bool, do or not the n-fold cross validation
- `n_folds_cv = 3`: Int,   n-fold of the cross validation

# Output:

"""
function downstream_decision_tree_regression(kimchi_results::Matrix{Any}, # output of kimchi results
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

  names_of_the_wells_res = kimchi_results[2, 1:end]
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
  
  output = convert.(Float64, kimchi_results[row_to_learn, index_res])

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
    kimchi_results,
    feature_matrix,
    row_to_learn;
    options = SymbolicRegression.Options(),
    )


Function that evalauates 

# Arguments:
- `kimchi_results::Matrix{Any}`: The matrix of results of fitting one file with Kimchi. Compatible functions `fit_file_ODE`,`fit_file_custom_ODE`, `ODE_model_selection_file`, `segmentation_ODE_file`, `fit_NL_model_file`, `fit_NL_model_selection_file`, and `fit_NL_segmentation_file`.
- `feature_matrix::Matrix{Any}`: Matrix of the features for the ML analysis. Important, the number of rows of this file should be te number of columns (minus one) of the kimchi_results, and the first column should have a name of the well in order to mach the feature with the names of the well in second row of kimchi_results. 
- `row_to_learn::Int`: which row of the matrix `kimchi_results` will be the target of the ML inference.
# Key Arguments:

-  'options = SymbolicRegression.Options()' the option class of the symbolic regression class, see example and https://astroautomata.com/SymbolicRegression.jl/stable/api/#SymbolicRegression.CoreModule.OptionsStructModule.Options for details.

# Outputs:
if `res =  downstream_symbolic_regression()`:
-`trees`: the trees representing the hall of fames results
-`res_output`: a matrix containing the hall_of_fame  of the inference, where first column is the Complexity score, second column MSE and third column the equation Equation
-`predictions`: For each equation we return the predicted value for each sample (in this matrix equations are columns and rows are the samples)
-`index_annotation`: Index on how to order the rows of features matrix to match the columns of kimchi results

"""
function downstream_symbolic_regression2(kimchi_results,
  feature_matrix,
  row_to_learn;
  options = SymbolicRegression.Options(),
  n_processor = 2,
  niterations = 5
)

#Distributed.addprocs(n_processor)



  names_of_the_wells_res = kimchi_results[2, 1:end]
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
  
  output = convert.(Float64, kimchi_results[row_to_learn, index_res])

  
  predictors = Matrix(transpose(  convert.(Float64, feature_matrix[index_annotation,2])))
  predictors = reshape(predictors, length(predictors), 1)

 # predictors = reduce(vcat, predictors)))
 # hall_of_fame = SymbolicRegression.equation_search(predictors,output)
  hall_of_fame = RunSR(predictors,output,niterations, options)
  dominating = calculateParetoFrontier(predictors,output,hall_of_fame,options)



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






"""
    downstream_symbolic_regression(
    kimchi_results,
    feature_matrix,
    row_to_learn;
    options = SymbolicRegression.Options(),
    )


Function that evalauates 

# Arguments:
- `kimchi_results::Matrix{Any}`: The matrix of results of fitting one file with Kimchi. Compatible functions `fit_file_ODE`,`fit_file_custom_ODE`, `ODE_model_selection_file`, `segmentation_ODE_file`, `fit_NL_model_file`, `fit_NL_model_selection_file`, and `fit_NL_segmentation_file`.
- `feature_matrix::Matrix{Any}`: Matrix of the features for the ML analysis. Important, the number of rows of this file should be te number of columns (minus one) of the kimchi_results, and the first column should have a name of the well in order to mach the feature with the names of the well in second row of kimchi_results. 
- `row_to_learn::Int`: which row of the matrix `kimchi_results` will be the target of the ML inference.
# Key Arguments:

-  'options = SymbolicRegression.Options()' the option class of the symbolic regression class, see example and https://astroautomata.com/SymbolicRegression.jl/stable/api/#SymbolicRegression.CoreModule.OptionsStructModule.Options for details.

# Outputs:
if `res =  downstream_symbolic_regression()`:
-`trees`: the trees representing the hall of fames results
-`res_output`: a matrix containing the hall_of_fame  of the inference, where first column is the Complexity score, second column MSE and third column the equation Equation
-`predictions`: For each equation we return the predicted value for each sample (in this matrix equations are columns and rows are the samples)
-`index_annotation`: Index on how to order the rows of features matrix to match the columns of kimchi results

"""
function downstream_symbolic_regression(kimchi_results,
  feature_matrix,
  row_to_learn;
  options = SymbolicRegression.Options(),
)





  names_of_the_wells_res = kimchi_results[2, 1:end]
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
  
  output = convert.(Float64, kimchi_results[row_to_learn, index_res])

  
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
