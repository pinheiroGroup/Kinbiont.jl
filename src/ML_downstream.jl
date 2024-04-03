

### downstream decision tree

function downstream_decision_tree_regression(jmaki_results::Matrix{Any}, # output of jmaki results
  feature_matrix::Matrix{Any},
  row_to_learn::Int;
  max_depth = -1,
  verbose = true,
  pruning_purity = 1.0,
  min_samples_leaf = 5,
  min_samples_split = 2,  
  min_purity_increase = 0.0, 
  n_subfeatures = 0,
  seed = 3,
  do_cross_validation = false,
  n_folds_cv = 3,
)
max_depth = convert(Int, max_depth)



  names_of_the_wells_res = jmaki_results[3, 2:end]
  names_of_the_wells_annotation = feature_matrix[1:end, 1]
  wells_to_use = intersect(names_of_the_wells_res, names_of_the_wells_annotation)


  index_res = Int
  index_annotation = Int




  for i in wells_to_use
    if i == wells_to_use[1]
      index_res = findfirst(names_of_the_wells_res .== i)
      index_annotation = findfirst(names_of_the_wells_annotation .== i)

    else

      index_res = hcat(index_res, findfirst(names_of_the_wells_res .== i))
      index_annotation = hcat(index_annotation, findfirst(names_of_the_wells_annotation .== i))

    end


  end

  # order results and annotation by well in the same order

  index_res[1, :] = index_res[1, :] .+ 1
  output = convert.(Float64, jmaki_results[row_to_learn, index_res[1, :]])

  predictors = convert.(Float64, feature_matrix[index_annotation[1, :], 2:end])





  model = build_tree(output, predictors)
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
  if verbose == true
    tree_out = DecisionTree.print_tree(model, max_depth)
  end
    return model, imp_1, imp_2, r2

end

function downstream_symbolic_regression(jmaki_results,
  feature_matrix,
  row_to_learn;
  do_plots=true,
  binary_operators=[+, *, /],
  unary_operators=nothing,
  batching=false,
  batch_size=1000,
  save_to_file=false,
  output_folder="/res/",
  output_file="output.txt",
  niterations=100,
  enable_autodiff=true,
  should_simplify=true,
  should_optimize_constants=false,
)



  if save_to_file == true
    mkpath(output_folder)
  end
  names_of_the_wells_res = jmaki_results[3, 2:end]
  names_of_the_wells_annotation = feature_matrix[1:end, 1]
  wells_to_use = intersect(names_of_the_wells_res, names_of_the_wells_annotation)


  index_res = Int
  index_annotation = Int




  for i in wells_to_use
    if i == wells_to_use[1]
      index_res = findfirst(names_of_the_wells_res .== i)
      index_annotation = findfirst(names_of_the_wells_annotation .== i)

    else

      index_res = hcat(index_res, findfirst(names_of_the_wells_res .== i))
      index_annotation = hcat(index_annotation, findfirst(names_of_the_wells_annotation .== i))

    end


  end

  # order results and annotation by well in the same order






  ML_model = SRRegressor(
    niterations=niterations,
    binary_operators=binary_operators,
    # unary_operators=unary_operators,
    should_simplify=should_simplify,
    should_optimize_constants=should_optimize_constants,
    batching=batching,
    batch_size=batch_size,
    output_file=output_file,
    save_to_file=save_to_file
  )


  index_res[1, :] = index_res[1, :] .+ 1
  output = convert.(Float64, jmaki_results[row_to_learn, index_res[1, :]])

  predictors = convert.(Float64, feature_matrix[index_annotation[1, :], 2:end])




  mach = MLJ.machine(ML_model, predictors, output) |> fit!

  # MLJ.fit!(mach)
  # report the fit 
  res_gr = report(mach)
  res_1 = res_gr.equation_strings
  res_2 = res_gr.equations
  res_3 = res_gr.equations[res_gr.best_idx]
  res_4 = res_gr.losses
  res_5 = res_gr.complexities
  res_6 = res_gr.scores
  index_min_score = findall(res_6 .== minimum(res_gr.scores))


  if size(feature_matrix)[2] == 2 && do_plots == true



    display(scatter(predictors, output, xlabel="predictor", ylabel="feature", label=["Data" nothing]))
    for k in 1:length(res_gr.equations)

      a1 = MLJ.predict(mach, (data=predictors, idx=k))
      display(scatter!(predictors, a1, xlabel="predictor", ylabel="feature", label=[res_1[k] nothing]))
    end

  elseif size(feature_matrix)[2] != 2 & do_plots == true

    println("Warning: the plots are possible only for 2D data (1 feature, 1 output)")

  end


  return res_1, res_2, res_3, res_4, res_5, res_6, res_gr.equations[index_min_score]

end

