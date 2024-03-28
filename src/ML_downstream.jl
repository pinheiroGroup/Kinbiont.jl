
using MLJ
using SymbolicRegression
using Zygote
using DelimitedFiles
using Plots
using CSV
using BetaML
using DecisionTree
using MLJDecisionTreeInterface




function downstream_decision_tree(jmaki_results,
  feature_matrix,
  row_to_learn;
  average_replicate=false,
  save_to_file=false,
  output_folder="/res/",
  output_file="output.txt",
  max_depth = 3,
  verbose = true
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

  index_res[1, :] = index_res[1, :] .+ 1
  output = convert.(Float64, jmaki_results[row_to_learn, index_res[1, :]])

  predictors = convert.(Float64, feature_matrix[index_annotation[1, :], 2:end])





  model = build_tree(output, predictors)
  imp_1 = impurity_importance(model)
  imp_2 = split_importance(model)

  if verbose == true
    tree_out = DecisionTree.print_tree(model, max_depth)
  end
    return model, imp_1, imp_2

end

function downstream_symbolic_regression(jmaki_results,
  feature_matrix,
  row_to_learn;
  average_replicate=false,
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
  should_optimize_constants=false)
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


jm_res_test = readdlm("/Users/fabrizio.angaroni/Documents/res_chem_isolate/res_to_test_ML_ode.csv", ',')
annotation_test = CSV.File("/Users/fabrizio.angaroni/Documents/JMAKi_utilities/real_dataset_tests/dataset/Monod_AA_detection/exp_4/annotation.csv")
names_of_annotation = propertynames(annotation_test)
feature_matrix = hcat(annotation_test[:V1], annotation_test[:V3])
jmaki_results = jm_res_test


# regression on gr

A = downstream_symbolic_regression(jmaki_results,
  feature_matrix,
  8;
  do_plots=true,)


# regression on N max

A2 = downstream_symbolic_regression(jmaki_results,
  feature_matrix,
  6;
  do_plots=true,
)

A3 = downstream_symbolic_regression(jmaki_results,
  feature_matrix,
  5;
  do_plots=true,
)




# testing decision trees on big dataset





jm_res_test = readdlm("/Users/fabrizio.angaroni/Documents/res_chem_isolate/res_to_test_ML_ode.csv", ',')
annotation_test = readdlm("/Users/fabrizio.angaroni/Documents/res_chem_isolate/annotation_to_test_ML_ode.csv", ',')
feature_matrix = annotation_test[:,2:(end-2)]
jmaki_results = jm_res_test



a = downstream_decision_tree(jmaki_results,
feature_matrix,
7;
average_replicate=false,
save_to_file=false,
output_folder="/res/",
output_file="output.txt",
max_depth = 3,
)

a = downstream_decision_tree(jmaki_results,
feature_matrix,
9;
average_replicate=false,
save_to_file=false,
output_folder="/res/",
output_file="output.txt",
max_depth = 3,
)


index_coli = findall(annotation_test[:,end].== "E. coli")
feature_matrix = annotation_test[index_coli,2:(end-2)]
jmaki_results = jm_res_test[:,index_coli]

a = downstream_decision_tree(jmaki_results,
feature_matrix,
7;
average_replicate=false,
save_to_file=false,
output_folder="/res/",
output_file="output.txt",
max_depth = 3,
)

a = downstream_decision_tree(jmaki_results,
feature_matrix,
9;
average_replicate=false,
save_to_file=false,
output_folder="/res/",
output_file="output.txt",
max_depth = 3,
)




index_mixture = findall(annotation_test[:,end].== "mixture")
feature_matrix = annotation_test[index_mixture,2:(end-2)]
jmaki_results = jm_res_test[:,index_mixture]

a = downstream_decision_tree(jmaki_results,
feature_matrix,
7;
average_replicate=false,
save_to_file=false,
output_folder="/res/",
output_file="output.txt",
max_depth = 3,
)

a = downstream_decision_tree(jmaki_results,
feature_matrix,
9;
average_replicate=false,
save_to_file=false,
output_folder="/res/",
output_file="output.txt",
max_depth = 3,
)