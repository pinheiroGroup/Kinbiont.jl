using Kimchi
using DelimitedFiles
using Plots
using AbstractTrees
using DecisionTree
using StatsPlots

Fit_res = Matrix(permutedims(readdlm("E:/Lavoro/JMAKi_utilities-main/real_dataset_tests/dataset/species-specific_activity_antibio/fit_results.csv", ',')))
features_matrix = readdlm("E:/Lavoro/JMAKi_utilities-main/real_dataset_tests/dataset/species-specific_activity_antibio/features_matrix.csv", ',')


#selecting e colu bw
splitted_strain=split.(Fit_res[2,:], "_")

list_strains =[splitted_strain[i][2] for i in eachindex(splitted_strain)]
splitted_strain2=split.(list_strains, "-")
list_strains =[splitted_strain2[i][1] for i in eachindex(splitted_strain2)]
last_BW = 100#findlast(list_strains.=="BW")
# removing no drug and single drug from df
features_matrix_bw =features_matrix[1:last_BW,2:(end-2)]
Fit_res_bw =  Fit_res[:,1:last_BW]



histogram(Fit_res_bw[6,2:end], xlabel = "lag [h]",ylabel="counts" )

@time dt_inference_lag = downstream_decision_tree_regression(Fit_res_bw, # output of jmaki results
features_matrix_bw,
6;
max_depth = -1,
verbose = true,
pruning_purity = 1.00,
min_samples_leaf = 5,
min_samples_split = 2,  
min_purity_increase = 0.0, 
n_subfeatures = 0,
do_pruning = false,
pruning_accuracy = 0.0,
seed = 3,
do_cross_validation = true,
n_folds_cv = 10,
)

CSV.write("E:/Lavoro/JMAKi_utilities-main/real_dataset_tests/dataset/species-specific_activity_antibio/impurity_lag.csv",dt_inference_lag[2])
CSV.write("E:/Lavoro/JMAKi_utilities-main/real_dataset_tests/dataset/species-specific_activity_antibio/rank_lag.csv",dt_inference_lag[3])

dt_inference_gr = downstream_decision_tree_regression(Fit_res_bw, # output of jmaki results
features_matrix_bw,
10;
max_depth = -1,
verbose = true,
pruning_purity = 1.0,
min_samples_leaf = 5,
min_samples_split = 2,  
min_purity_increase = 0.0, 
n_subfeatures = 0,
do_pruning = false,
pruning_accuracy = 1.0,
seed = 3,
do_cross_validation = true,
n_folds_cv = 10,
)
CSV.write("E:/Lavoro/JMAKi_utilities-main/real_dataset_tests/dataset/species-specific_activity_antibio/impurity_gr.csv",dt_inference_gr[2])
CSV.write("E:/Lavoro/JMAKi_utilities-main/real_dataset_tests/dataset/species-specific_activity_antibio/rank_gr.csv",dt_inference_gr[3])

# creating delta N
Fit_res_bw[8,2:end] = Fit_res_bw[8,2:end] .- Fit_res_bw[7,2:end]
histogram(Fit_res_bw[8,2:end], xlabel = "delta N [1/h]",ylabel="counts" )


dt_inference_delta_N = downstream_decision_tree_regression(Fit_res_bw, # output of jmaki results
features_matrix_bw,
8;
max_depth = -1,
verbose = true,
pruning_purity = 1.0,
min_samples_leaf = 5,
min_samples_split = 2,  
min_purity_increase = 0.0, 
n_subfeatures = 0,
do_pruning = false,
pruning_accuracy = 1.0,
seed = 3,
do_cross_validation = true,
n_folds_cv = 10,
)

CSV.write("E:/Lavoro/JMAKi_utilities-main/real_dataset_tests/dataset/species-specific_activity_antibio/impurity_N.csv",dt_inference_delta_N[2])
CSV.write("E:/Lavoro/JMAKi_utilities-main/real_dataset_tests/dataset/species-specific_activity_antibio/rank_N.csv",dt_inference_delta_N[3])
