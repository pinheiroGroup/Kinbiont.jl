# [Data formatting and outputs](@id data)
```@contents
Pages = ["index.md"]
Depth = 4
```

## Data and annotation formatting for fitting one dimensional data


Kinbiont can operate directly on data files or within a Julia notebook. The format of a single time series must be a `2 x n_time_points` matrix of `Float64`:
```
 0.0        2.0       4.0       6.0       8.0        10.0       10.0       12.0       14.0       16.0       18.0       20.0       22.0      24.0      26.0       28.0       30.0       32.0       34.0       36.0       …  
 0.0912154  0.107956  0.105468  0.101727  0.0931484   0.106318   0.103697   0.139821   0.173598   0.204888   0.251052   0.289018   0.31298   0.33752   0.359356   0.370861   0.376347   0.383732   0.398496   0.384511 …  
```
The first row should represent the time, and the second row should represent the quantity to be fitted (e.g., optical density or CFU).

If the user calls APIs that require a `.csv` input, they must provide Kinbiont.jl with the paths to the `.csv` data file and the optional `.csv` annotation file. In these cases, Kinbiont expects a data matrix where the first row contains the names of the wells, and the other columns contain the numerical values of the measurements:
```
Time,  A1,     A2,      A3, 
0.0,   0.09,   0.09,    0.087,
1.0,   0.08,   0.011,   0.012,
2.0,   0.011,  0.18,    0.1,
3.0,   0.012,  0.32,    0.22,
4.0,   0.008,  0.41,    0.122,
```

Kinbiont.jl expects a comma (,) as the separator between columns and  the first column will be used as time.

The annotation file is optional (but mandatory if blank subtraction and averaging of replicates are required) and should be a two-column .csv file where the number of rows corresponds to the number of wells. The first column should contain the name of the well (they should match the names of the wells in the data .csv), while the second column should contain a unique ID for each technical replicate. A `b` indicates that the well should be considered as a blank, an `X` indicates that the well should be discarded from the analysis, and if two wells have the same ID, they will be considered replicates:
```
A1, b
A2, X
A3, unique_ID
```

To provide a calibration curve of optical density (OD), that maps OD values obtained from a microplate reader to corresponding values obtained from an independent source, the file should be provided to KinBiont as a CSV file containing two columns:
- `Raw_OD`: Optical density values measured using a microplate reader.
- `Real_OD`: Optical density values measured using an independent source.

```csv
Raw_OD,Real_OD
1.9617333333333333,4.19666666666
1.57826666,2.813333333
1.1751333333,1.68333333333
0.87273,1.005
0.66826666666,0.74533333333
0.45426666666,0.492
0.2812,0.283
0.09426666,0.097
0.04726666,0.04933333333333334
0.0238,0.024
0.0,0.0
```


See the folder  `data_examples` for examples. 

## Data and Annotation Formatting for Fitting Multidimensional Data  

If the user wants to use the functions `fit_ODEs_System`, `RN_fit`, and `fit_Cybernetic_models`, the format of a single time series must be an `(m+1) \times n_{\text{time points}}` matrix of `Float64`, where \( m \) is the number of equations present in the system (e.g., the number of chemical species).  

The first row is assumed to represent time, while the remaining rows correspond to the quantities being fitted. For example, if there are two equations, an input could be:  

```
 0.0        2.0       4.0       6.0       8.0        10.0       12.0       14.0       16.0       18.0       20.0       22.0   …  
 0.0912154  0.107956  0.105468  0.101727  0.0931484   0.106318   0.103697   0.139821   0.173598   0.204888   0.251052   0.289018   0.31298  …  
 0.1        0.2       0.3       0.4       0.5        0.6        0.7        0.8        0.9        0.9        0.9        0.9        0.9   …  
```

These functions also allow users to exclude specific data points (e.g., if the data for the first species is not available). In such cases, the user should specify `set_of_equation_to_fit=[2]` in the fitting functions and provide the data in the following format:  

```
 0.0        2.0       4.0       6.0       8.0        10.0       12.0       14.0       16.0       18.0       20.0       22.0   …  
 0.1        0.2       0.3       0.4       0.5        0.6        0.7        0.8        0.9        0.9        0.9        0.9        0.9   …  
```



## Data and annotation formatting for downstream ML


All ML functions of Kinbiont take as input a matrix of results  (i.e., the outputs of fits `Kinbiont_results`) and a matrix of features  (e.g., the concentration of antibiotics present in all wells `feature_matrix`).
For example:
```julia
downstream_decision_tree_regression(Kinbiont_results, 
  feature_matrix,
  row_to_learn;
)
```

```julia
downstream_symbolic_regression(Kinbiont_results,
  feature_matrix,
  row_to_learn;
)
```

The first matrix is the standard output of any of the Kinbiont fit functions. In this case, each row represents a parameter, and each column a growth curve.

For example:

```
label_exp,       exp_2_no_corrections,    exp_2_no_corrections
well,            A1,                      A2
model,           HPM,                     HPM
gr,              0.008875566468779583,    0.010090369128600398
exit_lag_rate,   1.7249775833759684e-6,   1.4012949810152472e-6
N_max,           2.498999749717784,       1.6986904454789507
th_max_gr,       0.005778245042794245,    0.00599548534261212
emp_max_gr,      0.007951369027199616,    0.008096305651156249
loss,            0.0013005418069932683,   0.0013349159149782007
```

Note the first row is dedicated to the label of the experiment and will not be used by the functions. It is necessary that the second row (i.e., `well`) reports a unique ID for each curve. The functions will ask which is the target row of the regression (i.e., `row_to_learn`); please do not use the first row. The first column is dedicated to the names of the columns and will be discarded from the ML analysis.

Instead, the feature matrix specifies the conditions associated with each unique ID of the previous matrix. Only the wells where there is a match  between the first column of the feature matrix and the second row of the fitting results will be used. For example, suppose you have two different antibiotics, each with two different concentrations. Then the matrix could be:

```
ID_exp,   abx_1,   abx_2
A1,       0,       1,
A2,       2,       0,
A3,       1,       1,
A4,       1,       0,
A5,       0,       2,
```

Note that it is necessary to add one column for each new chemical/condition added to the experiment (even if in a specific well it is absent).  The first row will not be used and is dedicated to the feature names.

See the folder  `data_examples` for examples. 

## Outputs of Kinbiont

Kinbiont has different data structures as output. Here, you will find a brief summary of the main function; please consult the API section for the description of all the functions.

- `fitting_one_well_Log_Lin`

This structure stores results for a single well using a log-linear method.

1. `method:String` - The method used.
1. `params:Vector{Any}` - Parameters obtained from the fitting process.
1. `fit:Any` - The fitted function in the exponential window.
1. `times:Any` - The times at which measurements were taken.
1. `confidence_band:Any` - The confidence band of the fit.

- `fitting_one_well_ODE_constrained` and `fitting_one_well_custom_ODE`

This structure stores results for a single well.

1. `method:String` - The method used.
1. `params:Vector{Any}` - Parameters obtained from the fitting process.
1. `fit:Any` - The fitted function.
1. `times:Any` - The times at which measurements were taken.

- `ODE_Model_selection`

1. `method:String` - The method used.
1. `params:Vector{Any}` - Parameters obtained from the fitting process.
1. `fit:Any` - The fitted function.
1. `times:Any` - The times at which measurements were taken.

- `Kinbiont_res_bootstrap_NL`

This structure stores results of the bootstrap fitting of a NL function.

1. `method:String` - The method used.
1. `params:Matrix{Any}` - Parameters obtained from the fitting process.
1. `fit:Any` - The fitted function.
1. `times:Any` - The times at which measurements were taken.
1. `fin_param:Any` - The  parameters of each bootstrap fit.
1. `new_param_fin:Any` - The  parameters of each bootstrap fit after considering only the best $95/%$ of the losses.
1. `mean_param:Any` - Mean of the parameters.
1. `sd_param:Any` - Standard deviation of the parameters.

- `ODE_Model_selection`

This structure stores model selection results.

1. `method:String` - The method used.
1. `params:Vector{Any}` - Parameters obtained from the fitting process.
1. `fit:Vector{Float64}` - The fit result.
1. `times:Vector{Float64}` - The times at which measurements were taken.
1. `rss_array:Any` - Residual sum of squares array.
1. `min_rss_array:Any` - Minimum residual sum of squares array.
1. `param_min:Any` - Parameters corresponding to minimum RSS.
1. `min_AIC:Vector{Any}` - Minimum AIC values.
1. `selected_model:String` - The selected model.
1. `full_param:Vector{Any}` - Full parameter set.


- `NL_model_selection`

This structure stores non-linear model selection results.

1. `method:String` - The method used.
1. `params:Vector{Any}` - Parameters obtained from the fitting process.
1. `fit:Vector{Float64}` - The fit result.
1. `times:Vector{Float64}` - The times at which measurements were taken.
1. `score_res:Any` - Score results.
1. `top_loss:Any` - Top loss values.

- `sensitivity_NL`

This structure stores sensitivity analysis results using non-linear methods.

1. `method:String` - The method used.
1. `params:Matrix{Any}` - Parameters obtained from the fitting process.
1. `fit:Any` - The fit result.
1. `times:Any` - The times at which measurements were taken.
1. `combinations:Matrix{Any}` - The list of each of the starting hyperparameters  used in sensitivity analysis.

- `one_well_morris_sensitivity`

This structure stores sensitivity analysis results.

1. `method:String` - The method used.
1. `params:Matrix{Any}` - Parameters obtained from the fitting process.
1. `combinations:Matrix{Any}` - The list of each of the starting hyperparameters  used in sensitivity analysis.

- `segmentation_ODE`

This structure stores segmentation results using ODE methods.

1. `method:String` - The method used.
1. `params:Matrix{Any}` - Parameters obtained from the fitting process.
1. `fit:Array{Float64}` - The fitted functions.
1. `times:Array{Float64}` - The times at which measurements were taken.
1. `interval_cdp:Array{Any}` - Change point intervals.
1. `score_of_the_models:Any` - Scores of the models.

- `segmentation_NL`

This structure stores segmentation results using non-linear methods.

1. `method:String` - The method used.
1. `params:Matrix{Any}` - Parameters obtained from the fitting process.
1. `fit:Array{Float64}` - The fitted functions.
1. `times:Array{Float64}` - The times at which measurements were taken.
1. `interval_cdp:Array{Any}` - Change point intervals.

- `fit_ODEs_System`

1. `method:String` - The method used.
1. `params:Matrix{Any}` - Parameters obtained from the fitting process.
1. `fit:Array{Float64}` - The fitted functions in SciML format.

- `fit_Cybernetic_models`

1. `method:String` - The method used.
1. `params:Matrix{Any}` - Parameters obtained from the fitting process.
1. `fit:Array{Float64}` - The fitted functions in SciML format.

- `RN_fit`

1. `method:String` - The method used.
1. `params:Matrix{Any}` - Parameters obtained from the fitting process.
1. `fit:Array{Float64}` - The fitted functions in SciML format.


This structure stores results for log-linear fits across multiple curves in one file.

1. `method:String` - The method used.
1. `params:Matrix{Any}` - The matrix with the  parameters obtained from the fitting process of each curve of the file.
1. `fits:Tuple{Any}` - The fitted functions in each exponential window.
1. `data:Tuple{Any}` - The data used for fitting.
1. `confidence_bands:Tuple{Any}` - Confidence bands for the fits.

- `Kinbiont_res_one_file`

This structure stores results of the fit for all curves in a single file.

1. `method:String` - The method used.
1. `params:Matrix{Any}` - The matrix of the parameters obtained from the fitting process.
1. `fits:Tuple{Any}` - The fitted functions for each well.
1. `data:Tuple{Any}` - The data used for fitting.

- `Kinbiont_res_one_file_segmentation`

This structure stores segmentation results  for all curves in a single file.

1. `method:String` - The method used.
1. `params:Matrix{Any}` - The matrix of the parameters obtained from the fitting process.
1. `fits:Tuple{Any}` -  The fitted functions for each well.
1. `data:Tuple{Any}` - The data used for fitting.
1. `cp:Tuple{Any}` - Change points detected.
1. `vector_AIC:Any` - AIC (or AICc) values of the best model for each well.
 

 - `ODEs_system_fit`

This structure stores segmentation results  for all curves in a single file.

1. `method:String` - The method used.
1. `params:Matrix{Any}` - The matrix of the parameters obtained from the fitting process.
1. `fits:Tuple{Any}` -  The fitted functions for each well.
1. `data:Tuple{Any}` - The data used for fitting.
1. `cp:Tuple{Any}` - Change points detected.
1. `vector_AIC:Any` - AIC (or AICc) values of the best model for each well.
 
