# [Data formatting and outputs](@id data)
## Data and annotation formatting

KinBiont can operate directly on data files or inside the julia notebook.
When are in a julia notebook the  format of single time series that want to be analyzed is a 2 x n_time_points Matrix of FLoat64, e.g.,


```
 0.0        2.0       4.0       6.0       8.0        10.0       10.0       12.0       14.0       16.0       18.0       20.0       22.0      24.0      26.0       28.0       30.0       32.0       34.0       36.0       …  
 0.0912154  0.107956  0.105468  0.101727  0.0931484   0.106318   0.103697   0.139821   0.173598   0.204888   0.251052   0.289018   0.31298   0.33752   0.359356   0.370861   0.376347   0.383732   0.398496   0.384511 …  

```
The first row should be time and the second the quantity to be fitted (e.g., Optical Density or CFU)

Instead, three APIs call direclty the files: the user must input  the paths to  a .csv data file and a .csv annotation to the functions of KinBiont.jl
; In these cases KinBiont expect for data a matrix where the first row are the names of the wells and the columns the numerical value of the measurements. Note that the first one will be used as time:

```
Time,  A1,     A2,      A3, 
0.0,   0.09,   0.09,    0.087,
1.0,   0.08,   0.011,   0.012,
2.0,   0.011,  0.18,    0.1,
3.0,   0.012,  0.32,    0.22,
4.0,   0.008,  0.41,    0.122,
```
KinBiont expect a "," as separator between columns

The annotation file instead should be a two columns .csv file where the number of rows correspond to the number of wells, note that the name of the well should be the same between the data.csv and annotation.csv:

```
A1, b
A2, X
A3, unique_ID

```
as unique_ID the user can insert anything but consider that if two wells has the same ID the will be considered replicates. 'b' indicates that the well should be cosidered a blank and 'X' that the well should be discarded from any analysis






To provide a calibration curve of OD, that maps optical density (OD) values obtained from a microplate reader to corresponding values obtained from a spectrophotometer, the file should be provided to KinBiont as a CSV file contains two columns:
- `Raw_Microplate_OD`: Optical density values measured using a microplate reader.
- `Real_Spectrophotometer_OD`: Optical density values measured using a real spectrophotometer.

```csv
Raw_Microplate_OD,Real_Spectrophotometer_OD
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


See the folders  XXXXX for examples. 

## Data and annotation formatting for downstream ML


## Outputs of KinBiont

KinBiont has different data struct as output

- `KinBiont_res_one_well_log_lin`

This structure stores results for a single well using a log-linear method.

1. `method:String` - The method used.
1. `params:Vector{Any}` - Parameters obtained from the fitting process.
1. `fit:Any` - The fit result.
1. `times:Any` - The times at which measurements were taken.
1. `confidence_band:Any` - The confidence band for the fit.

- `KinBiont_res_one_well`

This structure stores results for a single well.

1. `method:String` - The method used.
1. `params:Vector{Any}` - Parameters obtained from the fitting process.
1. `fit:Any` - The fit result.
1. `times:Any` - The times at which measurements were taken.

- `KinBiont_res_bootstrap_NL`

This structure stores results from a bootstrap process using non-linear methods.

1. `method:String` - The method used.
1. `params:Matrix{Any}` - Parameters obtained from the fitting process.
1. `fit:Any` - The fit result.
1. `times:Any` - The times at which measurements were taken.
1. `fin_param:Any` - Final parameters after bootstrapping.
1. `new_param_fin:Any` - New final parameters.
1. `mean_param:Any` - Mean of the parameters.
1. `sd_param:Any` - Standard deviation of the parameters.

- `KinBiont_res_Log_Lin_files`

This structure stores results for log-linear fits across multiple files.

1. `method:String` - The method used.
1. `params:Matrix{Any}` - Parameters obtained from the fitting process.
1. `fits:Tuple{Any}` - The fit results.
1. `data:Tuple{Any}` - The data used for fitting.
1. `confidence_bands:Tuple{Any}` - Confidence bands for the fits.

- `KinBiont_res_one_file`

This structure stores results for a single file.

1. `method:String` - The method used.
1. `params:Matrix{Any}` - Parameters obtained from the fitting process.
1. `fits:Tuple{Any}` - The fit results.
1. `data:Tuple{Any}` - The data used for fitting.

- `KinBiont_res_one_file_segmentation`

This structure stores segmentation results for a single file.

1. `method:String` - The method used.
1. `params:Matrix{Any}` - Parameters obtained from the fitting process.
1. `fits:Tuple{Any}` - The fit results.
1. `data:Tuple{Any}` - The data used for fitting.
1. `cp:Tuple{Any}` - Change points detected.
1. `vector_AIC:Any` - AIC values for model selection.

- `KinBiont_res_model_selection`

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


- `KinBiont_res_NL_model_selection`

This structure stores non-linear model selection results.

1. `method:String` - The method used.
1. `params:Vector{Any}` - Parameters obtained from the fitting process.
1. `fit:Vector{Float64}` - The fit result.
1. `times:Vector{Float64}` - The times at which measurements were taken.
1. `score_res:Any` - Score results.
1. `top_loss:Any` - Top loss values.

- `KinBiont_res_sensitivity_NL`

This structure stores sensitivity analysis results using non-linear methods.

1. `method:String` - The method used.
1. `params:Matrix{Any}` - Parameters obtained from the fitting process.
1. `fit:Any` - The fit result.
1. `times:Any` - The times at which measurements were taken.
1. `combinations:Matrix{Any}` - Combinations of parameters used in sensitivity analysis.

- `KinBiont_res_sensitivity`

This structure stores sensitivity analysis results.

1. `method:String` - The method used.
1. `params:Matrix{Any}` - Parameters obtained from the fitting process.
1. `combinations:Matrix{Any}` - Combinations of parameters used in sensitivity analysis.

- `KinBiont_res_segmentation_ODE`

This structure stores segmentation results using ODE methods.

1. `method:String` - The method used.
1. `params:Matrix{Any}` - Parameters obtained from the fitting process.
1. `fit:Array{Float64}` - The fit result.
1. `times:Array{Float64}` - The times at which measurements were taken.
1. `interval_cdp:Array{Any}` - Change point intervals.
1. `score_of_the_models:Any` - Scores of the models.

- `KinBiont_res_segmentation_NL`

This structure stores segmentation results using non-linear methods.

1. `method:String` - The method used.
1. `params:Matrix{Any}` - Parameters obtained from the fitting process.
1. `fit:Array{Float64}` - The fit result.
1. `times:Array{Float64}` - The times at which measurements were taken.
1. `interval_cdp:Array{Any}` - Change point intervals.
