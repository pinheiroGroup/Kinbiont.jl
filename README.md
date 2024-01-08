<p align="center">
  <img width="300" height="800" src="https://github.com/ang-one/JMAKi.jl/blob/main/static/jmaki_logo.png">
</p>
# JMAKi.jl

JMAKi, a versatile software tool that utilizes Ordinary Differential Equations (ODEs) to fit bacterial growth data from plate reader experiments. 
With JMAKi it is possible to simulate, fit, perform model selection, and conduct sensitivity analysis for multi-well plate reader experiments.
The parameter fitting in JMAKi is defined as a constrained optimization problem, which is solved using a differential evolution non-linear optimizer.

To address complex cases,  JMAKi preprocesses data by detecting change points in the differential operator within the time series. 
It then automatically assembles and fits a segmented ODE model, resulting in a fully interpretable representation of the time series.

1. [Installation](#installation)
2. [Data and annotation formatting](#data)
3. [API](#API)
4. [The mathematical models](#models)
5. [Examples and Tutorial](#examples)

<a name="installation"></a>
# Installation & requirements
## Manual installation
Download the folder from Github. 
First, it is necessary to install the Julia from https://julialang.org/.   
Next, the user need to copy the project folder in the chosen working directory. 

1. Using REPL or the COMMAND LINE move to the working directory.  
2. If you use the COMMAND LINE, to start a Julia session run the command:

> julia

3. To enter in the Pkg REPL  type 

>]  

4. Type the command 
> activate .

5. To activate the JMAKi project, type
> instantiate

6. at the start of the code or notebook (where you are going to do analyses) you should write 

```
using DifferentialEquations, Optimization, Plots, Random, CSV,  DataFrames, Statistics, Optim, OptimizationBBO, NaNMath, CurveFit, StatsBase, Tables, Distributions,Interpolations,GaussianProcesses,Peaks,ChangePointDetection

include("your_path_to_JMAKi_main_folder/src/functions.jl")


```
this last step is Temporary before the official realese
## Package installation
## Requirements
### Dependencies

1. Julia (1.7,1.8,1.9)
2. DifferentialEquations.jl
3. Optimization.jl
4. Plots.jl
5. Random.jl
6. CSV.jl
7. DataFrames.jl
8. Statistics.jl
9. Optim.jl
10. OptimizationBBO.jl
11. NaNMath.jl
12. CurveFit.jl
13. StatsBase.jl
14. Tables.jl
15. Distributions.jl
16. Interpolations.jl
17. GaussianProcesses.jl 
18. Peaks.jl
19. ChangePointDetection.jl


<a name="data"></a>
# Data and annotation formatting
JMAKi can operate directly on data files or inside the julia notebook.
When are in a julia notebook the  format of single time series that want to be analyzed is a 2 x n_time_points Matrix of FLoat64, e.g.,


```
 0.0        2.0       4.0       6.0       8.0        10.0       10.0       12.0       14.0       16.0       18.0       20.0       22.0      24.0      26.0       28.0       30.0       32.0       34.0       36.0       …  
 0.0912154  0.107956  0.105468  0.101727  0.0931484   0.106318   0.103697   0.139821   0.173598   0.204888   0.251052   0.289018   0.31298   0.33752   0.359356   0.370861   0.376347   0.383732   0.398496   0.384511 …  

```
The first row should be time and the second the quantity to be fitted (e.g., Optical Density or CFU)

Instead, three APIs call direclty the files: the user must input  the paths to  a .csv data file and a .csv annotation to the functions of JMAKi.jl
; In these cases JMAKi expect for data a matrix where the first row are the names of the wells and the columns the numerical value of the measurements. Note that the first one will be used as time:

```
Time,  A1,     A2,      A3, 
0.0,   0.09,   0.09,    0.087,
1.0,   0.08,   0.011,   0.012,
2.0,   0.011,  0.18,    0.1,
3.0,   0.012,  0.32,    0.22,
4.0,   0.008,  0.41,    0.122,
```
JMAKi expect a "," as separator between columns

The annotation file instead should be a two columns .csv file where the number of rows correspond to the number of wells, note that the name of the well should be the same between the data.csv and annotation.csv:

```
A1, b
A2, X
A3, unique_ID

```
as unique_ID the user can insert anything but consider that if two wells has the same ID the will be considered replicates. 'b' indicates that the well should be cosidered a blank and 'X' that the well should be discarded from any analysis


See the folders  XXXXX for some examples. 


If a OD calibration curved is provided it should have the following format XXXXXXX

<a name="API"></a>
# The main functions of JMAKi
1. [Simulate ODE](#simulating-ODE)
2. [Stochastic simulation](#simulating-stochastic)
3. [Plotting a dataset from file](#plot-file)
4. [Specific growth rate evaluation](#Specific_growth_rate_evaluation)
5. [Fitting growth rate with log-lin fitting for one well](#log-lin-one-well)
6. [Fitting growth rate with log-lin fitting for one well](#log-lin-file)
7. [Fitting ODE function for one well](#ODE-one-well)
8. [Fitting ODE function for one file](#ODE-file)
9. [Fitting custom ODE function](#custom-ODE)
10. [Sensitivity analysis](#Sensitivity-analysis)
11. [Model selection](#Model-selection)
12. [Change points detection](#cdp)
13. [Fitting segmented ODE with fixed change-point number](#cdp-fixed)
14. [Fitting segmented ODE with direct search for a maximum number of change-points](#cdp-search)


<a name="simulating-ODE"></a>
## Simulate ODE
```

   ODE_sim(model::String, 
    n_start::Vector{Float64}, 
    tstart::Float64, 
    tmax::Float64, 
    delta_t::Float64, 
    param_of_ode::Vector{Float64}; 
    integrator = KenCarp4() 
)
```
This function performs an ODE simulation of a model, considering the initial conditions, time range, and integration parameters.

Arguments:
- `model::String`: The model to simulate.
- `n_start::Vector{Float64}`: The starting conditions.
- `tstart::Float64`: The start time of the simulation.
- `tmax::Float64`: The final time of the simulation.
- `delta_t::Float64`: The time step of the output.
- `param_of_ode::Vector{Float64}`: The parameters of the ODE model.

  
Key argument:
- `integrator=KenCarp4() `: The chosen solver from the SciML ecosystem for ODE integration, default KenCarp4 algorithm.

<a name="simulating-stochastic"></a>
##  Stochastic simulation
```
    stochastic_sim(model::String,
         n_start::Int,
         n_mol_start::Float64,
        tstart::Float64,
        tmax::Float64,
        delta_t::Float64,
        k_1_val::Float64,
        k_2_val::Float64,
        alpha_val::Float64,
        lambda::Float64,
        n_mol_per_birth::Float64,
        volume::Float64)
```
This function performs a stochastic simulation of a model, considering cell growth and nutrient consumption over time.

Arguments:

- `model::String`: The model to simulate. PUT the options
- `n_start::Int`: The number of starting cells.
- `n_mol_start::Float64`: The starting concentration of the limiting nutrient.
- `tstart::Float64`: The start time of the simulation.
- `tmax::Float64`: The final time of the simulation.
- `delta_t::Float64`: The time step for the Poisson approximation.
- `k_1_val::Float64`: The value of parameter k1.
- `k_2_val::Float64`: The value of the Monod constant.
- `alpha_val::Float64`: The maximum possible growth rate.
- `lambda::Float64`: The lag time.
- `n_mol_per_birth::Float64`: The nutrient consumed per division (mass).
- `volume::Float64`: The volume.

<a name="plot-file"></a>
## Plotting a dataset from file
```
plot_data( label_exp::String, 
    path_to_data::String, 
    path_to_annotation::String;
    path_to_plot="NA", 
    display_plots=true ,
    save_plot=false, 
    overlay_plots=true, 
    blank_subtraction="NO", 
    average_replicate=false, 
    correct_negative="thr_correction", 
    thr_negative=0.01,
    )
    
```

This function plot all the data from .csv file
Arguments:

- `path_to_data::String`: The path to the .csv of data
-  `path_to_annotation::String`: The path to the .csv of annotation 
- `name_well::String`: The name of the well.
- `label_exp::String`: The label of the experiment.

Key Arguments:
- `path_to_plot= "NA"`: Path to save the plots.
-  `save_plot=false` : save the plot or not
- ` display_plots=true`: Whether or not diplay the plot in julia
-    `overlay_plots =true` : true on plot for all dataset false one plot per well
- `verbose=false`: Whether to enable verbose output.
- `pt_avg=7`: Number of points to use for smoothing average.
- ` blank_subtraction="NO"`: 
- ` average_replicate=false`
- `multiple_scattering_correction=false`: Whether or not correct the data qith a calibration curve.
- `calibration_OD_curve="NA"`: The path where the .csv calibration data are located, used only if `multiple_scattering_correction=true`.




<a name="Specific-growth-rate-evaluation"></a>
## Specific growth rate evaluation
```
specific_gr_evaluation(data_smooted::Matrix{Float64},
     pt_smoothing_derivative::Int)
```


<a name="log-lin-one-well"></a>
## Fitting growth rate with log-lin fitting for one well
``` fitting_one_well_Log_Lin(data::Matrix{Float64}, 
    name_well::String, 
    label_exp::String; 
    do_plot=false,
    path_to_plot="NA", 
    type_of_smoothing="rolling_avg",
    pt_avg=7, 
    pt_smoothing_derivative=7, 
    pt_min_size_of_win=7, 
    type_of_win="maximum", 
    threshold_of_exp=0.9, 
    multiple_scattering_correction=false,
    calibration_OD_curve ="NA" 
    )
```
This function fits a logarithmic-linear model to a single well's data and performs analysis such as plotting and error calculation.
Arguments:

- `data::Matrix{Float64}`: The dataset of OD/fluorescence values.
- `name_well::String`: The name of the well.
- `label_exp::String`: The label of the experiment.

Key Arguments:
- `do_plot=true`: Whether to generate and save plots.
- `path_to_plot="NA"`: The path to save the plots, used only if `do_plot=true`.
- `pt_avg=7`: The number of points to do rolling average smoothing.
- `pt_smoothing_derivative=7`: The number of points of the window to evaluate specific growth rate.
- `pt_min_size_of_win=7`: The minimum size of the exponential windows in the number of smoothed points.
- `type_of_win="maximum"`: How the exponential phase window is selected ("maximum" or "global_thr").
- `threshold_of_exp=0.9`: The threshold of the growth rate in quantile to define the exponential windows.
- `multiple_scattering_correction=false`: Whether or not correct the data qith a calibration curve.
- `calibration_OD_curve="NA"`: The path where the .csv calibration data are located, used only if `multiple_scattering_correction=true`.

<a name="log-lin-file"></a>
## Fitting growth rate with log-lin fitting for one file
```
    fit_one_file_Log_Lin(
    label_exp::String, 
    path_to_data::String,
    path_to_annotation::String;
    path_to_results = "NA",
    path_to_plot= "NA",
    do_plot=false, 
    verbose=false,
    write_res=false, 
    type_of_smoothing="rolling_avg", 
    pt_avg=7,
    pt_smoothing_derivative=7, 
    pt_min_size_of_win=7, 
    type_of_win="maximum", 
    threshold_of_exp=0.9, 
    blank_subtraction="avg_blank", 
    fit_replicate=false, 
    correct_negative="thr_correction",
    thr_negative=0.01, 
    multiple_scattering_correction=false, 
    calibration_OD_curve="NA" 
    )
```
This function fits a logarithmic-linear model to a single file's data. It performs model fitting, error analysis, and provides various options for customization.
Arguments:

- `label_exp::String`: Label of the experiment.
- `path_to_data::String`: Path to the folder containing the data.
- `path_to_annotation::String`: Path to the annotation of the wells.

Key Arguments:



- `path_to_results= "NA"`: Path to save the results.
- `path_to_plot= "NA"`: Path to save the plots.
- `do_plot=false`: Whether to generate and visualize plots of the data.
- `verbose=false`: Whether to enable verbose output.
- `write_res= false`: Whether to write results.
- `pt_avg=7`: Number of points to use for smoothing average.
- `pt_smoothing_derivative=7`: The number of points of the window to evaluate specific growth rate.
- `pt_min_size_of_win=7`: Minimum size of the exponential windows in number of smoothed points.
- `type_of_win="maximum`: How the exponential phase window is selected ("maximum" or "global_thr").
- `threshold_of_exp=0.9`: Threshold of growth rate in quantile to define the exponential windows.
- `blank_subtraction="avg_blank"`: How to use blank data for subtraction (options: "NO", "avg_subtraction", "time_avg").
- `fit_replicate=false`: If `true`, fit the average between replicates; if `false`, fit all replicates independently.
- `correct_negative="thr_correction`: Method to correct negative values (options: "thr_correction", "blank_correction").
- `thr_negative=0.01`: Threshold value used only if `correct_negative == "thr_correction"`.
- `multiple_scattering_correction=false`: Whether or not correct the data qith a calibration curve.
- `calibration_OD_curve="NA"`: The path where the .csv calibration data are located, used only if `multiple_scattering_correction=true`.



<a name="ODE-one-well"></a>
## Fitting ODE function for one well
```
 fitting_one_well_ODE_constrained(data::Matrix{Float64},
    name_well::String, 
    label_exp::String,
    model::String,
    lb_param::Vector{Float64}, 
    ub_param::Vector{Float64}; 
    param= lb_param .+ (ub_param.-lb_param)./2,
    optmizator =   BBO_adaptive_de_rand_1_bin_radiuslimited(), 
    integrator =KenCarp4(autodiff=true), 
    do_plot=false, 
    path_to_plot="NA", 
    pt_avg=1, 
    pt_smooth_derivative=7,
    smoothing=false,
    type_of_loss="RE",
    blank_array=zeros(100),
    multiple_scattering_correction=false, 
    calibration_OD_curve="NA"  
    )
```
This function performs constrained parameter fitting on a single well's dataset using an ordinary differential equation (ODE) model. It estimates the model parameters within specified lower and upper bounds.

Arguments:

- `data::Matrix{Float64}`: Dataset 
-  `model::String`: ODE model to use
- `name_well::String`: Name of the well.
- `label_exp::String`: Label of the experiment.
- `lb_param::Vector{Float64}`: Lower bounds of the model parameters.
- `ub_param::Vector{Float64}`: Upper bounds of the model parameters.


 Key Arguments:

- `param= lb_param .+ (ub_param.-lb_param)./2`: Initial guess for the model parameters.
- `integrator =KenCarp4(autodiff=true)' sciML integrator
- `optmizator =   BBO_adaptive_de_rand_1_bin_radiuslimited()` optimizer from optimizationBBO
- `do_plot=true`: Whether to generate plots or not.
- `path_to_plot="NA"`: Path to save the plots.
- `pt_avg=7`: Number of points to generate the initial condition or do the rolling avg smoothing.
- `smoothing=false`: Whether to apply smoothing to the data or not.
- `type_of_loss:="RE" `: Type of loss function to be used. (options= "RE", "L2", "L2_derivative" and "blank_weighted_L2")
- `blank_array=zeros(100)`: Data of all blanks in single array.
- `verbose=false`: Whether to enable verbose output.
- `write_res=true`: Whether to write results.
- `multiple_scattering_correction=false`: Whether or not correct the data qith a calibration curve.
- `calibration_OD_curve="NA"`: The path where the .csv calibration data are located, used only if `multiple_scattering_correction=true`.


<a name="ODE-file"></a>
## Fitting ODE function for one file
```
   fit_file_ODE(
    label_exp::String,
    path_to_data::String,
    path_to_annotation::String,
    model::String, 
    lb_param::Vector{Float64},
    ub_param::Vector{Float64}; 
    optmizator =   BBO_adaptive_de_rand_1_bin_radiuslimited(),
    integrator = KenCarp4(autodiff=true), 
    path_to_results="NA",
    path_to_plot="NA", 
    loss_type="RE", 
    smoothing=false, 
    do_plot=false,
    verbose=false, 
    write_res=false, 
    pt_avg=1, 
    pt_smooth_derivative=7,
    blank_subtraction="avg_blank", 
    fit_replicate=false, 
    correct_negative="thr_correction",
    thr_negative=0.01, 
    multiple_scattering_correction=false, 
    calibration_OD_curve="NA"  
    )
```
This function fits an ordinary differential equation (ODE) model to a single file's data. 

Arguments:


- `path_to_data::String`: path to the csv file of data
- `path_to_annotation::String` path to the annotation of the dataset
-  `model::String`: ODE model to use
- `label_exp::String`: Label of the experiment.
- `lb_param::Vector{Float64}`: Lower bounds of the model parameters.
- `ub_param::Vector{Float64}`: Upper bounds of the model parameters.


 Key Arguments:

- `param= lb_param .+ (ub_param.-lb_param)./2`: Initial guess for the model parameters.
- `integrator =KenCarp4(autodiff=true)' sciML integrator
- `optmizator =   BBO_adaptive_de_rand_1_bin_radiuslimited()` optimizer from optimizationBBO
- `do_plot=true`: Whether to generate plots or not.
- `path_to_plot="NA"`: Path to save the plots.
- `pt_avg=7`: Number of points to generate the initial condition or do the rolling avg smoothing.
- `smoothing=false`: Whether to apply smoothing to the data or not.
- `type_of_loss:="RE" `: Type of loss function to be used. (options= "RE", "L2", "L2_derivative" and "blank_weighted_L2")
- `blank_array=zeros(100)`: Data of all blanks in single array.
- `verbose=false`: Whether to enable verbose output.
- `write_res=true`: Whether to write results.
- `multiple_scattering_correction=false`: Whether or not correct the data qith a calibration curve.
- `calibration_OD_curve="NA"`: The path where the .csv calibration data are located, used only if `multiple_scattering_correction=true`.
- fit_replicate=false,  if true the average between replicates is fitted.
  
<a name="custom-ODE"></a>
## Fitting custom ODE function for one well
```
fitting_one_well_custom_ODE(data::Matrix{Float64},
    name_well::String, 
    label_exp::String,
    model::Any, 
    lb_param::Vector{Float64}, 
    ub_param::Vector{Float64},
    n_equation::Int; 
    param= lb_param .+ (ub_param.-lb_param)./2,
    optmizator =   BBO_adaptive_de_rand_1_bin_radiuslimited(), 
    integrator =KenCarp4(autodiff=true),
    do_plot=false, 
    path_to_plot="NA", 
    pt_avg=1, 
    pt_smooth_derivative=7,
    smoothing=false, 
    type_of_loss="RE", 
    blank_array=zeros(100), 
    multiple_scattering_correction=false,
    calibration_OD_curve="NA"  
    )
```

This function is designed for fitting an ordinary differential equation (ODE) model to a dataset representing the growth curve of a microorganism in a well. It utilizes a customizable ODE model, optimization methods, and integration techniques for parameter estimation.

 Arguments:

- `data::Matrix{Float64}`: The dataset with the growth curve, where the first row represents times, and the second row represents optical density (OD).
- `name_well::String`: The name of the well.
- `label_exp::String`: The label of the experiment.
- `model::Any`: The ODE model to use.
- `lb_param::Vector{Float64}`: Lower bounds for the parameters.
- `ub_param::Vector{Float64}`: Upper bounds for the parameters.
- `n_equation::Int`: The number of ODEs in the system.
  
Key   Arguments:

- `param= lb_param .+ (ub_param.-lb_param)./2`: Initial guess for the parameters.
- `optmizator=BBO_adaptive_de_rand_1_bin_radiuslimited()`: The optimization method to use.
- `integrator=KenCarp4(autodiff=true)`: The integrator for solving the ODE.
- `do_plot=false`: Whether to generate plots or not.
- `path_to_plot="NA"`: Path to save the generated plots.
- `pt_avg=1`: Number of points to generate the initial condition or do the rolling avg smoothing.
- `pt_smooth_derivative=7`: Number of points for evaluation of specific growth rate. If <2 it uses interpolation algorithm otherwise a sliding window approach.
- `smoothing=false`: Determines whether smoothing is applied to the data.
- `type_of_loss="RE"`: Type of loss used for optimization (options= "RE", "L2", "L2_derivative" and "blank_weighted_L2")
- `blank_array=zeros(100)`: Data representing blanks for correction.
- `multiple_scattering_correction=false`: If `true`, uses a given calibration curve to correct the data.
- `calibration_OD_curve="NA"`: The path to the calibration curve used for data correction.


<a name="Sensitivity-analysis"></a>
## Sensitivity analysis
```
 one_well_morris_sensitivity(data::Matrix{Float64}, 
    name_well::String,
    label_exp::String, 
    model::String, 
    lb_param::Vector{Float64}, 
    ub_param::Vector{Float64}; 
    N_step_morris =7,
    optmizator =   BBO_adaptive_de_rand_1_bin_radiuslimited(), 
    integrator =KenCarp4(autodiff=true), 
    pt_avg=1, 
    pt_smooth_derivative=7,
    write_res=false,
    smoothing=false,
    type_of_loss="RE", 
    blank_array=zeros(100),
    multiple_scattering_correction=false, 
    calibration_OD_curve="NA"  
    )
```

This function is designed to perform Morris sensitivity analysis on a dataset representing the growth curve of a microorganism in a well. It assesses the sensitivity of the model to variations in input parameters of the optimization.

Arguments:

- `data::Matrix{Float64}`: The dataset with the growth curve, where the first row represents times, and the second row represents optical density (OD).
- `name_well::String`: The name of the well.
- `label_exp::String`: The label of the experiment.
- `model::String`: The ODE model to use.
- `lb_param::Vector{Float64}`: Lower bounds for the parameters.
- `ub_param::Vector{Float64}`: Upper bounds for the parameters.


Key Arguments:

- `N_step_morris=7`: Number of steps for the Morris sensitivity analysis.
- `optmizator=BBO_adaptive_de_rand_1_bin_radiuslimited()`: The optimization method to use.
- `integrator=KenCarp4(autodiff=true)`: The integrator for solving the ODE.
- `pt_avg=1`: Number of points to generate the initial condition or do the rolling avg smoothing.
- `pt_smooth_derivative=7`: Number of points for evaluation of specific growth rate. If <2 it uses interpolation algorithm otherwise a sliding window approach.
- `write_res=false`: If `true`, writes the sensitivity analysis results to a file.
- `smoothing=false`: Determines whether smoothing is applied to the data.
- `type_of_loss="RE"`: Type of loss used for optimization (options= "RE", "L2", "L2_derivative" and "blank_weighted_L2")
- `blank_array=zeros(100)`: Data representing blanks for correction.
- `multiple_scattering_correction=false`: If `true`, uses a given calibration curve to correct the data.
- `calibration_OD_curve="NA"`: The path to the calibration curve used for data correction.


<a name="Model-selection"></a>

## Model selection
```
ODE_Model_selection(data::Matrix{Float64}, 
    name_well::String, 
    label_exp::String,
    models_list::Vector{String}, 
    lb_param_array::Any,
    ub_param_array::Any; 
    optmizator =   BBO_adaptive_de_rand_1_bin_radiuslimited(), 
    integrator = KenCarp4(autodiff=true), 
    pt_avg = 1 , 
    beta_penality = 2.0, 
    smoothing= false, 
    type_of_loss="L2",
    blank_array=zeros(100), 
    plot_best_model=false, 
    path_to_plot="NA",
    pt_smooth_derivative=7,
    multiple_scattering_correction=false, 
    calibration_OD_curve="NA", 
    verbose=false
)
```

This function performs model selection based on a dataset representing the growth curve of a microorganism in a well. It evaluates multiple ODE models and selects the best-fitting model using the Akaike Information Criterion (AIC).

Arguments:

- `data::Matrix{Float64}`: The dataset with the growth curve, where the first row represents times, and the second row represents optical density (OD).
- `name_well::String`: The name of the well.
- `label_exp::String`: The label of the experiment.
- `models_list::Vector{String}`: A vector of ODE models to evaluate.
- `lb_param_array::Any`: Lower bounds for the parameters (compatible with the models).
- `ub_param_array::Any`: Upper bounds for the parameters (compatible with the models).

Key Arguments:

- `optmizator=BBO_adaptive_de_rand_1_bin_radiuslimited()`: The optimization method to use.
- `integrator=KenCarp4(autodiff=true)`: The integrator for solving the ODE.
- `pt_avg=1`: Number of points to generate the initial condition or do the rolling avg smoothing.
- `beta_penality=2.0`: Penalty for AIC evaluation.
- `smoothing=false`: Determines whether smoothing is applied to the data.
- `type_of_loss="L2"`: Type of loss used for optimization (options= "RE", "L2", "L2_derivative" and "blank_weighted_L2")
- `blank_array=zeros(100)`: Data representing blanks for correction.
- `plot_best_model=false`: If `true`, the results of the best-fit model will be plotted.
- `path_to_plot="NA"`: Path to save the generated plots.
- `pt_smooth_derivative=7`: Number of points for evaluation of specific growth rate. If <2 it uses interpolation algorithm otherwise a sliding window approach.
- `multiple_scattering_correction=false`: If `true`, uses a given calibration curve to correct the data.
- `calibration_OD_curve="NA"`: The path to the calibration curve used for data correction.
- `verbose=false`: If `true`, enables verbose output.

<a name="cdp"></a>
## Change point detection
```
cpd_local_detection(data::Matrix{Float64},
    n_max_cp::Int;
    type_of_detection="lsdd",
    type_of_curve="original", 
    pt_derivative = 0,
    size_win =2)

```
This function performs change point detection on a dataset, identifying local changes in the growth curve. It uses various algorithms based on user-defined parameters.

Arguments:

- `data::Matrix{Float64}`: The dataset with the growth curve, where the first row represents times, and the second row represents optical density (OD).
- `n_max_cp::Int`: The maximum number of change points to detect.

Key Arguments:

- `type_of_detection="lsdd"`: Type of change point detection algorithm. Options are "lsdd" or piecewise linear fitting 
- `type_of_curve="deriv"`: Type of curve used for the change point detection. Options are "deriv" for the  derivative/specific gr or "original" for growth curve.
- `pt_derivative=0`: Number of points to evaluate the derivative or specific growth rate. If 0, numerical derivative is used; if >1, specific growth rate is calculated with the given window size.
- `size_win=2`: Size of the sliding window used in all detection methods.
  
<a name="cdp-fixed"></a>
## Fitting segmented ODE with fixed change-point number
```
selection_ODE_fixed_change_points(data_testing::Matrix{Float64}, 
    name_well::String,
    label_exp::String,
    list_of_models::Vector{String}, 
    list_lb_param::Any,
    list_ub_param::Any, 
    n_change_points::Int;
    type_of_loss="L2", 
    optmizator =   BBO_adaptive_de_rand_1_bin_radiuslimited(), 
    integrator = KenCarp4(autodiff=true), 
    type_of_detection =  "lsdd",
    type_of_curve = "original", 
    smoothing=false,
    pt_avg=1,
    do_plot=false, 
    path_to_plot="NA", 
    win_size=2, 
    pt_smooth_derivative=0,
    multiple_scattering_correction=false, 
    calibration_OD_curve="NA",
    beta_smoothing_ms = 2.0 
    )
```

This function performs model selection for ordinary differential equation (ODE) models while considering fixed change points in a growth curve dataset. It allows for the evaluation of multiple ODE models and considers a specified number of change points.

Arguments:

- `data_testing::Matrix{Float64}`: The dataset with the growth curve, where the first row represents times, and the second row represents optical density (OD).
- `name_well::String`: The name of the well.
- `label_exp::String`: The label of the experiment.
- `list_of_models::Vector{String}`: A vector of ODE models to evaluate.
- `list_lb_param::Any`: Lower bounds for the parameters (compatible with the models).
- `list_ub_param::Any`: Upper bounds for the parameters (compatible with the models).

 Key Arguments:

- `n_change_points::Int`: The number of fixed change points to consider.
- `type_of_loss="L2"`: Type of loss used for optimization(options= "RE", "L2", "L2_derivative" and "blank_weighted_L2").
- `optmizator=BBO_adaptive_de_rand_1_bin_radiuslimited()`: The optimization method to use.
- `integrator=KenCarp4(autodiff=true)`: The integrator for solving the ODE.
- `type_of_detection="lsdd"`: Type of change point detection algorithm. Options are "lsdd" or piecewise linear fitting 
- `type_of_curve="original"`: Type of curve used for the change point detection. Options are "deriv" for the  derivative/specific gr or "original" for growth curve.
- `smoothing=false`: Determines whether smoothing is applied to the data.
- `pt_avg=1`: Number of points to generate the initial condition or do the rolling avg smoothing.
- `do_plot=false`: Whether to generate plots or not.
- `path_to_plot="NA"`: Path to save the generated plots.
- `win_size=2`: Number of points for the  window of the change point detection algorithm.
- `pt_smooth_derivative=0`: Number of points for evaluation of specific growth rate. If <2 it uses interpolation algorithm otherwise a sliding window approach.
- `multiple_scattering_correction=false`: If `true`, uses a given calibration curve to correct the data.
- `calibration_OD_curve="NA"`: The path to the calibration curve used for data correction.
- `beta_smoothing_ms=2.0`: Penality parameter of the Akaike Information Criterion (AIC) penalty.

<a name="cdp-search"></a>
## Fitting segmented ODE with direct search for a maximum number of change-points 
```
ODE_selection_NMAX_change_points(data_testing::Matrix{Float64}, 
    name_well::String, 
    label_exp::String, 
    list_lb_param::Any, 
    list_ub_param::Any, 
    list_of_models::Vector{String}, 
    n_max_change_points::Int; 
    optmizator =   BBO_adaptive_de_rand_1_bin_radiuslimited(),  
    integrator = KenCarp4(autodiff=true),
    type_of_loss="L2", # type of used loss 
    type_of_detection =  "lsdd",
    type_of_curve = "deriv", 
    pt_avg = pt_avg , 
    smoothing= true, 
    do_plot=false, 
    path_to_plot="NA", 
    path_to_results="NA",
    win_size=2, 
    pt_smooth_derivative=7,
    penality_parameter=2.0,
    multiple_scattering_correction="false", 
    calibration_OD_curve="NA",  
   save_all_model=false )
```
This function fits segmented ordinary differential equation (ODE) models to a growth curve dataset using direct search for a maximum number of change-points. It allows for the evaluation of multiple ODE models with a varying number of change-points.

Arguments:

- `data_testing::Matrix{Float64}`: The dataset with the growth curve, where the first row represents times, and the second row represents optical density (OD) or fluorescence.
- `name_well::String`: The name of the well.
- `label_exp::String`: The label of the experiment.
- `list_lb_param::Any`: Lower bounds for the parameters (compatible with the models).
- `list_ub_param::Any`: Upper bounds for the parameters (compatible with the models).
- `list_of_models::Vector{String}`: A vector of ODE models to evaluate.
- `n_max_change_points::Int`: The maximum number of change-points to consider.

Key Arguments:

  
- `optmizator=BBO_adaptive_de_rand_1_bin_radiuslimited()`: The optimization method to use.
- `integrator=KenCarp4(autodiff=true)`: The integrator for solving the ODE.
- `type_of_loss="L2"`: Type of loss used for optimization (options: "L2" for squared loss).
- `type_of_detection="lsdd"`: Type of change point detection algorithm. Options are "lsdd" for piecewise linear fitting on the specific growth rate.
- `type_of_curve="original"`: Type of curve used for the change point detection. Options are "deriv" for the  derivative/specific gr or "original" for growth curve.
- `pt_avg=pt_avg`: Number of points to generate the initial condition or do the rolling avg smoothing.
- `smoothing=true`: Determines whether smoothing is applied to the data.
- `do_plot=false`: Whether to generate plots or not.
- `path_to_plot="NA"`: Path to save the generated plots.
- `path_to_results="NA"`: Path to save the fitting results.
- `win_size=2`: Number of points for the  window of the change point detection algorithm.
- `pt_smooth_derivative=0`: Number of points for evaluation of specific growth rate. If <2 it uses interpolation algorithm otherwise a sliding window approach.
- `penality_parameter=2.0`: Parameter for penalizing the change in the number of parameters.
- `multiple_scattering_correction=false`: If `true`, uses a given calibration curve to correct the data.
- `calibration_OD_curve="NA"`: The path to the calibration curve used for data correction.
- `save_all_model=false`: If `true`, saves fitting results for all evaluated models.


<a name="models"></a>
# The mathematical models
## ODEs for bacterial growth

TO DO

The models and their parameters are sumarrized in the following table
| Model                                  |  Parameters                                       |
| --------------------------------------- | -------------------------------------------------------- |
| Diauxic_piecewise_damped_logistic      | gr_1, N_max, shape_1, lag, linear_const, t_shift, gr_2, N_max_2, shape_2, end_second_lag, lag_2_gr |
| Diauxic_replicator_1                    | gr, N_max, lag, arbitrary_const, linear_const |
| Diauxic_replicator_2                    | gr, N_max, lag, arbitrary_const, linear_const, growth_stationary |
| HPM                                    | gr, exit_lag_rate, N_max |
| HPM_3_death                            | gr, exit_lag_rate, inactivation_rate, death_rate |
| HPM_3_death_resistance                 | gr, exit_lag_rate, inactivation_rate, death_rate, n_res, n_max |
| HPM_3_inhibition                       | gr, exit_lag_rate, inactivation_rate |
| HPM_inhibition                         | gr, inhibition_rate, gr_inhibition, N_max |
| HPM_exp                                | gr, exit_lag_rate |
| ODEs_HPM_SR                            | gr, gr_phage, scale, death_rate, resistance_rate |
| baranyi_exp                            | gr, lag_time, shape |
| baranyi_richards                       | gr, N_max, lag_time, shape |
| baranyi_roberts                        | gr, N_max, lag_time, shape_1, shape_2 |
| bertalanffy_richards                   | gr, N_max, shape |
| exponential                             | gr |
| four_piecewise                          | gr, gr_2, gr_3, gr_4, lag, t_decay_gr, t_stationary |
| gbsm_piecewise                         | gr, a_1, b_1, c, a_2, b_2 |
| gompertz                               | gr, N_max |
| hyper_gompertz                         | gr, N_max, shape |
| hyper_logistic                         | doubling_time, gr, N_max, shape |
| huang                                  | gr, N_max, lag |
| logistic                               | gr, N_max |
| logistic                               | gr, N_max, shape |
| ode_von_bertalanffy                    | alpha, beta, a, b |
| piecewise_damped_logistic              | gr, N_max, lag, shape, linear_const |
| triple_piecewise                       | gr, gr_2, gr_3, lag, t_stationary |
| triple_piecewise_bertalanffy_richards  | gr, gr_lag, t_lag, t_stationary, gr_stat, shape, N_max |
| triple_piecewise_damped_logistic       | gr, gr_2, gr_3, lag, t_stationary, N_max |
| triple_piecewise_sublinear             | gr, gr_2, gr_3, lag, t_stationary, N_max |


## Stochastic models for bacterial growth


In the stochastic version of the growth models, the growth rate of each population component (denoted as $\mu_i$) is evaluated based on the concentration of the limiting nutrient. The user is required to specify the starting amount of nutrients and the volume of the solution. Various kinetic growth models are considered.

### Monod Model

The Monod model is described by the following equation:

\[
\mu([\text{Nut.}]; k_1, \mu_\text{max}) = \mu_\text{max} \frac{[\text{Nut.}]}{k_1 + [\text{Nut.}]}.
\]

### Haldane Model

The Haldane model is expressed as:

\[
\mu([\text{Nut.}]; k_1, k_2, \mu_\text{max}) = \mu_\text{max} \frac{[\text{Nut.}]}{k_1 + [\text{Nut.}] + \frac{k_2}{[\text{Nut.}]^2}}.
\]

### Blackman Model

The Blackman model is given by:

\[
\mu([\text{Nut.}]; k_1, \mu_\text{max}) = \mu_\text{max} \frac{[\text{Nut.}]}{k_1 + [\text{Nut.}]}.
\]

### Tessier Model

The Tessier model is represented as:

\[
\mu([\text{Nut.}]; k_1, \mu_\text{max}) = \mu_\text{max} (1 - e^{k_1[\text{Nut.}] }).
\]

### Moser Model

The Moser model is defined by:

\[
\mu([\text{Nut.}]; k_1, k_2, \mu_\text{max}) = \mu_\text{max} \frac{[\text{Nut.}]^{k_2}}{k_1 + [\text{Nut.}]^{k_2}}.
\]

### Aiba-Edwards Model

The Aiba-Edwards model is given by:

\[
\mu([\text{Nut.}]; k_1, k_2, \mu_\text{max}) = \mu_\text{max} \frac{[\text{Nut.}]}{k_1 + [\text{Nut.}]} e^{-\frac{[\text{Nut.}]}{k_2}}.
\]

### Verhulst Model

The Verhulst model is defined as:

\[
\mu(N; N_\text{max}, \mu_\text{max}) = \mu_\text{max} \left(1 - \frac{N}{N_\text{max}}\right).
\]

Where $[\text{Nut.}]$ is the limiting nutrient concentration, $\mu_\text{max}$ is the maximum possible growth rate, $k_1$ and $k_2$ are numerical constants with meanings depending

## Type of loss functions

`type_of_loss = "L2"`: Minimize the L2 norm of the difference between the numerical solution of an ODE and the given data.

`type_of_loss = "L2_derivative"`: Minimize the L2 norm of the difference between the derivatives of the numerical solution of an ODE and the corresponding derivatives of the data.

`type_of_loss ="RE" `: Minimize the relative error between the solution and data 

`type_of_loss = "blank_weighted_L2"` : Minimize a weighted version of the L2 norm, where the difference between the solution and data is weighted based on a distribution obtained from empirical blank data. 

<a name="examples"></a>
# Examples and Tutorial


## Table of Contents
1. [Simulating Data with J-MAKi](#simulating-data)
2. [Data Preprocessing](#data-preprocessing)
3. [Smoothing and Multiple Scattering Correction](#smoothing-correction)
4. [Model Fitting](#model-fitting)
    - [Fitting ODE Models](#fitting-ode)
    - [Custom ODE Fitting](#custom-ode-fitting)
    - [Sensitivity Analysis](#sensitivity-analysis)
    - [Model Selection](#model-selection)
5. [Change Point Detection](#change-point-detection)
6. [Conclusion](#conclusion)


##  Simulating Data with J-MAKi

To demonstrate the package, we'll begin by simulating data using a specific ODE model. The following code snippet simulates a growth curve using a triple-piecewise damped logistic model:

```julia
<a name="simulating-data"></a>

# Simulating data with a triple-piecewise damped logistic model
model = "triple_piecewise_damped_logistic"
n_start = [0.1]
tstart = 0.0
tmax = 600.0
delta_t = 10.0
integrator = KenCarp4()
param_of_ode = [0.06, 1.0, 200, 0.5, 0.001, 450, -0.0002]

sim = ODE_sim(model, n_start, tstart, tmax, delta_t, integrator, param_of_ode)

# Plotting scatterplot of simulated data
Plots.scatter(sim, xlabel="Time", ylabel="Arb. Units", label=["Simulated Data" nothing], color=:blue, size=(300, 300))
```

In the above code, we use the `ODE_sim` function to simulate data based on a specified ODE model and parameters. The resulting data is then plotted.

<a name="data-preprocessing"></a>
##  Data Preprocessing

Next, we'll explore data preprocessing steps, including adding noise, smoothing, and multiple scattering correction:

```julia
# ... (Previous code)

# Adding uniform random noise to the simulated data
noise_uniform = rand(Uniform(-0.05, 0.05), length(sim.t))
data_t = reduce(hcat, sim.t)
data_o = reduce(hcat, sim.u)
data_OD = vcat(data_t, data_o)
data_OD[2, :] = data_OD[2, :] .+ noise_uniform

# Plotting scatterplot of data with noise
Plots.scatter(data_OD[1, :], data_OD[2, :], xlabel="Time", ylabel="Arb. Units", label=["Noisy Data" nothing], color=:blue, markersize=2, size=(300, 300))

# Smoothing of the data with rolling average
data_OD_smooth = smoothing_data(data_OD, 7)
data_OD_smooth = Matrix(data_OD_smooth)

# Plotting scatterplot of smoothed data
Plots.scatter(data_OD_smooth[1, :], data_OD_smooth[2, :], xlabel="Time", ylabel="Arb. Units", label=["Smoothed Data" nothing], markersize=2, color=:blue, size=(300, 300))

# Multiple scattering correction (requires an external file)
# Comment out the following line if correction is not needed
data_OD_smooth2 = correction_OD_multiple_scattering(data_OD_smooth, "/path/to/calibration/file.csv")

# Plotting scatterplot of pre-processed data
Plots.scatter(data_OD_smooth2[1, :], data_OD_smooth2[2, :], xlabel="Time", ylabel="Arb. Units", label=["Pre-processed Data" nothing], markersize=2, color=:blue, size=(300, 300))
```

In this section, we add noise to the simulated data, smooth it using a rolling average, and perform multiple scattering correction if applicable.

<a name="model-fitting"></a>
## Model Fitting

Now, let's explore model fitting using J-MAKi. We'll cover fitting ODE models, custom ODEs, sensitivity analysis, and model selection.

<a name="fitting-ode"></a>
### Fitting ODE Models

We start by fitting the simulated data with a predefined ODE model:

```julia
# ... (Previous code)

# Fitting the simulated data with a log-linear model
path_to_plotting = "/path/to/save/plots/"
res_log_lin = fitting_one_well_Log_Lin(data_OD_smooth2, "test", "test", do_plot=true, path_to_plot=path_to_plotting)
```

In this snippet, the `fitting_one_well_Log_Lin` function is used to fit the data with a log-linear model.

<a name="custom-ode-fitting"></a>
### Custom ODE Fitting

We can also fit data using a custom ODE. Here, we demonstrate fitting with a custom two-equation ODE:

```julia
# ... (Previous code)

# Defining the custom ODE function
function ODE_custom(du, u, param, t)
    du[1] = u[1] * (1 - u[1]) * param[2] + param[1] * u[1]
    du[2] = +u[1] * param[2] + param[4] * u[2] * (1 - (u[1] + u[2]) / param[3])
end

custom_ub = [1.2, 1.1, 2.0, 20]
custom_lb = [0.0001, 0.00000001, 0.00, 0]

# Fitting data with the custom ODE
results_custom = fitting_one_well_custom_ODE(data_OD_smooth2, "test", "test_model_custom", ODE_custom, custom_lb, custom_ub, 2, do_plot=true, path_to_plot=path_to_plotting)
```

In this example, `ODE_custom` represents a custom ODE function with two equations. The `fitting_one_well_custom_ODE` function is used for fitting.

<a name="sensitivity-analysis"></a>
### Sensitivity Analysis

We can perform sensitivity analysis to evaluate the impact of parameter variations:

```julia
# ... (Previous code)

# Performing sensitivity analysis using Morris method
n_step_sensitivity = 3
sensitivity_test = one_well_morris_sensitivity(data_OD_smooth2
