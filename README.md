# JMAKi

JMAKi, a versatile software tool that utilizes Ordinary Differential Equations (ODEs) to fit bacterial growth data from plate reader experiments. 
With JMAKi it is possible to simulate, fit, perform model selection, and conduct sensitivity analysis for multi-well plate reader experiments.
The parameter fitting in JMAKi is defined as a constrained optimization problem, which is solved using a differential evolution non-linear optimizer.


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

# The main functions of JMAKi


## Data and annotation formatting
JMAKi can operate directly on data files or inside the julia notebook.

In particular two API call direclty the files. 

In these cases the user must input  the paths to  a .csv data file and a .csv annotation

Instead for all other API the standard file formatting is a 2xn_times Matrix of Float64 .,e.g.:

```
 0.0        2.0       4.0       6.0       8.0        10.0       10.0       12.0       14.0       16.0       18.0       20.0       22.0      24.0      26.0       28.0       30.0       32.0       34.0       36.0       …  
 0.0912154  0.107956  0.105468  0.101727  0.0931484   0.106318   0.103697   0.139821   0.173598   0.204888   0.251052   0.289018   0.31298   0.33752   0.359356   0.370861   0.376347   0.383732   0.398496   0.384511 …  

```
The first row should be time and the second the quantity to be fitteted (e.g., Optical Density or CFU)



## Simulate ODE
```

    ODE_sim(model::String, n_start::Vector{Float64}, tstart::Float64, tmax::Float64,
        delta_t::Float64, integrator::Any, param_of_ode::Vector{Float64})
```
This function performs an ODE simulation of a model, considering the initial conditions, time range, and integration parameters.

- `model::String`: The model to simulate.
- `n_start::Vector{Float64}`: The starting conditions.
- `tstart::Float64`: The start time of the simulation.
- `tmax::Float64`: The final time of the simulation.
- `delta_t::Float64`: The time step of the output.
- `integrator::Any`: The chosen solver from the SciML ecosystem for ODE integration.
- `param_of_ode::Vector{Float64}`: The parameters of the ODE model.




##  Stochastic simulation
```
    stochastic_sim(model::String, n_start::Int, n_mol_start::Float64, tstart::Float64, tmax::Float64,
        delta_t::Float64, k_1_val::Float64, k_2_val::Float64, alpha_val::Float64, lambda::Float64,
        n_mol_per_birth::Float64, volume::Float64)
```
This function performs a stochastic simulation of a model, considering cell growth and nutrient consumption over time.

- `model::String`: The model to simulate.
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


## Fitting growth rate with log-lin fitting for one well
```
    fitting_one_well_Log_Lin_Fit(data::Matrix{Float64}, name_well::String, label_exp::String,
        do_plot::Bool, path_to_plot::String, pt_avg::Int, pt_smoothing_derivative::Int,
        pt_min_size_of_win::Int, type_of_win::String, threshold_of_exp::Float64, do_error_plot::Bool)
```
This function fits a logarithmic-linear model to a single well's data and performs analysis such as plotting and error calculation.

- `data::Matrix{Float64}`: The dataset of OD/fluorescence values.
- `name_well::String`: The name of the well.
- `label_exp::String`: The label of the experiment.
- `do_plot::Bool`: Whether to generate and save plots.
- `path_to_plot::String`: The path to save the plots.
- `pt_avg::Int`: The number of points to generate the initial condition.
- `pt_smoothing_derivative::Int`: The number of points to smooth the derivative.
- `pt_min_size_of_win::Int`: The minimum size of the exponential windows in the number of smoothed points.
- `type_of_win::String`: How the exponential phase window is selected ("maximum" or "global_thr").
- `threshold_of_exp::Float64`: The threshold of the growth rate in quantile to define the exponential windows.
- `do_error_plot::Bool`: Whether to generate a distribution plot of the goodness of fit in the interval of the fitted growth rate.

## Fitting growth rate with log-lin fitting for one file
```
    fit_one_file_Log_Lin(label_exp::String, temp_file::String, path_to_data::String,
        path_to_annotation::String, path_to_results::String, path_to_plot::String,
        do_plot::Bool, verbose::Bool, write_res::Bool, pt_avg::Int, pt_smoothing_derivative::Int,
        pt_min_size_of_win::Int, type_of_win::String, threshold_of_exp::Float64,
        blank_subtraction::String, fit_replicate::Bool, correct_negative::String,
        thr_negative::Float64, do_error_plot::Bool)
```
This function fits a logarithmic-linear model to a single file's data. It performs model fitting, error analysis, and provides various options for customization.

- `label_exp::String`: Label of the experiment.
- `temp_file::String`: Name of the file to analyze.
- `path_to_data::String`: Path to the folder containing the data.
- `path_to_annotation::String`: Path to the annotation of the wells.
- `path_to_results::String`: Path to save the results.
- `path_to_plot::String`: Path to save the plots.
- `do_plot::Bool`: Whether to generate and visualize plots of the data.
- `verbose::Bool`: Whether to enable verbose output.
- `write_res::Bool`: Whether to write results.
- `pt_avg::Int`: Number of points to use for smoothing average.
- `pt_smoothing_derivative::Int`: Number of points to smooth the derivative.
- `pt_min_size_of_win::Int`: Minimum size of the exponential windows in number of smoothed points.
- `type_of_win::String`: How the exponential phase window is selected ("maximum" or "global_thr").
- `threshold_of_exp::Float64`: Threshold of growth rate in quantile to define the exponential windows.
- `blank_subtraction::String`: How to use blank data for subtraction (options: "NO", "avg_subtraction", "time_avg").
- `fit_replicate::Bool`: If `true`, fit the average between replicates; if `false`, fit all replicates independently.
- `correct_negative::String`: Method to correct negative values (options: "thr_correction", "blank_correction").
- `thr_negative::Float64`: Threshold value used only if `correct_negative == "thr_correction"`.
- `do_error_plot::Bool`: Generate a distribution plot of goodness of fit in the interval of fitted growth rate.




## Fitting ODE function for one well
```
    fitting_one_well_constrained(data::Matrix{Float64}, name_well::String, label_exp::String,
        lb_param::Vector{Float64}, ub_param::Vector{Float64}, param::Vector{Float64},
        model::String, do_plot::Bool, path_to_plot::String, pt_avg::Int, smoothing::Bool,
        type_of_loss::String, blank_array::Vector{Float64}, verbose::Bool,
        error_analysis::Bool, write_res::Bool)
```
This function performs constrained parameter fitting on a single well's dataset using an ordinary differential equation (ODE) model. It estimates the model parameters within specified lower and upper bounds.

- `data::Matrix{Float64}`: Dataset of x times y OD/fluorescence.
- `name_well::String`: Name of the well.
- `label_exp::String`: Label of the experiment.
- `lb_param::Vector{Float64}`: Lower bounds of the model parameters.
- `ub_param::Vector{Float64}`: Upper bounds of the model parameters.
- `param::Vector{Float64}`: Initial guess for the model parameters.
- `model::String`: ODE model to use.
- `do_plot::Bool`: Whether to generate plots or not.
- `path_to_plot::String`: Path to save the plots.
- `pt_avg::Int`: Number of points to generate the initial condition.
- `smoothing::Bool`: Whether to apply smoothing to the data or not.
- `type_of_loss::String`: Type of loss function to be used.
- `blank_array::Vector{Float64}`: Data of all blanks.
- `verbose::Bool`: Whether to enable verbose output.
- `error_analysis::Bool`: Perform error analysis.
- `write_res::Bool`: Whether to write results.



## Fitting ODE function for one file
```
    fit_one_file_ODE(label_exp::String, temp_file::String, path_to_data::String,
        path_to_annotation::String, path_to_results::String, path_to_plot::String,
        model::String, loss_type::String, smoothing::Bool, do_plot::Bool, verbose::Bool,
        write_res::Bool, pt_avg::Int, lb_param::Vector{Float64}, ub_param::Vector{Float64},
        blank_subtraction::String, fit_replicate::Bool, error_analysis::Bool,
        correct_negative::String, thr_negative::Float64)
```
This function fits an ordinary differential equation (ODE) model to a single file's data. It performs model fitting, error analysis, and provides various options for customization.

- `label_exp::String`: Label of the experiment.
- `temp_file::String`: Name of the file to analyze.
- `path_to_data::String`: Path to the folder containing the data.
- `path_to_annotation::String`: Path to the annotation of the wells.
- `path_to_results::String`: Path to save the results.
- `path_to_plot::String`: Path to save the plots.
- `model::String`: string of the model.
- `loss_type::String`: Type of the used loss function.
- `smoothing::Bool`: Whether to perform data smoothing with rolling average or not.
- `do_plot::Bool`: Whether to generate and visualize plots of the data.
- `verbose::Bool`: Whether to enable verbose output.
- `write_res::Bool`: Whether to write results.
- `pt_avg::Int`: Number of points to use for smoothing average.
- `lb_param::Vector{Float64}`: Array of the lower bounds of the parameters.
- `ub_param::Vector{Float64}`: Array of the upper bounds of the parameters.
- `blank_subtraction::String`: How to use blank data for subtraction (options: "NO", "avg_subtraction", "time_avg").
- `fit_replicate::Bool`: If `true`, fit the average between replicates; if `false`, fit all replicates independently.
- `error_analysis::Bool`: Perform error analysis.
- `correct_negative::String`: Method to correct negative values (options: "thr_correction", "blank_correction").
- `thr_negative::Float64`: Threshold value used only if `correct_negative == "thr_correction"`.



## Fitting PINN fucntion

```
```
## Sensitivity analysis
```

    Sensitivity_ODE_Morris(data::Matrix{Float64}, name_well::String, label_exp::String,
        lb_param::Vector{Float64}, ub_param::Vector{Float64}, model::String, do_plot::Bool,
        path_to_plot::String, pt_avg::Int, smoothing::Bool, type_of_loss::String,
        blank_array::Vector{Float64}, N_step_morris::Int, write_res::Bool)
```
This function performs sensitivity analysis using the Morris method on an ordinary differential equation (ODE) model. It calculates the sensitivity indices of the model parameters with respect to the observed data.

- `data::Matrix{Float64}`: Dataset of x times y OD/fluorescence.
- `name_well::String`: Name of the well.
- `label_exp::String`: Label of the experiment.
- `lb_param::Vector{Float64}`: Lower bounds of the model parameters.
- `ub_param::Vector{Float64}`: Upper bounds of the model parameters.
- `model::String`: ODE model to use.
- `do_plot::Bool`: Whether to generate plots or not.
- `path_to_plot::String`: Path to save the plots.
- `pt_avg::Int`: Number of points to generate the initial condition.
- `smoothing::Bool`: Whether to apply smoothing to the data or not.
- `type_of_loss::String`: Type of loss function to be used.
- `blank_array::Vector{Float64}`: Data of all blanks.
- `N_step_morris::Int`: Number of steps in the Morris method.
- `write_res::Bool`: Results to be written.




# The mathematical models
## ODEs for bacterial growth
The models and their parameters are sumarrized in the following table

| Model              | Parameters                     |
|--------------------|-------------------------------|
| "hyper_gompertz"     | gr, N_max, shape              |
| "hyper_logistic"     |  gr, N_max, shape  TO REDO |
| "gbsm_piecewise     | gr, a_1, b_1, c, a_2, b_2     |
| "bertalanffy_richards | gr, N_max, shape              |
| "logistic"           | gr, N_max                     |
| "gompertz"           | gr, N_max                     |
| "baranyi_richards"   | gr, N_max, lag_time, shape    |
| "baranyi_roberts"    | gr, N_max, lag_time, shape_1, shape_2 |
| "huang"              | gr, N_max, lag                |
| "piecewise_damped_logistic" | gr, N_max, lag, shape, linear_const |
| "triple_piecewise_damped_logistic" | gr, N_max, lag, shape, linear_const, t_stationary, linear_lag |
| "triple_piecewise"   | gr, gr_2, gr_3, lag, t_stationary |
| "four_piecewise"     | gr, gr_2, gr_3, gr_4, lag, t_decay_gr, t_stationary |
| "Diauxic_replicator_1" | gr, N_max, lag, arbitrary_const, linear_const TO REDO|
| "Diauxic_replicator_2" | gr, N_max, lag, arbitrary_const, linear_const, growth_stationary TO REDO |
| "Diauxic_piecewise_damped_logistic" | gr, N_max, lag, shape, linear_const, t_stationary, growth_stationary TO REDO|

## Stochastic models for bacterial growth

Monod,Haldane,Blackman,Tesseir,Moser,Aiba-Edwards,Verhulst

## Type of loss functions

`type_of_loss = "L2"`: Minimize the L2 norm of the difference between the numerical solution of an ODE and the given data.

`type_of_loss = "L2_derivative"`: Minimize the L2 norm of the difference between the derivatives of the numerical solution of an ODE and the corresponding derivatives of the data.

`type_of_loss ="RE" `: Minimize the relative error between the solution and data 

`type_of_loss = "blank_weighted_L2"` : Minimize a weighted version of the L2 norm, where the difference between the solution and data is weighted based on a distribution obtained from empirical blank data. 

##  Numerical integration and optimization options

# Examples and benchmark functions
## Simulation example
### ODE simulation
```
model= "hyper_gompertz";
ic=  [0.01];
tstart = 0.0;
tend = 100.0;
delta_t = 0.5;
param =   [0.02 ,  2  , 1]  ;

sim = ODE_sim( model, 
        ic ,
        tstart,
        tend, 
        delta_t,
        Tsit5(),
        param  
  )
  
  ```
  The output is the same of a SciML outputs
  ```
  plot(sim)
  ```
  It is useful to use stiff integrator for piecewise and the more complicated models
  
  ```
model= "triple_piecewise_damped_logistic";
ic=  [0.02];
tstart = 0.0;
tend = 1000.0;
delta_t = 0.5;
param =   [0.013 ,  1.1  ,100.0 , 1.5, 0.01, 800.0, 00.0001]  ;

sim = ODE_sim( model, 
        ic ,
        tstart,
        tend, 
        delta_t,
        KenCarp4(),
        param  
  )


plot(sim)
```
### Stochastic simulation
Call the function `stochastic_sim`
```
sim = stochastic_sim("Monod", #string of the model
    2000, # number of starting cells
    0.900, # starting concentration of the limiting nutrients
    0.0, # start time of the sim
    3000.0, # final time of the sim
    1.0, # delta t for poisson approx
    0.2,
    21.3, # monod constant
    0.02, # massimum possible growth rate
    342.0, # lag time
    0.001,# nutrient consumed per division (conc)
    100.0
)
```
I the dimension 3 of the output you will find the times 

In the dimension 2 of the output you will find the concentration of the limiting nutrients

In the dimension 1 of the output you will find the number of cells  


```
plot(sim[3],sim[2],xlabel="Time", ylabel="# of cells Arb. Units")
plot(sim[3],sim[1],xlabel="Time", ylabel="[conc of nutrient]Arb. Units")
```

## Fitting one file
### Fitting ODE

Download the folder data
Choose the ODE model 
```
model= "piecewise_damped_logistic"
```
Set the lower bound of parameter given the choosen model
```
ub_piece_wise_log =[ 0.03 , 200000.0 , 1200.0 , 30.0 ,  0.01    ]
lb_piece_wise_log =[ 0.0001 , 100.1, 50.00, 0.0 , -0.01   ]
```
Set the paths of data, results, and plots
```
path_to_data = "your_path_to_main_JMAKi_folder/data/" ;
path_to_plot ="your_path_to_plots/";
path_to_results ="your_path_to_results/";
path_to_annotation ="your_path_to_main_JMAKi_folder/data/annotation_S8R_green.csv"
```
Call the function to fit ODE
```
results_green  = fit_one_file_ODE(  "Green_S8R" , # label of the exp
  "exp_S8R_Green.csv",# name of the file to analyze
  path_to_data, # path to the folder to analyze
  path_to_annotation,# path to the annotation of the wells
  path_to_results, # path where save results
  path_to_plot, # path where to save Plots
  model, # string of the used model to do analysis 
  "RE", # type of the loss, relative error
  false, # do smoothing of data with rolling average
  true, #  do and visulaze the plots of data
  true, # verbose
  true , #  write the results
  7, # number of points to do smoothing average
  lb_piece_wise_log,# array of the array of the lower bound of the parameters
  ub_piece_wise_log ,# array of the array of the upper bound of the parameters
  "avg_blank", # type of blank subtraction
  false, # avg of replicate
  false, # error analysis 
  "blank_correction" , # if "thr_correction" it put a thr on the minimum value of the data with blank subracted, if "blank_correction" uses blank distrib to impute negative values
  100.0 # ignored in this case
) 
```
### Fitting Log Lin

Download the folder data

Set the paths of data, results, and plots
```
path_to_data = "your_path_to_main_JMAKi_folder/data/" ;
path_to_plot ="your_path_to_plots/";
path_to_results ="your_path_to_results/";
path_to_annotation ="your_path_to_main_JMAKi_folder/data/annotation_S8R_green.csv"
```
Call the function to Log Lin model
```

path_to_annotation ="G:/JMAKi.jl-main/data/annotation_S8R_red.csv"

results = fit_one_file_Log_Lin(
    "test" , # label of the exp
    "exp_S8R_Red.csv",# name of the file to analyze
    path_to_data, # path to the folder to analyze
    path_to_annotation,# path to the annotation of the wells
    path_to_results, # path where save results
    path_to_plot, # path where to save Plots
    true, # 1 do and visulaze the plots of data
    true, # 1 true verbose
    true, # write results
    7, # number of points to do smoothing average
    9, # number of poits to smooth the derivative
    4, # minimum size of the exp windows in number of smooted points
    "maximum" ,  # how the exp. phase win is selected, "maximum" of "global_thr"
    0.9,# threshold of growth rate in quantile to define the exp windows
    "avg_blank", # string on how to use blank (NO,avg_blank,time_blank)
    false, # if true the average between replicates is fitted. If false all replicate are fitted indipendelitly
    "thr_correction", # if "thr_correction" it put a thr on the minimum value of the data with blank subracted, if "blank_correction" uses blank distrib to impute negative values
    200.0,  # used only if correct_negative == "thr_correction"
    true #  distribution of goodness of fit in the interval of growth rate fitted
)  

```

## Fitting one Data
### ODE
### Log Lin
## Sensitivity analysis
## PINN Usage
