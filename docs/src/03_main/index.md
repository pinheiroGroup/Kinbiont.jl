# [The main functions of Kimchi](@id API)

1. [Simulate ODE](#simulating-ODE)
2. [Stochastic simulation](#simulating-stochastic)
3. [Plotting a dataset from file](#plot-file)
4. [Specific growth rate evaluation](#Specific_growth_rate_evaluation)
5. [Fitting growth rate with log-lin fitting for one well](#log-lin-one-well)
6. [Fitting growth rate with log-lin fitting for one file](#log-lin-file)
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
```julia

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

Output:

- it returns a standard SciML output (i.e., if `sim =ODE_sim(...)`, then `sim.t` is the array of times and `sim.u` is the array of the simulation)

<a name="simulating-stochastic"></a>
##  Stochastic simulation
```julia
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


Output (if `sim =stochastic_sim(...)`):

- `sim[1]`: array of the number of individuals in the population.
- `sim[2]`: array of the number of biomass equivalent mass of the limiting nutrient concentration.
- `sim[3]`: array of the times of the simulation. 

<a name="plot-file"></a>
## Plotting a dataset from file
```julia
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


Output:

- For this function the output are saved or displayed depending on the values of key arguments.
  
<a name="Specific-growth-rate-evaluation"></a>
## Specific growth rate evaluation
```
specific_gr_evaluation(data_smooted::Matrix{Float64},
     pt_smoothing_derivative::Int)
```
This function evaluate the specific growth rate of a growth curve

Arguments:
- `data_smooted::Matrix{Float64}`: The dataset with the growth curve, where the first row represents times, and the second row represents the variable to fit (e.g., OD), see [data formatting](#data).
-  `pt_smoothing_derivative::Int` Number of points for evaluation of specific growth rate. If <2 it uses interpolation algorithm otherwise a sliding window approach.
<a name="log-lin-one-well"></a>
## Fitting growth rate with log-lin fitting for one well
```julia
 fitting_one_well_Log_Lin(data::Matrix{Float64}, 
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

- `data::Matrix{Float64}`: The dataset with the growth curve, where the first row represents times, and the second row represents the variable to fit (e.g., OD), see [data formatting](#data).
- `name_well::String`: The name of the well.
- `label_exp::String`: The label of the experiment.

Key Arguments:
- `do_plot=true`: Whether to generate and save plots.
- `path_to_plot="NA"`: The path to save the plots, used only if `do_plot=true`.
- `pt_avg=7`: The number of points to do rolling average smoothing.
- `pt_smoothing_derivative=7`:  Number of points for evaluation of specific growth rate. If <2 it uses interpolation algorithm otherwise a sliding window approach.
- `pt_min_size_of_win=7`: The minimum size of the exponential windows in the number of smoothed points.
- `type_of_win="maximum"`: How the exponential phase window is selected ("maximum" or "global_thr").
- `threshold_of_exp=0.9`: The threshold of the growth rate in quantile to define the exponential windows.
- `multiple_scattering_correction=false`: Whether or not correct the data qith a calibration curve.
- `calibration_OD_curve="NA"`: The path where the .csv calibration data are located, used only if `multiple_scattering_correction=true`.


Output:

- an array with the following contents:
`results_lin_log_fit = [label_exp, name_well, start of exp win,  end of exp win,  start of exp win, Maximum specific GR ,specific GR,  2 sigma  CI of GR, doubling time,doubling time - 2 sigma ,doubling time + 2 sigma  , intercept log-lin fitting, 2 sigma intercept ,R^2]`
- The plots of the log-linear fitting and of the dynamics of specific growth rate if `do_plot=true`

<a name="log-lin-file"></a>
## Fitting growth rate with log-lin fitting for one file
```julia
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
This function fits a logarithmic-linear model to a single file's data. 
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
- `pt_smoothing_derivative=7`:  Number of points for evaluation of specific growth rate. If <2 it uses interpolation algorithm otherwise a sliding window approach.
- `pt_min_size_of_win=7`: Minimum size of the exponential windows in number of smoothed points.
- `type_of_win="maximum`: How the exponential phase window is selected ("maximum" or "global_thr").
- `threshold_of_exp=0.9`: Threshold of growth rate in quantile to define the exponential windows.
- `blank_subtraction="avg_blank"`: How to use blank data for subtraction (options: "NO", "avg_subtraction", "time_avg").
- `fit_replicate=false`: If `true`, fit the average between replicates; if `false`, fit all replicates independently.
- `correct_negative="thr_correction`: Method to correct negative values (options: "thr_correction", "blank_correction").
- `thr_negative=0.01`: Threshold value used only if `correct_negative == "thr_correction"`.
- `multiple_scattering_correction=false`: Whether or not correct the data with a calibration curve.
- `calibration_OD_curve="NA"`: The path where the .csv calibration data are located, used only if `multiple_scattering_correction=true`.

Output:

- a matrix wich each column has  the following contents:
`results_lin_log_fit[:,1] = [label_exp, name_well, start of exp win,  end of exp win,  start of exp win, Maximum specific GR ,specific GR,  2 sigma  CI of GR, doubling time,doubling time - 2 sigma ,doubling time + 2 sigma  , intercept log-lin fitting, 2 sigma intercept ,R^2]`. It can be saved into a .csv if `write_res=true`.
- The plots of the log-linear fitting and of the dynamics of specific growth rate if `do_plot=true`

<a name="ODE-one-well"></a>
## Fitting ODE function for one well
```julia
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
    calibration_OD_curve="NA" ,
 PopulationSize =100,
          maxiters = 10000,
           abstol = 0.001
    )
```
This function performs constrained parameter fitting on a single well's dataset using an ordinary differential equation (ODE) model. It estimates the model parameters within specified lower and upper bounds.

Arguments:

- `data::Matrix{Float64}`: The dataset with the growth curve, where the first row represents times, and the second row represents the variable to fit (e.g., OD), see [data formatting](#data).
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
- ` PopulationSize =100`: Size of the population of the optimization
-  ` maxiters = 10000`: stop criterion, the optimization is stopped when the number of iteration is bigger than `abstol`
- `abstol = 0.001`: stop criterion, the optimization is stopped when the loss is lesser than `abstol`


Output (if `results_ODE_fit =fitting_one_well_ODE_constrained(...)`:

- `results_ODE_fit[1]` an array with the following contents: `["name of model", "well", "param_1","param_2",..,"param_n","maximum specific gr using ode","maximum specific gr using data", "objective function value (i.e. loss of the solution)"]`
where ' "param_1","param_2",..,"param_n" ' are the parameter of the selected ODE as in this [table](#ODE_list)
- `results_ODE_fit[2]` the times of the fitted ODE
- `results_ODE_fit[3]` the numerical solution of the fitted ODE

- The plot of the  fitting  if `do_plot=true`

<a name="ODE-file"></a>
## Fitting ODE function for one file
```julia
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
    calibration_OD_curve="NA" ,
   PopulationSize =100,
          maxiters = 10000,
           abstol = 0.001 
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
- ` PopulationSize =100`: Size of the population of the optimization
-  ` maxiters = 10000`: stop criterion, the optimization is stopped when the number of iteration is bigger than `abstol`
- `abstol = 0.001`: stop criterion, the optimization is stopped when the loss is lesser than `abstol`

Output (if `results_ODE_fit =fit_file_ODE(...)`:

- `results_ODE_fit[1]` a matrix with the following contents for each column: `["name of model", "well", "param_1","param_2",..,"param_n","maximum specific gr using ode","maximum specific gr using data", "objective function value (i.e. loss of the solution)"]`. It can be saved into a .csv if `write_res=true`.
where ' "param_1","param_2",..,"param_n" ' are the parameter of the selected ODE as in this [table](#ODE_list)
- The plot of the  fitting  if `do_plot=true`
  
<a name="custom-ODE"></a>
## Fitting custom ODE function for one well
```julia
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
    calibration_OD_curve="NA"  ,
   PopulationSize =100,
          maxiters = 10000,
           abstol = 0.001
    )
```

This function is designed for fitting an ordinary differential equation (ODE) model to a dataset representing the growth curve of a microorganism in a well. It utilizes a customizable ODE model, optimization methods, and integration techniques for parameter estimation.

 Arguments:

- `data::Matrix{Float64}`: The dataset with the growth curve, where the first row represents times, and the second row represents the variable to fit (e.g., OD), see [data formatting](#data).
- `name_well::String`: The name of the well.
- `label_exp::String`: The label of the experiment.
- `model::Any`: The ODE model to use.
- `lb_param::Vector{Float64}`: Lower bounds for the parameters.
- `ub_param::Vector{Float64}`: Upper bounds for the parameters.
- `n_equation::Int`: The number of ODEs in the system.

Output (if `results_ODE_fit =fitting_one_well_custom_ODE(...)`:

- `results_ODE_fit[1]` a matrix with the following contents for each column: `["name of model", "well", "param_1","param_2",..,"param_n","maximum specific gr using ode","maximum specific gr using data", "objective function value (i.e. loss of the solution)"]`. It can be saved into a .csv if `write_res=true`.
where ' "param_1","param_2",..,"param_n" ' are the parameter of the selected ODE as in this [table](#ODE_list)
- The plot of the  fitting  if `do_plot=true`
  
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
- ` PopulationSize =100`: Size of the population of the optimization
-  ` maxiters = 10000`: stop criterion, the optimization is stopped when the number of iteration is bigger than `abstol`
- `abstol = 0.001`: stop criterion, the optimization is stopped when the loss is lesser than `abstol`



<a name="Sensitivity-analysis"></a>
## Sensitivity analysis
```julia
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
    calibration_OD_curve="NA"  ,
   PopulationSize =100,
          maxiters = 10000,
           abstol = 0.001
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
- ` PopulationSize =100`: Size of the population of the optimization
-  ` maxiters = 10000`: stop criterion, the optimization is stopped when the number of iteration is bigger than `abstol`
- `abstol = 0.001`: stop criterion, the optimization is stopped when the loss is lesser than `abstol`

Output (if `results_ODE_morris_sensitivity =one_well_morris_sensitivity(...)`:
- `results_ODE_morris_sensitivity[1]` a with in each column the initial guess for the parameters of the optimization in the same order of [table](#ODE_list)
- `results_ODE_morris_sensitivity[2]` a matrix with the following contents for each column: `["name of model", "well", "param_1","param_2",..,"param_n","maximum specific gr using ode","maximum specific gr using data", "objective function value (i.e. loss of the solution)"]`. It can be saved into a .csv if `write_res=true`.
where ' "param_1","param_2",..,"param_n" ' are the parameter of the selected ODE as in this [table](#ODE_list)
- The plot of the  fitting  if `do_plot=true`

<a name="Model-selection"></a>

## Model selection
```julia
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
    verbose=false,
   PopulationSize =100,
          maxiters = 10000,
           abstol = 0.001
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
- ` PopulationSize =100`: Size of the population of the optimization
-  ` maxiters = 10000`: stop criterion, the optimization is stopped when the number of iteration is bigger than `abstol`
- `abstol = 0.001`: stop criterion, the optimization is stopped when the loss is lesser than `abstol`

Output (if `Model_selection =ODE_Model_selection(...)`:

- `Model_selection[1]` a Matrix containing the loss and the AIC score for each model.
- `Model_selection[2]` a Tuple containing all the fitted models.
- `Model_selection[3]` the AIC score of the best model
- `Model_selection[4]` , the loss of the best model
- `Model_selection[5]` , the parameter of the best model
- `Model_selection[6]` , the string of the best model
- `Model_selection[7]` , the numerical solution of the fitted ODE
- The plot of the  fitting of the best model if `do_plot=true`


<a name="cdp"></a>
## Change point detection
```julia
cpd_local_detection(data::Matrix{Float64},
    n_max_cp::Int;
    type_of_detection="lsdd",
    type_of_curve="original", 
    pt_derivative = 0,
    size_win =2,
method= "peaks_prominence",
number_of_bin = 40)

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
- `method = "peaks_prominence"` : method to detect peak on the dissimilarity curve. Option "peaks_prominence" use prominece of peaks to score them. `"thr_scan"` grid scan with a threshold to detect peaks.
- `number_of_bin = 40`: number of bins for the grid search. used only if `method = "thr_scan"`

Output (if `cdps =cpd_local_detection(...)`:

- `cdps[1]` Indexes of the detected change points

<a name="cdp-fixed"></a>
## Fitting segmented ODE with fixed change-point number
```julia
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
    beta_smoothing_ms = 2.0,
method_peaks_detection= "peaks_prominence",
n_bins = 40,
   PopulationSize =100,
          maxiters = 10000,
           abstol = 0.001
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
- `method_peaks_detection = "peaks_prominence"` : method to detect peak on the dissimilarity curve. Option "peaks_prominence" use prominece of peaks to score them. `"thr_scan"` grid scan with a threshold to detect peaks.
- `n_bins = 40`: number of bins for the grid search. used only if `method_peaks_detection = "thr_scan"`
- ` PopulationSize =100`: Size of the population of the optimization
-  ` maxiters = 10000`: stop criterion, the optimization is stopped when the number of iteration is bigger than `abstol`
- `abstol = 0.001`: stop criterion, the optimization is stopped when the loss is lesser than `abstol`

Output (if `res =selection_ODE_fixed_change_points(...)`:

- `res[1]`. Parameters of each segment
- `res[2]`. Interval of the ODE segment
- `res[3]`. Time of the fitted solution
- `res[4]`. Numerical fitted solution
- The plot of the  fitting if `do_plot=true`


<a name="cdp-search"></a>
## Fitting segmented ODE with direct search for a maximum number of change-points 
```julia
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
    type_of_curve = "original", 
    pt_avg = 1 , 
    smoothing= true, 
    do_plot=false, 
    path_to_plot="NA", 
    path_to_results="NA",
    win_size=2, 
    pt_smooth_derivative=7,
    penality_parameter=2.0,
    multiple_scattering_correction="false", 
    calibration_OD_curve="NA",  
   save_all_model=false,
    method_peaks_detection= "peaks_prominence",
    n_bins = 40,
   PopulationSize =100,
          maxiters = 10000,
           abstol = 0.001 )
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
- `pt_avg=1`: Number of points to generate the initial condition or do the rolling avg smoothing.
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
- `method_peaks_detection = "peaks_prominence"` : method to detect peak on the dissimilarity curve. Option "peaks_prominence" use prominece of peaks to score them. `"thr_scan"` grid scan with a threshold to detect peaks.
- `n_bins = 40`: number of bins for the grid search. used only if `method_peaks_detection = "thr_scan"`
- ` PopulationSize =100`: Size of the population of the optimization
-  ` maxiters = 10000`: stop criterion, the optimization is stopped when the number of iteration is bigger than `abstol`
- `abstol = 0.001`: stop criterion, the optimization is stopped when the loss is lesser than `abstol`


Output (if `res =ODE_selection_NMAX_change_points(...)`:

- `res[1]`. The  parameters of each segment of the top model
- `res[2]`. Time of the fitted solution
- `res[3]`. Numerical fitted solution
- The plot of the  fitting if `do_plot=true`
- If `save_all_model=true` the best segnmented model is saved for each number of change points between 0 and `n_max_change_points`


