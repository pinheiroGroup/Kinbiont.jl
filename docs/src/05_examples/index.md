# [Examples and Tutorial](@id examples)

This section provides some copy-and-paste examples of JMAKi.jl

1. [ Simulating/Loading single kinetics data](#simulating-data)
   -[Loading  data from a .csv](#loading-data-ODE)
   -[Simulating Data with ODEs](#simulating-data-ODE)
   -[Simulating Data with stochastic simulations](#simulating-data-stochastic)
2. [Data Preprocessing](#data-preprocessing)
3. [Fitting a single kinetics](#model-fitting)
    - [Log-Lin fitting](#fitting-log-lin)
    - [Fitting ODE Models](#fitting-ode)
    - [Custom ODE Fitting](#custom-ode-fitting)
    - [Sensitivity Analysis](#sensitivity-analysis)
    - [ODE Model Selection](#model-selection)
    - [Fitting NL Models](#fitting-nl)
    - [Custom NL Fitting](#custom-nl-fitting)
    - [NL Sensitivity Analysis](#nl-sensitivity-analysis)
    - [NL Model Selection](#nl-model-selection)
    - [ODE segmentation](#ODE-segmented-single)
    - [NL segmentation](#NL-segmented-single)
4. [Fitting one file (a plate)](#model-fitting-plate)
    - [Plot one file](#plot-file)
    - [Log-Lin fitting](#fitting-log-lin-file)
    - [Fitting ODE Models](#fitting-ode-file)
    - [Fitting NL Models](#fitting-NL-file)
    - [ODE segmentation](#ODE-segmented)
    - [NL segmentation](#NL-segmented)

## Simulating/Loading single kinetics data
### Loading  data from a .csv
It is possible to load single curves using CSV package in Julia and format them to be compatible with JMAKi. In this example we suppose that the .csv file is formatted as required in JMAKi documentation.

```julia
df_data  =CSV.file("your_path/data.csv")
names_of_cols = propertynames(df_data)  
# assuming first column is time and we want to fit the second one 

data_OD = Matrix(hcat(df_data[names_of_cols[1]],df_data[names_of_cols[2]],))
```
### Simulating Data with ODEs

To simulate data using Ordinary Differential Equations (ODEs):

```julia
# Simulating data with an ODE
model = "triple_piecewise_adjusted_logistic"
n_start = [0.1]
tstart = 0.0
tmax = 600.0
delta_t = 10.0

param_of_ode = [0.06, 1.0, 200, 0.5, 0.001, 450, -0.0002]

# Calling the simulation function
sim = ODE_sim(model, n_start, tstart, tmax, delta_t, param_of_ode)

# Plotting scatterplot of data without noise
Plots.scatter(sim, xlabel="Time", ylabel="Arb. Units", label=["Data " nothing], color=:blue, size=(300, 300))

#adding uniform random noise
noise_unifom = rand(Uniform(-0.05,0.05),length(sim.t))


data_t = reduce(hcat,sim.t)
data_o = reduce(hcat,sim.u)
data_OD = vcat(data_t,data_o)
data_OD[2,:] = data_OD[2,:] .+ noise_unifom
# ploting scatterplot of data with noise

Plots.scatter!(data_OD[1,:],data_OD[2,:], xlabel="Time", ylabel="Arb. Units", label=["Data " nothing],color=:red,markersize =2 ,size = (300,300))

```

### Simulating Data with stochastic simulations

To simulate data using  stochastic models:

```julia
sim =stochastic_sim("Monod", #string of the model
   1, # number of starting cells
   10000.1, # starting # molecules (in "biomass units") of the limiting nutrients
    0.0, # start time of the sim
    2000.0, # final time of the sim
    0.1, # delta t for poisson approx
    11.1,
    10.1, # monod constant
    0.06, # massimum possible growth rate
   10.0, # lag time
    0.0000001,# nutrient consumed per division (conc)
    1000.0 #volume
)
plot(sim[3],sim[1], xlabel="Time",ylabel="# of indivuals")
plot(sim[3],sim[2], xlabel="Time",ylabel="nutrients/volume")

data_OD = Matrix(transpose(hcat(sim[3],sim[1])))
```

## Data Preprocessing
We start applying a rolling average smoothing to the data. In the example, a rolling window of size 7 is applied to the original data (data_OD generated in the previous examples). 
```julia
data_ODsmooth = smoothing_data(data_OD, 7)
data_ODsmooth = Matrix(data_ODsmooth)

# Plotting scatterplot of smoothed data
Plots.scatter(data_ODsmooth[1, :], data_ODsmooth[2, :], xlabel="Time", ylabel="Arb. Units", label=["Smoothed data " nothing], markersize=2, color=:blue, size=(300, 300))
```
Furthermore, to address potential external influences, a correction for multiple scattering is applied to the smoothed data. This correction is executed through the correction_OD_multiple_scattering function, requiring an external file (calibration_curve.csv).  it is  optional in the provided example. 
```julia

# Multiple scattering correction (optional, comment out if not needed)
data_ODsmooth = correction_OD_multiple_scattering(data_ODsmooth, "/your_path/calibration_curve.csv")

# Plotting scatterplot of preprocessed data
Plots.scatter(data_ODsmooth[1, :], data_ODsmooth[2, :], xlabel="Time", ylabel="Arb. Units", label=["Pre-processed data" nothing], markersize=2, color=:blue, size=(300, 300))
```

## Fitting single well
When data are store in a Matrix of Float64 with 2 rows and the number of columns 
### Log-Lin fitting

To perform the log-linear fitting  (of data generated in the previous examples) it is suffucient to run. 

```julia
res_log_lin = fitting_one_well_Log_Lin(
    data_OD, # dataset first row times second row OD
   "test", # name of the well
    "test log-lin fitting"; #label of the experiment
    display_plots=true, # do plots or no
    type_of_smoothing="rolling_avg", # type of smoothing
    pt_avg=7, # number of the point for rolling avg not used in the other cases
    pt_smoothing_derivative=7, # number of poits to smooth the derivative
    pt_min_size_of_win=7, # minimum size of the exp windows in number of smooted points
)
```
 The results are stored in the res_log_lin variable. 


###    Fitting ODE Models

Before fitting, the model and upper and lower bounds for the ODE parameters are defined 
```julia


model ="aHPM"
# Upper bounds of the parameters of the ODE
ub_ahpm = [1.2, 1.1, 2.0, 20]

# Lower bounds of the parameters of the ODE
lb_ahpm = [0.0001, 0.00000001, 0.00, 0]
```
The actual fitting is accomplished through the fitting_one_well_ODE_constrained function. This function takes the 
dataset (data_OD generated in the previous examples), the name and label of the well, the ODE model to use ("dHPM" in this case), as well as the upper and lower bounds for the ODE parameters. Additionally, the function allows for plotting the results (do_plot=true) and specifying the path to save the generated plots (path_to_plot=path_to_plotting).
```julia
# Performing ODE fitting
results_ODE_fit = fitting_one_well_ODE_constrained(
    data_OD, 
    "test",
    "test_ODE",
    model,
    lb_dhpm,
    ub_dhpm;
    display_plots=true, # display plots in julia or not
    smoothing=true, # the smoothing is done or not?
    pt_avg=3, # number of the points to do rolling avg
    PopulationSize=300,
    maxiters=2000000,

)

```
The results are stored in 'results_ODE_fit' with the following format
```julia
 results_ODE_fit = ["name of model", "well", "param_1","param_2",..,"param_n","maximum specific gr using ode","maximum specific gr using data", "objective function value (i.e. loss of the solution)"]
```

where ' "param_1","param_2",..,"param_n" ' are the parameter of the selected ODE.


###   Custom ODE Fitting

Using the custom ODE fitting, users can fit any ordinary differential equation. First it is necessary to define the differential equation. This is done by declaring a new function, e.g.:
```julia

# Custom ODE function
function ODE_custom(du, u, param, t)
    du[1] = u[1] * (1 - u[1]) * param[2] + param[1] * u[1]end

```
Now, we set upper and lower bounds of the parameters of the ODE:

```julia
# Bounds for the custom ODE parameters
custom_ub = [1.2, 1.1]
custom_lb = [0.0001, 0.00000001]
```
Finally, we perform the fit:
```julia

# Performing custom ODE fitting
results_ODE_fit = fitting_one_well_custom_ODE(
    data_OD, # dataset first row times second row OD
   "test", # name of the well
    "test custom ode", #label of the experiment
   ODE_custom, # ode model to use
    custom_lb, # lower bound param
    custom_ub, # upper bound param
    1; # number ode in the system
    optmizator=BBO_adaptive_de_rand_1_bin_radiuslimited(), # selection of optimization method
    integrator=Tsit5(), # selection of sciml integrator
    display_plots=true, # do plots or no
    pt_avg=3, # numebr of the point to generate intial condition
    smoothing=true, # the smoothing is done or not?
    PopulationSize=300,
    maxiters=2000000,
)
```
The results are stored in 'results_ODE_fit' with the same format of the previous examples.

###   Sensitivity Analysis

The sensitivity analysis is performed with the one_well_morris_sensitivity function. This function takes the preprocessed dataset (data_OD generated in the previous examples), the name and label of the well, the ODE model to use ("dHPM" in this case), as well as the lower and upper bounds for the ODE parameters. The number of steps in the Morris method (n_step_sensitivity) is specified to control the granularity of the analysis.

```julia
# Number of steps for Morris sensitivity analysis
n_step_sensitivity = 3

# Performing Morris sensitivity analysis
sensitivity_test = one_well_morris_sensitivity(
    data_OD, "test", "test_sensitivity", "dHPM", lb_dhpm, ub_dhpm,
    N_step_morris=n_step_sensitivity
)

```

###   ODE Model Selection


Several ODE models are considered for the selection process, each defined by upper (`list_ub`) and lower (`list_lb`) bounds for their respective parameters. The models include "dHPM," "piecewise_damped_logistic," "triple_piecewise," and "baranyi_roberts."

```julia
# Model candidates and their parameter bounds
list_of_models = ["dHPM", "piecewise_damped_logistic", "triple_piecewise", "baranyi_roberts"]

ub_piece_wise_logistic =[ 0.06 , 2.0 , 500.0 , 10.0 ,  0.001    ]
lb_piece_wise_logistic =[ 0.0001 , 0.001,0.0  , 0.001 ,  - 0.001  ]
ub_baranyi_roberts =[ 0.06 , 2.0 , 500.0 , 10.0,  10   ]
lb_baranyi_roberts =[ 0.0001 , 0.001, 0.0 ,  0.01 , 0  ]
ub_triple_exp =[ 1.2 , 0.001 ,  0.2 , 500.0  , 2000   ]
lb_triple_exp =[ 0.0001 , -0.001, 0.0  , 00.0 ,  200.0   ]
ub_dhpm =[ 1.2 , 1.1 , 2.0  ,20  ]
lb_dhpm =[ 0.0001 , 0.00000001, 0.00 ,0 ]

list_ub = [ub_dhpm, ub_piece_wise_logistic, ub_triple_exp, ub_baranyi_roberts]
list_lb = [lb_dhpm, lb_piece_wise_logistic, lb_triple_exp, lb_baranyi_roberts]
```


The model selection process is initiated with the `ODE_Model_selection` function. This function takes the preprocessed dataset (`data_OD` generated in the previous examples), the name and label of the well, the list of models, and their respective upper and lower bounds. Additionally, the function allows for plotting the results of the best-fit model (`plot_best_model=true`) and specifies the path to save the generated plots (`path_to_plot=path_to_plotting`). The `verbose=true` option provides detailed output during the model selection process.

```julia
# Performing model selection
results_ms = ODE_Model_selection(
    data_OD, "test", "test_model_selection",
    list_of_models, list_lb, list_ub,
    plot_best_model=true, path_to_plot=path_to_plotting, verbose=true
)
```




The results of the model selection process are stored in the `results_ms` variable.


## Fitting one file (a plate)
The next three functions work directly on a file. So in this case are mandatory the .csv of data and annotation (see  [Data and annotation formatting](#data) ). Aslo in XXXXX the user can download an examples of data and annotation.


### Plot one file
The provided code  is an example of plotting experimental data in a .csv file:


```julia

# Paths to data, annotation, results, and plots
path_to_data = "/example/data_channel_1.csv"
path_to_annotation = "/example/annotation_channel_1_media_M9 + 0.2% Glucose.csv"
path_to_plot = "/example/plots/"

plot_data(   "example", #label of the experiment
    path_to_data, # path to the folder to analyze
    path_to_annotation;# path to the annotation of the wells
    path_to_plot=path_to_plot, # path where to save Plots
    display_plots=true ,# display plots in julia or not
    save_plot=true, # save the plot or not
    overlay_plots=true, # true a single plot for all dataset false one plot per well
    blank_subtraction="avg_blank" # string on how to use blank (NO,avg_blank,time_avg)
)
```

### Log-Lin fitting

 The provided code  is an example of Log-Lin fitting of experimental data in a .csv file:


```julia

# Paths to data, annotation, results, and plots
path_to_data = "/example/data_channel_1.csv"
path_to_annotation = "/example/annotation_channel_1_media_M9 + 0.2% Glucose.csv"
path_to_results = "/example/results/"
path_to_plot = "/example/plots/"

res = fit_one_file_Log_Lin(
    "log_lin_WT_CHL_dose_reponse", #label of the experiment
    path_to_data, # path to the folder to analyze
    path_to_annotation;
  path_to_results = path_to_results,  # path where to save results
    path_to_plot = path_to_plot,        # path where to save plots
    do_plot = true,          # do and visualize the plots of data
    write_res = true,        # write results
    )
```

###   Fitting ODE Models
 The provided code  is an example of fitting a differential equation model to experimental data in a .csv file:

```julia
# Define upper and lower bounds for the parameters of the ODE model
ub_dhpm = [0.1, 0.1, 2.0, 5.0]
lb_dhpm = [0.001, 0.00001, 0.01, 0.5]

# Paths to data, annotation, results, and plots
path_to_data = "/example/data_channel_1.csv"
path_to_annotation = "/example/annotation_channel_1_media_M9 + 0.2% Glucose.csv"
path_to_results = "/example/results/"
path_to_plot = "/example/plots/"

# Fit the ODE model to the experimental data
res = fit_file_ODE(
    "WT_CHL_dose_reponse",  # label of the experiment
    path_to_data,            # path to the data
    path_to_annotation,      # path to the annotation of the wells
    "dHPM",                   # string of the used model
    lb_dhpm,                  # array of the lower bound of the parameters
    ub_dhpm;                 # array of the upper bound of the parameters
    path_to_results = path_to_results,  # path where to save results
    path_to_plot = path_to_plot,        # path where to save plots
    do_plot = true,          # do and visualize the plots of data
    write_res = true,        # write results
    pt_avg = 2,              # number of points to do smoothing average
    PopulationSize = 500,    # population size for optimization
    maxiters = 500000,       # maximum number of iterations
    abstol = 0.00000000001   # absolute tolerance for optimization
)
```

This example is fitting an ODE model (specifically the "dHPM" model) to experimental data provided in CSV files.

## ODE segmentation with fixed number of change points

In this example, we demonstrate the process of fitting a dataset with a sequence of ODEs using a segmentation approach. The dataset is generated with three segments, each modeled by a different ODE.
Then we fit it with the 'selection_ODE_fixed_change_points' function

 ```julia
# First segment ODE
model = "exponential"
n_start = [0.1]
tstart = 0.0
tmax = 10.0
delta_t = 2.0
param_of_ode = [0.00]
sim_1 = ODE_sim(model, n_start, tstart, tmax, delta_t, integrator, param_of_ode)
sol_1 = reduce(vcat, sim_1)

# Second segment ODE
model = "logistic"
n_start = [sol_1[end]]
tstart = 10.0
tmax = 80.0
delta_t = 2.0
param_of_ode = [0.2, 0.4]
sim_2 = ODE_sim(model, n_start, tstart, tmax, delta_t, integrator, param_of_ode)
sol_2 = reduce(vcat, sim_2)

# Third segment ODE
model = "logistic"
n_start = [sol_2[end]]
tstart = 80.0
tmax = 200.0
delta_t = 2.0
param_of_ode = [0.1, 0.8]
sim_3 = ODE_sim(model, n_start, tstart, tmax, delta_t, integrator, param_of_ode)
sol_3 = reduce(vcat, sim_3)

# Concatenating simulations
times_sim = vcat(sim_1.t, sim_2.t)
times_sim = vcat(times_sim, sim_3.t)
sol_sim = vcat(sol_1, sol_2)
sol_sim = vcat(sol_sim, sol_3)

# Plotting the generated dataset
Plots.scatter(sol_sim, xlabel="Time", ylabel="Arb. Units", label=["Data " nothing], color=:blue, size=(300, 300))

# Computing and visualizing the first derivative
data_OD = Matrix(transpose(hcat(times_sim, sol_sim)))
deriv = specific_gr_evaluation(data_OD, 0)
Plots.scatter(deriv, xlabel="Time", ylabel="Arb. Units", label=["Data " nothing], color=:blue, size=(300, 300))

# Adding uniform noise to the dataset
noise_uniform = rand(Uniform(-0.01, 0.01), length(sol_sim))
data_OD = Matrix(transpose(hcat(times_sim, sol_sim)))
data_OD[2, :] = data_OD[2, :] .+ noise_uniform

# Plotting the noisy dataset
Plots.scatter(data_OD[1, :], data_OD[2, :], xlabel="Time", ylabel="Arb. Units", label=["Data " nothing], color=:blue, markersize=2, size=(300, 300))

# Initializing all the models for selection
ub_exp = [0.1]
lb_exp = [-0.01]
ub_logistic = [0.9, 5.0]
lb_logistic = [0.0001, 0.001]
ub_hpm = [0.1, 20.0, 50.001]
lb_hpm = [0.0001, 0.000001, 0.001]
ub_hpm_exp = [0.1, 20.0]
lb_hpm_exp = [0.0001, 0.0000001]

list_of_models = ["exponential", "HPM", "HPM_exp", "logistic"]
list_ub_param = [ub_exp, ub_hpm, ub_hpm_exp, ub_logistic]
list_lb_param = [lb_exp, lb_hpm, lb_hpm_exp, lb_logistic]

# Fitting with a fixed number of change points
test_fixed_cdp = selection_ODE_fixed_change_points(
    data_OD, "test", "test", list_of_models, list_lb_param, list_ub_param, 2;
    do_plot=true, path_to_plot=path_to_plotting, pt_smooth_derivative=0
)
```

results are stored in test_fixed_cdp.

## ODE segmentation
Using the same code as the previous example to generate the data the  fit  is performed with 
```julia
# Fitting with direct search on the number of change points
 test_cdp = ODE_selection_NMAX_change_points(data_OD,
    "test",
    "test",
    list_lb_param,
    list_ub_param,
    list_of_models,
   3;
    do_plot=true,
 path_to_plot=path_to_plotting,
 pt_smooth_derivative=0)
```
results are stored in test_cdp.

