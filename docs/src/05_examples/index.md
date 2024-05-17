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
    - [ODE segmentation](#ODE-segmented-single)
    - [Fitting NL Models](#fitting-nl)
    - [Custom NL Fitting](#custom-nl-fitting)
    - [NL Sensitivity Analysis](#nl-sensitivity-analysis)
    - [NL Model Selection](#nl-model-selection)
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

## Fitting a single kinetics
When data are store in a Matrix of Float64 with 2 rows it is possible to perfor various analyses

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
The actual fitting is accomplished through the fitting_one_well_ODE_constrained function. 
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

###   ODE Sensitivity Analysis

The sensitivity analysis is performed with the one_well_morris_sensitivity function. This function takes the preprocessed dataset (data_OD generated in the previous examples), the name and label of the well, the ODE model to use ("aHPM" in this case), as well as the lower and upper bounds for the ODE parameters. The number of steps in the Morris method (n_step_sensitivity) should be specified.

```julia
# Number of steps for Morris sensitivity analysis
n_step_sensitivity = 5

# Performing Morris sensitivity analysis
sensitivity_test = one_well_morris_sensitivity(
    data_OD, "test", "test_sensitivity", "dHPM", lb_dhpm, ub_dhpm,
    N_step_morris=n_step_sensitivity ;
    optmizator=BBO_adaptive_de_rand_1_bin_radiuslimited(), # selection of optimization method
    integrator=Tsit5(), # selection of sciml integrator
    pt_avg=3, # numebr of the point to generate intial condition
    smoothing=true, # the smoothing is done or not?
    type_of_smoothing="rolling_avg",
    PopulationSize=300,
    maxiters=2000000,
    abstol=0.00001,
)

```
In this case the user will obtain the results of the fit for each of the Morris step and a matrix containg the coordinate of the hyperspace of the starting guess assicioted to each fit.
###   ODE Model Selection

The user can use the model selection function to use AIC (or AICc) to identify the best model to fit a specific kinetics.
We start defining the list of the models and the upper and lower bounds of the parameters of each one:

```julia
# Model candidates and their parameter bounds
list_of_models = ["aHPM", "piecewise_damped_logistic", "triple_piecewise", "baranyi_roberts"]

ub_piece_wise_logistic =[ 0.06 , 2.0 , 500.0 , 10.0 ,  0.001    ]
lb_piece_wise_logistic =[ 0.0001 , 0.001,0.0  , 0.001 ,  - 0.001  ]
ub_baranyi_roberts =[ 0.06 , 2.0 , 500.0 , 10.0,  10   ]
lb_baranyi_roberts =[ 0.0001 , 0.001, 0.0 ,  0.01 , 0  ]
ub_triple_exp =[ 1.2 , 0.001 ,  0.2 , 500.0  , 2000   ]
lb_triple_exp =[ 0.0001 , -0.001, 0.0  , 00.0 ,  200.0   ]
ub_ahpm =[ 1.2 , 1.1 , 2.0  ,20  ]
lb_ahpm =[ 0.0001 , 0.00000001, 0.00 ,0 ]

list_ub = [ub_ahpm, ub_piece_wise_logistic, ub_triple_exp, ub_baranyi_roberts]
list_lb = [lb_ahpm, lb_piece_wise_logistic, lb_triple_exp, lb_baranyi_roberts]
```


The model selection process is runned with the `ODE_Model_selection` function. 

```julia
# Performing model selection
results_ms = ODE_Model_selection(
    data_OD,
    "test", 
    "test_model_selection",
    list_of_models,
    list_lb,
    list_ub;
    optmizator=BBO_adaptive_de_rand_1_bin_radiuslimited(), # selection of optimization method
    integrator=KenCarp4(autodiff=true), # selection of sciml integrator
    pt_avg=3, # number of the point to generate intial condition
    beta_penality=2.0, # penality for AIC evaluation
    smoothing=true, # the smoothing is done or not?
    type_of_smoothing="rolling_avg",
    type_of_loss="L2", # type of used loss
    display_plot_best_model=true, # results of the best fit to be plotted
    PopulationSize=300,
    maxiters=2000000,
    abstol=0.00001,
    correction_AIC=false,)
```


The results of the model selection process are stored in the `results_ms` variable.

### ODE segmentation
For a single kinetics it is possible to run the segmentaion in two different ways. A manual selection of the change points or a using a change points detection algorithm to find them.

First of all we generate a synthetic daset composed by more than one model  (skip this part if you import a real dataset in the data variable).


```julia

# first segment ODE

model = "logistic"
n_start =[0.1]
tstart =0.0

tmax = 0100.0
delta_t=2.0
param_of_ode= [0.1,0.2]

sim_1 = ODE_sim(model, #string of the model
    n_start, # starting condition
    tstart, # start time of the sim
    tmax, # final time of the sim
    delta_t, # delta t for poisson approx
    param_of_ode # parameters of the ODE model
)

sol_1 =reduce(vcat,sim_1)

# second segment ODE


model = "logistic"
n_start =[sol_1[end]]
tstart =100.0
tmax = 0200.0
delta_t=2.0
param_of_ode= [0.2,0.5]


sim_2= ODE_sim(model, #string of the model
    n_start, # starting condition
    tstart, # start time of the sim
    tmax, # final time of the sim
    delta_t, # delta t for poisson approx
    param_of_ode # parameters of the ODE model
)

sol_2 =reduce(vcat,sim_2)
# third segment ODE


model = "logistic"
n_start =[sol_2[end]]
tstart =200.0
tmax = 0300.0
delta_t=2.0
param_of_ode= [0.1,0.9]

sim_3= ODE_sim(model, #string of the model
    n_start, # starting condition
    tstart, # start time of the sim
    tmax, # final time of the sim
    delta_t, # delta t for poisson approx
    param_of_ode # parameters of the ODE model
)

sol_3 =reduce(vcat,sim_3)
times_sim =vcat(sim_1.t,sim_2.t)
times_sim =vcat(times_sim,sim_3.t)

# binding the simulatios
sol_sim =vcat(sol_1,sol_2)
sol_sim =vcat(sol_sim,sol_3)


Plots.scatter(sol_sim, xlabel="Time", ylabel="Arb. Units", label=["Data " nothing],  color=:blue, size = (300,300))

Plots.scatter(sol_sim[40:60], xlabel="Time", ylabel="Arb. Units", label=["Data " nothing],  color=:blue, size = (300,300))


data_OD = Matrix(transpose(hcat(times_sim,sol_sim)))

# Adding uniform noise to the dataset
noise_uniform = rand(Uniform(-0.01, 0.01), length(sol_sim))
data_OD = Matrix(transpose(hcat(times_sim, sol_sim)))
data_OD[2, :] = data_OD[2, :] .+ noise_uniform

# Plotting the noisy dataset
Plots.scatter(data_OD[1, :], data_OD[2, :], xlabel="Time", ylabel="Arb. Units", label=["Data " nothing], color=:blue, markersize=2, size=(300, 300))

```
Then we generate the list of models that will be used and their parameters's lower/upper bounds:
```julia
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

```

First, we fit giving to JMAKi the list of change points:

```julia
cdp_list = [100.0, 200.0]

res = selection_ODE_fixed_intervals(
   data_OD, # dataset first row times second row OD
    "test", # name of the well
    "test segmentation ODE", #label of the experiment
    list_of_models, # ode models to use
    list_lb_param, # lower bound param
    list_ub_param, # upper bound param
    cdp_list;
    type_of_loss="L2", # type of used loss
    optmizator=BBO_adaptive_de_rand_1_bin_radiuslimited(), # selection of optimization method
    integrator=Tsit5(), # selection of sciml integrator
    smoothing=true,
    type_of_smoothing="rolling_avg",
    pt_avg=3,
    display_plots=true,
    path_to_plot="NA", # where save plots
    pt_smooth_derivative=0,
    PopulationSize=300,
    maxiters=2000000,
    abstol=0.0000000001,
    correction_AIC=false)
```
Finally we can run a cpd algorithm and perfom the fitting:
```julia
n_change_points =2
segmentation_ODE(
    data_OD, # dataset first row times second row OD
    "test", # name of the well
    "test segmentation ODE", #label of the experiment
    list_of_models, # ode model to use
    list_lb_param, # lower bound param
    list_ub_param, # upper bound param
    n_change_pointst;
    detect_number_cpd=false,
    fixed_cpd=true,
    optmizator=BBO_adaptive_de_rand_1_bin_radiuslimited(), # selection of optimization method
    integrator=Tsit5(), # selection of sciml integrator
    type_of_loss="L2", # type of used loss
    type_of_detection="slinding_win",
    type_of_curve="original",
    pt_avg=3, # number of the point to generate intial condition
    smoothing=true, # the smoothing is done or not?
    display_plot=true, # do plots or no
    win_size=10, #  
    pt_smooth_derivative=0,
    PopulationSize=300,
    maxiters=2000000,
    abstol=0.00001,
    type_of_smoothing="rolling_average",
    thr_lowess=0.05,
    correction_AIC=true)

```
### Fitting NL Models
(we should discuss about this in theory the model selection functio can do all the stuffs except segmentation) With JKMAKi it is possible to fit any non-linear model this can be done by calling the function NL_model_selection in different ways.


First we declare upper and lower bound and the model (note that in this case we use array of array because the input can be more than one model)

```julia
nl_lb_1 = [0.0001 , 0.00000001, 0.00,0.0 ]
nl_ub_1 = [2.0001 , 10.00000001, 5.00,5.0]
list_models_f = ["NL_Richards"]
list_lb =[nl_lb_1]
list_ub = [nl_ub_1]
```


To perform a single fit on a time series (i.e., data_OD) then we run the following specifying method_of_fitting = "single_fit":
```julia
 NL_model_selection(data_OD, # dataset first row times second row OD
  "test", # name of the well
    "test NL fit", #label of the experiment
    list_models_f, #  model to use
    list_lb, # lower bound param
    list_ub; # upper bound param
    method_of_fitting="single_fit",
    list_u0=list_lb .+ (list_ub .- list_lb) ./ 2,# initial guess param
    optmizator=BBO_adaptive_de_rand_1_bin_radiuslimited(),
    display_plots=true, # display plots in julia or not
    pt_avg=3, # numebr of the point to generate intial condition
    smoothing=true, # the smoothing is done or not?
    type_of_smoothing="rolling_avg",
    PopulationSize=300,
    maxiters=2000000,
    abstol=0.00001,
)

```


The user can specify the intial guess list_u0 to improve the convergence of the fit. Otherwise it is possible to automatically find a good guess using a Markov Chain Monte Carlo restart (method_of_fitting ="MCMC"). 
This it is done running the following:

```julia
 NL_model_selection(data_OD, # dataset first row times second row OD
  "test", # name of the well
    "test NL fit", #label of the experiment
    list_models_f, #  model to use
    list_lb, # lower bound param
    list_ub; # upper bound param
    method_of_fitting="MCMC",
    nrep = 100,
    optmizator=BBO_adaptive_de_rand_1_bin_radiuslimited(),
    display_plots=true, # display plots in julia or not
    pt_avg=3, # numebr of the point to generate intial condition
    smoothing=true, # the smoothing is done or not?
    type_of_smoothing="rolling_avg",
    PopulationSize=300,
    maxiters=2000000,
    abstol=0.00001,
)

```

Alternately, the user can  opt to a Bootstrap approach (method_of_fitting ="Bootstrap"): 



```julia
 NL_model_selection(data_OD, # dataset first row times second row OD
  "test", # name of the well
    "test NL fit", #label of the experiment
    list_models_f, #  model to use
    list_lb, # lower bound param
    list_ub; # upper bound param
    method_of_fitting="Bootstrap",
    nrep = 100,
    optmizator=BBO_adaptive_de_rand_1_bin_radiuslimited(),
    display_plots=true, # display plots in julia or not
    pt_avg=3, # numebr of the point to generate intial condition
    smoothing=true, # the smoothing is done or not?
    type_of_smoothing="rolling_avg",
    PopulationSize=300,
    maxiters=2000000,
    abstol=0.00001,
    size_bootstrap=0.7,
)

```
If one wants to specify a custom NL function,he should start declaring the function:

```julia
function NL_model_exp(p, times)
    model = p[1] .* exp.(times .* p[2])
    return model
end


nl_ub =  [2.0001 , 10.00000001]
nl_lb =  [0.0001 , 0.00000001 ]


list_models_f = [NL_model_exp]
list_lb =[nl_lb]
list_ub = [nl_ub]
```

after this you can just supply this variables to the previous function, e.g.:

```julia
 NL_model_selection(data_OD, # dataset first row times second row OD
  "test", # name of the well
    "test NL fit", #label of the experiment
    list_models_f, #  model to use
    list_lb, # lower bound param
    list_ub; # upper bound param
    method_of_fitting="single_fit",
    list_u0=list_lb .+ (list_ub .- list_lb) ./ 2,# initial guess param
    optmizator=BBO_adaptive_de_rand_1_bin_radiuslimited(),
    display_plots=true, # display plots in julia or not
    pt_avg=3, # numebr of the point to generate intial condition
    smoothing=true, # the smoothing is done or not?
    type_of_smoothing="rolling_avg",
    PopulationSize=300,
    maxiters=2000000,
    abstol=0.00001,
)
```


### NL Sensitivity Analysis

As for the ODE model, also for the NL fit it is possible to perform a sensitivity analysis with respect to the initial starting guess of the parameters, in this case we just use method_of_fitting= "Morris_sensitivity" and nrep as the number of Morris steps


```julia
 NL_model_selection(data_OD, # dataset first row times second row OD
    "test", # name of the well
    "test NL fit", #label of the experiment
    list_models_f, #  model to use
    list_lb, # lower bound param
    list_ub; # upper bound param
    method_of_fitting="Morris_sensitivity",
    optmizator=BBO_adaptive_de_rand_1_bin_radiuslimited(),
    display_plots=true, # display plots in julia or not
    pt_avg=3, # numebr of the point to generate intial condition
    smoothing=true, # the smoothing is done or not?
    type_of_smoothing="rolling_avg",
    PopulationSize=300,
    maxiters=2000000,
    abstol=0.00001,
    nrep =20,
)
```



### NL Model Selection
To perform model selection we just specify an array with the list of all the used models:
```julia

nl_ub_1 =  [2.0001 , 10.00000001, 500.00]
nl_lb_1 =  [0.0001 , 0.00000001, 0.00 ]

nl_ub_2 =  [2.0001 , 10.00000001, 5.00,5.0]
nl_lb_2 =  [0.0001 , 0.00000001, 0.00,0.0 ]

nl_ub_3 =  [2.0001 , 10.00000001]
nl_lb_3 =  [0.0001 , 0.00000001]

nl_ub_4 =  [2.0001 , 10.00000001, 500.00]
nl_lb_4 =  [0.0001 , 0.00000001, 0.00 ]

list_models_f = ["NL_Gompertz","NL_Bertalanffy","NL_exponential","NL_Gompertz"]
list_lb =[nl_lb_1,nl_lb_2,nl_lb_3,nl_lb_4]
list_ub = [nl_ub_1,nl_ub_2,nl_ub_3,nl_ub_4]

```
and the we perform the fit, tuning the AIC parameter (beta_param) :
```julia

 NL_model_selection(data_OD, # dataset first row times second row OD
  "test", # name of the well
    "test NL fit", #label of the experiment
    list_models_f, #  model to use
    list_lb, # lower bound param
    list_ub; # upper bound param
    method_of_fitting="MCMC",
    nrep = 100,
    optmizator=BBO_adaptive_de_rand_1_bin_radiuslimited(),
    display_plots=true, # display plots in julia or not
    pt_avg=3, # numebr of the point to generate intial condition
    smoothing=true, # the smoothing is done or not?
    type_of_smoothing="rolling_avg",
    PopulationSize=300,
    maxiters=2000000,
    abstol=0.00001,
    beta_param=2.0,
    correction_AIC =false,
)
```

### NL segmentation

As the ODE case, for NL segmentation we can perform a segmente fit with manual selection of the change points or a using a change points detection algorithm to find them.
Note that in the case of the NL fit, the fit of the segment goes backward and since it is not possible to specificy the bounduary condition of the fit the user can force the optimization problem to mantain the continuty between segment with the parameter  penality_CI=8.0. 

we specify the list of used models and parameters bounds:

```julia

nl_ub_1 =  [2.0001 , 10.00000001, 500.00]
nl_lb_1 =  [0.0001 , 0.00000001, 0.00 ]

nl_ub_2 =  [2.0001 , 10.00000001, 5.00,5.0]
nl_lb_2 =  [0.0001 , 0.00000001, 0.00,0.0 ]

nl_ub_3 =  [2.0001 , 10.00000001]
nl_lb_3 =  [0.0001 , 0.00000001]

nl_ub_4 =  [2.0001 , 10.00000001, 500.00]
nl_lb_4 =  [0.0001 , 0.00000001, 0.00 ]

list_models_f = ["NL_Gompertz","NL_Bertalanffy","NL_exponential","NL_Gompertz"]
list_lb =[nl_lb_1,nl_lb_2,nl_lb_3,nl_lb_4]
list_ub = [nl_ub_1,nl_ub_2,nl_ub_3,nl_ub_4]

```


As before we specify the intervals of the segment and we proceed to fit:


```julia
cdp_list = [100.0, 200.0]

selection_NL_fixed_interval(
    data_OD, # dataset first row times second row OD
   "test ", # name of the well
    "test NL segmentation" , #label of the experiment
    list_models_f, # ode models to use
    list_lb, # lower bound param
    list_ub, # upper bound param
    cdp_list;
    type_of_loss="L2", # type of used loss
    optmizator=BBO_adaptive_de_rand_1_bin_radiuslimited(), # selection of optimization method
    method_of_fitting="MCMC", # selection of sciml integrator
    nrep=100,
    display_plots=true,
    beta_smoothing_ms=2.0, #  parameter of the AIC penality
    PopulationSize=300,
    maxiters=2000000,
    abstol=0.000000001,
    penality_CI=8.0,
    correction_AIC=true,
)

```
Another option could be to  run a cpd algorithm and perfom the fitting:




```julia
n_change_points = 2

selection_NL_max_change_points(
    data_OD, # dataset first row times second row OD
   "test ", # name of the well
    "test NL segmentation" , #label of the experiment
    list_models_f, # ode models to use
    list_lb, # lower bound param
    list_ub, # upper bound param
    n_change_points;
    optmizator=BBO_adaptive_de_rand_1_bin_radiuslimited(), # selection of optimization method
    method_of_fitting="MCMC", # selection of sciml integrator
    nrep=100,
    display_plots=true,
    win_size=7, # numebr of the point to generate intial condition
    pt_smooth_derivative=0,
    PopulationSize=300,
    maxiters=2000000,
    abstol=0.000000001,
    detect_number_cpd=false,
    fixed_cpd=true,
    penality_CI=8.0,
)

```
## Fitting a .csv file

Instead fitting a single kinetics the user can supply a ".csv" file (formatted as described in the section), and JMAKi will proceed to perform all the analysis on all the wells of the experiment. Note that the user can supply a annotation .csv in this case becomes possible to subtract the blanks and fit the average of replicates.

To use the following functions the user should input to JMAKi variables that contains the string of the path to the .csv files:

```julia
path_to_data = "your_path_to_data/data.csv"
path_to_annotation ="your_path_to_annotation/annotation.csv"
```
In the following examples we assume that these two variables have a value.


### Plot one file

It is possible to plot or display the plot of an experiment with the following:
```julia

 plot_data(
   "test", #label of the experiment
    path_to_data; # path to the folder to analyze
    path_to_annotation = path_to_annotation,# path to the annotation of the wells
    display_plots=true,# display plots in julia or noy
    save_plot = false,
    overlay_plots=true, # true a single plot for all dataset false one plot per well
    do_blank_subtraction="avg_blank", # string on how to use blank (NO,avg_subtraction,time_avg)
    )
```

### Log-Lin fitting

```julia
fit_one_file_Log_Lin(
    label_exp::String, #label of the experiment
    path_to_data::String; # path to the folder to analyze
    path_to_annotation::Any = missing,# path to the annotation of the wells
    path_to_results="NA",# path where save results
    path_to_plot="NA",# path where to save Plots
    display_plots=true,# display plots in julia or not
    save_plots=false, # save the plot or not    verbose=false, # 1 true verbose
    write_res=false, # write results
    type_of_smoothing="rolling_avg", # option, NO, gaussian, rolling avg
    pt_avg=7, # number of points to do smoothing average
    pt_smoothing_derivative=7, # number of poits to smooth the derivative
    pt_min_size_of_win=7, # minimum size of the exp windows in number of smooted points
    type_of_win="maximum", # how the exp. phase win is selected, "maximum" of "global_thr"
    threshold_of_exp=0.9, # threshold of growth rate in quantile to define the exp windows
    do_blank_subtraction="avg_blank", # string on how to use blank (NO,avg_subtraction,time_avg)
    avg_replicate=false, # if true the average between replicates is fitted. If false all replicate are fitted indipendelitly
    correct_negative="thr_correction", # if "thr_correction" it put a thr on the minimum value of the data with blank subracted, if "blank_correction" uses blank distrib to impute negative values
    thr_negative=0.01, # used only if correct_negative == "thr_correction"
    multiple_scattering_correction=false, # if true uses the given calibration curve to fix the data
    method_multiple_scattering_correction="interpolation",
    calibration_OD_curve="NA", #  the path to calibration curve to fix the data
    thr_lowess=0.05, # keyword argument of lowees smoothing
    verbose=false,
    blank_value = 0.0,
    blank_array = [0.0],)
```
### Fitting ODE Models
```julia

 fit_file_ODE(
    label_exp::String, #label of the experiment
    path_to_data::String, # path to the folder to analyze
    model::String, # string of the used model
    lb_param::Vector{Float64},# array of the array of the lower bound of the parameters
    ub_param::Vector{Float64}; # array of the array of the upper bound of the parameters
    path_to_annotation::Any = missing,# path to the annotation of the wells
    optmizator=BBO_adaptive_de_rand_1_bin_radiuslimited(), # selection of optimization method 
    integrator=Tsit5(), # selection of sciml integrator
    path_to_results="NA", # path where save results
    path_to_plot="NA", # path where to save Plots
    loss_type="RE", # string of the type of the used loss
    smoothing=false, # 1 do smoothing of data with rolling average
    type_of_smoothing="lowess",
    display_plots=true,# display plots in julia or not
    save_plots=false,
    verbose=false, # 1 true verbose
    write_res=false, # write results
    pt_avg=1, # number of points to do smoothing average
    pt_smooth_derivative=7, # number of points to do ssmooth_derivative
    do_blank_subtraction="avg_blank", # string on how to use blank (NO,avg_subtraction,time_avg)
    avg_replicate=false, # if true the average between replicates is fitted. If false all replicate are fitted indipendelitly
    correct_negative="thr_correction", # if "thr_correction" it put a thr on the minimum value of the data with blank subracted, if "blank_correction" uses blank distrib to impute negative values
    thr_negative=0.01,  # used only if correct_negative == "thr_correction"
    multiple_scattering_correction=false, # if true uses the given calibration curve to fix the data
    method_multiple_scattering_correction="interpolation",
    calibration_OD_curve="NA",  #  the path to calibration curve to fix the data
    PopulationSize=300,
    maxiters=2000000,
    abstol=0.00001,
    thr_lowess=0.05,
    blank_value = 0.0,
    blank_array = [0.0],
)
```

```julia

 fit_file_custom_ODE(
    label_exp::String, #label of the experiment
    path_to_data::String, # path to the folder to analyze
    model::Any, # string of the used model
    lb_param::Vector{Float64},# array of the array of the lower bound of the parameters
    ub_param::Vector{Float64}, # array of the array of the upper bound of the parameters
    n_equation::Int;
    path_to_annotation::Any = missing,# path to the annotation of the wells
    optmizator=BBO_adaptive_de_rand_1_bin_radiuslimited(), # selection of optimization method 
    integrator=Tsit5(), # selection of sciml integrator
    path_to_results="NA", # path where save results
    path_to_plot="NA", # path where to save Plots
    loss_type="RE", # string of the type of the used loss
    smoothing=false, # 1 do smoothing of data with rolling average
    type_of_smoothing="lowess",
    display_plots=true,# display plots in julia or not
    save_plots=false,
    verbose=false, # 1 true verbose
    write_res=false, # write results
    pt_avg=1, # number of points to do smoothing average
    pt_smooth_derivative=7, # number of points to do ssmooth_derivative
    do_blank_subtraction="avg_blank", # string on how to use blank (NO,avg_subtraction,time_avg)
    avg_replicate=false, # if true the average between replicates is fitted. If false all replicate are fitted indipendelitly
    correct_negative="thr_correction", # if "thr_correction" it put a thr on the minimum value of the data with blank subracted, if "blank_correction" uses blank distrib to impute negative values
    thr_negative=0.01,  # used only if correct_negative == "thr_correction"
    multiple_scattering_correction=false, # if true uses the given calibration curve to fix the data
    method_multiple_scattering_correction="interpolation",
    calibration_OD_curve="NA",  #  the path to calibration curve to fix the data
    PopulationSize=300,
    maxiters=2000000,
    abstol=0.00001,
    thr_lowess=0.05,
    blank_value = 0.0,
    blank_array = [0.0],
)
```
```julia

 ODE_model_selection_file(
    label_exp::String, #label of the experiment
    path_to_data::String, # path to the folder to analyze
    models_list::Vector{String}, # ode model to use 
    lb_param_array::Any, # lower bound param
    ub_param_array::Any; # upper bound param
    path_to_annotation::Any = missing,# path to the annotation of the wells
    optmizator=BBO_adaptive_de_rand_1_bin_radiuslimited(), # selection of optimization method 
    integrator=Tsit5(), # selection of sciml integrator
    path_to_results="NA", # path where save results
    path_to_plot="NA", # path where to save Plots
    loss_type="L2", # string of the type of the used loss
    smoothing=false, # 1 do smoothing of data with rolling average
    type_of_smoothing="lowess",
    display_plot_best_model=false, # one wants the results of the best fit to be plotted
    save_plot_best_model=false,
    beta_penality=2.0, # penality for AIC evaluation
    verbose=false, # 1 true verbose
    write_res=false, # write results
    pt_avg=1, # number of points to do smoothing average
    pt_smooth_derivative=7, # number of points to do ssmooth_derivative
    do_blank_subtraction="avg_blank", # string on how to use blank (NO,avg_subtraction,time_avg)
    avg_replicate=false, # if true the average between replicates is fitted. If false all replicate are fitted indipendelitly
    correct_negative="thr_correction", # if "thr_correction" it put a thr on the minimum value of the data with blank subracted, if "blank_correction" uses blank distrib to impute negative values
    thr_negative=0.01,  # used only if correct_negative == "thr_correction"
    multiple_scattering_correction=false, # if true uses the given calibration curve to fix the data
    method_multiple_scattering_correction="interpolation",
    calibration_OD_curve="NA",  #  the path to calibration curve to fix the data
    PopulationSize=300,
    maxiters=2000000,
    abstol=0.00001,
    thr_lowess=0.05,
    correction_AIC=true,
    blank_value = 0.0,
    blank_array = [0.0],
)
```

### Fitting NL Models
```julia

fit_NL_model_file(
    label_exp::String, #label of the experiment
    path_to_data::String, # path to the folder to analyze
    model::Any, # string of the used model
    lb_param::Vector{Float64},# array of the array of the lower bound of the parameters
    ub_param::Vector{Float64}; # array of the array of the upper bound of the parameters
    path_to_annotation::Any = missing,# path to the annotation of the wells
    u0=lb_param .+ (ub_param .- lb_param) ./ 2,# initial guess param
    method_of_fitting="MCMC",
    nrep=100,
    errors_estimation=false,
    optmizator=BBO_adaptive_de_rand_1_bin_radiuslimited(), # selection of optimization method 
    path_to_results="NA", # path where save results
    path_to_plot="NA", # path where to save Plots
    loss_type="RE", # string of the type of the used loss
    smoothing=false, # 1 do smoothing of data with rolling average
    type_of_smoothing="lowess",
    display_plots=true,# display plots in julia or not
    save_plots=false,
    verbose=false, # 1 true verbose
    write_res=false, # write results
    pt_avg=1, # number of points to do smoothing average
    pt_smooth_derivative=7, # number of points to do ssmooth_derivative
    do_blank_subtraction="avg_blank", # string on how to use blank (NO,avg_subtraction,time_avg)
    avg_replicate=false, # if true the average between replicates is fitted. If false all replicate are fitted indipendelitly
    correct_negative="thr_correction", # if "thr_correction" it put a thr on the minimum value of the data with blank subracted, if "blank_correction" uses blank distrib to impute negative values
    thr_negative=0.01,  # used only if correct_negative == "thr_correction"
    multiple_scattering_correction=false, # if true uses the given calibration curve to fix the data
    method_multiple_scattering_correction="interpolation",
    calibration_OD_curve="NA",  #  the path to calibration curve to fix the data
    PopulationSize=300,
    maxiters=2000000,
    abstol=0.00001,
    thr_lowess=0.05,
    penality_CI=8.0,
    size_bootstrap=0.7,
    blank_value = 0.0,
    blank_array = [0.0],
)
```

```julia
 fit_NL_model_selection_file(
    label_exp::String, #label of the experiment
    path_to_data::String, # path to the folder to analyze
    list_model_function::Any, # ode model to use
    list_lb_param::Vector{Float64}, # lower bound param
    list_ub_param::Vector{Float64}; # upper bound param
    path_to_annotation::Any = missing,# path to the annotation of the wells
    method_of_fitting="MCMC",
    nrep=100,
    list_u0=lb_param .+ (ub_param .- lb_param) ./ 2,# initial guess param
    optmizator=BBO_adaptive_de_rand_1_bin_radiuslimited(), # selection of optimization method 
    path_to_results="NA", # path where save results
    path_to_plot="NA", # path where to save Plots
    loss_type="RE", # string of the type of the used loss
    smoothing=false, # 1 do smoothing of data with rolling average
    type_of_smoothing="lowess",
    display_plots=true,# display plots in julia or not
    save_plots=false,
    verbose=false, # 1 true verbose
    write_res=false, # write results
    pt_avg=1, # number of points to do smoothing average
    pt_smooth_derivative=7, # number of points to do ssmooth_derivative
    do_blank_subtraction="avg_blank", # string on how to use blank (NO,avg_subtraction,time_avg)
    avg_replicate=false, # if true the average between replicates is fitted. If false all replicate are fitted indipendelitly
    correct_negative="thr_correction", # if "thr_correction" it put a thr on the minimum value of the data with blank subracted, if "blank_correction" uses blank distrib to impute negative values
    thr_negative=0.01,  # used only if correct_negative == "thr_correction"
    multiple_scattering_correction=false, # if true uses the given calibration curve to fix the data
    method_multiple_scattering_correction="interpolation",
    calibration_OD_curve="NA",  #  the path to calibration curve to fix the data
    PopulationSize=300,
    maxiters=2000000,
    abstol=0.00001,
    thr_lowess=0.05,
    beta_param=2.0,
    penality_CI=8.0,
    size_bootstrap=0.7,
    correction_AIC=true,
    blank_value = 0.0,
    blank_array = [0.0],
)
```

### ODE segmentation
```julia

 segmentation_ODE_file(
    label_exp::String, #label of the experiment
    path_to_data::String, # path to the folder to analyze
    list_of_models::Vector{String}, # ode model to use 
    lb_param_array::Any, # lower bound param
    ub_param_array::Any,# upper bound param
    n_max_change_points::Int;
    path_to_annotation::Any = missing,# path to the annotation of the wells
    detect_number_cpd=true,
    fixed_cpd=false,
    optmizator=BBO_adaptive_de_rand_1_bin_radiuslimited(), # selection of optimization method 
    integrator=Tsit5(), # selection of sciml integrator
    type_of_loss="L2", # type of used loss 
    type_of_detection="sliding_win",
    type_of_curve="original",
    do_blank_subtraction="avg_blank",
    correct_negative="thr_correction",
    thr_negative=0.01,
    pt_avg=1, # number of the point to generate intial condition
    smoothing=true, # the smoothing is done or not?
    save_plots=false, # do plots or no
    display_plots=false, # do plots or no
    path_to_plot="NA", # where save plots
    path_to_results="NA",
    win_size=7, # numebr of the point to generate intial condition
    pt_smooth_derivative=0,
    penality_parameter=2.0,
    avg_replicate=false,
    multiple_scattering_correction="false", # if true uses the given calibration curve to fix the data
    method_multiple_scattering_correction="interpolation",
    calibration_OD_curve="NA",  #  the path to calibration curve to fix the data
    write_res=false,
    save_all_model=false,
    method_peaks_detection="peaks_prominence",
    n_bins=40,
    PopulationSize=300,
    maxiters=2000000,
    abstol=0.00001,
    type_of_smoothing="lowess",
    thr_lowess=0.05,
    verbose=false,
    correction_AIC=true,
    blank_value = 0.0,
    blank_array = [0.0],)
```

### NL segmentation
```julia

fit_NL_segmentation_file(
    label_exp::String, #label of the experiment
    path_to_data::String, # path to the folder to analyze
    list_model_function::Any, # ode model to use
    list_lb_param::Vector{Vector{Float64}}, # lower bound param
    list_ub_param::Vector{Vector{Float64}}, # upper bound param
    n_change_points::Int;
    path_to_annotation::Any = missing,# path to the annotation of the wells
    method_of_fitting="MCMC",
    nrep=100,
    list_u0=lb_param .+ (ub_param .- lb_param) ./ 2,# initial guess param
    optmizator=BBO_adaptive_de_rand_1_bin_radiuslimited(), # selection of optimization method 
    path_to_results="NA", # path where save results
    path_to_plot="NA", # path where to save Plots
    loss_type="RE", # string of the type of the used loss
    smoothing=false, # 1 do smoothing of data with rolling average
    type_of_smoothing="lowess",
    display_plots=true,# display plots in julia or not
    save_plots=false,
    verbose=false, # 1 true verbose
    write_res=false, # write results
    pt_avg=1, # number of points to do smoothing average
    pt_smooth_derivative=7, # number of points to do ssmooth_derivative
    do_blank_subtraction="avg_blank", # string on how to use blank (NO,avg_subtraction,time_avg)
    avg_replicate=false, # if true the average between replicates is fitted. If false all replicate are fitted indipendelitly
    correct_negative="thr_correction", # if "thr_correction" it put a thr on the minimum value of the data with blank subracted, if "blank_correction" uses blank distrib to impute negative values
    thr_negative=0.01,  # used only if correct_negative == "thr_correction"
    multiple_scattering_correction=false, # if true uses the given calibration curve to fix the data
    method_multiple_scattering_correction="interpolation",
    calibration_OD_curve="NA",  #  the path to calibration curve to fix the data
    PopulationSize=300,
    maxiters=2000000,
    abstol=0.00001,
    size_bootstrap=0.7,
    thr_lowess=0.05,
    detect_number_cpd=true,
    type_of_detection="sliding_win",
    type_of_curve="original",
    fixed_cpd=false,
    penality_CI=8.0,
    beta_smoothing_ms=2.0,
    win_size=7, # number of the point of cpd sliding win
    n_bins=40,
    correction_AIC=true,
    blank_value = 0.0,
    blank_array = [0.0],
)
```