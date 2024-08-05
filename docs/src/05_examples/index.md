# [Examples and Tutorial](@id examples)

This section provides some copy-and-paste examples of KinBiont.jl

```@contents
Pages = ["index.md"]
Depth = 3
```

## Simulating/Loading single kinetics data
### Loading  data from a .csv
It is possible to load single curves using CSV package in Julia and format them to be compatible with KinBiont. In this example we suppose that the .csv file is formatted as required in KinBiont documentation.

```julia
df_data  =CSV.file("your_path/data.csv")
names_of_cols = propertynames(df_data)  
# assuming first column is time and we want to fit the second one 

data_OD = Matrix(hcat(df_data[names_of_cols[1]],df_data[names_of_cols[2]]))
```
### Simulating Data with ODEs

To simulate data using Ordinary Differential Equations (ODEs):

```julia
# Simulating data with an ODE
model = "triple_piecewise_adjusted_logistic"
n_start = [0.1]
tstart = 0.0
tmax = 500.0
delta_t = 15.0


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
We start applying a rolling average smoothing to the data. In the example, a rolling window of size 7 is applied to the original data (```data_OD``` generated in the previous examples). 
```julia
data_ODsmooth = smoothing_data(data_OD, 7)
data_ODsmooth = Matrix(data_ODsmooth)

# Plotting scatterplot of smoothed data
Plots.scatter(data_ODsmooth[1, :], data_ODsmooth[2, :], xlabel="Time", ylabel="Arb. Units", label=["Smoothed data " nothing], markersize=2, color=:blue, size=(300, 300))
```
Furthermore, to address potential external influences, a correction for multiple scattering is applied to the smoothed data. This correction is executed through the ```correction_OD_multiple_scattering``` function, requiring an external file (```calibration_curve.csv```).  it is  optional in the provided example. 
```julia

# Multiple scattering correction (optional, comment out if not needed)
data_ODsmooth = correction_OD_multiple_scattering(data_ODsmooth, "/your_path/calibration_curve.csv")
using Plots
# Plotting scatterplot of preprocessed data
Plots.scatter(data_ODsmooth[1, :], data_ODsmooth[2, :], xlabel="Time", ylabel="Arb. Units", label=["Pre-processed data" nothing], markersize=2, color=:blue, size=(300, 300))
```

## Fitting a single kinetics
When data are store in a Matrix of Float64 with 2 rows it is possible to perfor various analyses
### Evaluation of the dynamics of specific growth rate
The user can evaluate   of the specific growth rate  during all the curve. This can be done by running the following:


```julia

pt_win = 7
specific_gr_array = specific_gr_evaluation(data_OD,pt_win )
specific_gr_times = [
    (data_OD[1, r] + data_OD[1, 	(r+pt_smoopt_winthing_derivative)]) / 2 for
    r = 1:1:(eachindex(data_OD[2, :])[end].- pt_win)
 	]
Plots.scatter(specific_gr_times,specific_gr_array, xlabel="Time", ylabel="Arb. Units", label=["Data " nothing],  color=:blue, size = (300,300))


```
### Log-Lin fitting

To perform the log-linear fitting  (of data generated in the previous examples) it is suffucient to run. 

```julia
res_log_lin = fitting_one_well_Log_Lin(
    data_OD, # dataset first row times second row OD
   "test", # name of the well
    "test log-lin fitting"; #label of the experiment
    type_of_smoothing="rolling_avg", # type of smoothing
    pt_avg=7, # number of the point for rolling avg not used in the other cases
    pt_smoothing_derivative=7, # number of poits to smooth the derivative
    pt_min_size_of_win=7, # minimum size of the exp windows in number of smooted points
)
```
 The results are stored in the ```res_log_lin```  data struct. 


###    Fitting ODE Models

Before fitting, the model, the initial guess of parameters for the optimization and upper and lower bounds for the ODE parameters are defined 
```julia


model ="aHPM"

P_GUESS = [0.01, 0.001, 1.00, 1]
ub_ahpm = P_GUESS.*4
lb_ahpm = P_GUESS./4
```
The actual fitting is accomplished through the ```fitting_one_well_ODE_constrained``` function.  In this case wed do not use the lb and ub they will be generated automatically.
```julia
# Performing ODE fitting
 results_ODE_fit = fitting_one_well_ODE_constrained(
    data_OD, 
    "test",
    "test_ODE",
    model,
    P_GUESS;
)
Plots.scatter(data_OD[1,:],data_OD[2,:], xlabel="Time", ylabel="Arb. Units", label=["Data " nothing],color=:red,markersize =2 ,size = (300,300))
Plots.plot!(results_ODE_fit[4],results_ODE_fit[3], xlabel="Time", ylabel="Arb. Units", label=["fit " nothing],color=:red,markersize =2 ,size = (300,300))


```
The results are stored in  the second entry of ```results_ODE_fit``` with the following format. 

To add lower and upper bounds one should use the ```the opt_params...```,  in addition we put also ```multistart = true```:
```julia
# Performing ODE fitting

@time results_ODE_fit = fitting_one_well_ODE_constrained(
    data_OD, 
    "test",
    "test_ODE",
    model,
    P_GUESS;
   lb = lb_ahpm,
   ub = ub_ahpm
)

Plots.scatter(data_OD[1,:],data_OD[2,:], xlabel="Time", ylabel="Arb. Units", label=["Data " nothing],color=:red,markersize =2 ,size = (300,300))
Plots.plot!(results_ODE_fit[4],results_ODE_fit[3], xlabel="Time", ylabel="Arb. Units", label=["fit " nothing],color=:red,markersize =2 ,size = (300,300))


```

To change the optimization method one should call the desired library from Optimization.jl and use the keyword argument ```optimizer``` and if required by the selected optimizer it necessary to specify the differentiation methods, e.g.:


- using Broyden–Fletcher–Goldfarb–Shanno algorithm
```julia
using Optimization

@time results_ODE_fit = fitting_one_well_ODE_constrained(
    data_OD, 
    "test",
    "test_ODE",
    model,
    P_GUESS;
   lb = lb_ahpm,
   ub = ub_ahpm,
    multistart = true,
   optimizer=BFGS(), 
)
```

- using PRAXIS algorithm, with Tik-Tak restart (from Benchmarking global optimizers, Arnoud et al 2019)

```julia
using OptimizationNLopt

@time results_ODE_fit = fitting_one_well_ODE_constrained(
    data_OD, 
    "test",
    "test_ODE",
    model,
    P_GUESS;
   lb = lb_ahpm,
   ub = ub_ahpm,
   optimizer= BFGS(), 
)
```

- Changing the number of iterations and absolute tolerance:

```julia

@time results_ODE_fit = fitting_one_well_ODE_constrained(
    data_OD, 
    "test",
    "test_ODE",
    model,
    P_GUESS;
   lb = lb_ahpm,
   ub = ub_ahpm,
   optimizer= BFGS(), 
   abstol = 0.0001,
   maxiters= 20000,
)
```
###   Custom ODE Fitting

Using the custom ODE fitting, users can fit any ordinary differential equation. First it is necessary to define the differential equation. This is done by declaring a new function, e.g.:
```julia

# Custom ODE function
function ODE_custom(du, u, param, t)
    du[1] = u[1] * (1 - u[1]) * param[2] + param[1] * u[1]
end

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
@time results_ODE_fit = fitting_one_well_custom_ODE(
    data_OD, # dataset first row times second row OD
   "test", # name of the well
    "test custom ode", #label of the experiment
   ODE_custom, # ode model to use
   param_guess,
    1; # number ode in the system
  )


```
The results are stored in ```results_ODE_fit``` .

###   ODE Sensitivity Analysis

The sensitivity analysis is performed with the ```one_well_morris_sensitivity``` function. This function takes the preprocessed dataset (```data_OD``` generated in the previous examples), the name and label of the well, the ODE model to use ("aHPM" in this case), as well as the lower and upper bounds for the ODE parameters. The number of steps in the Morris method (```n_step_sensitivity```) should be specified.

```julia
# Number of steps for Morris sensitivity analysis
n_step_sensitivity = 5
P_GUESS = [0.01, 0.001, 1.00, 1]
ub_ahpm = P_GUESS.*5
lb_ahpm = P_GUESS./5
# Performing Morris sensitivity analysis
@time sensitivity_test = one_well_morris_sensitivity(
    data_OD, 
    "test",
     "test_sensitivity",
      "aHPM", 
      lb_ahpm,
       ub_ahpm,
    N_step_morris=n_step_sensitivity ;
)


```
In this case the user will obtain the results of the fit for each of the Morris step and a matrix containg the coordinate of the hyperspace of the starting guess assicioted to each fit.
###   ODE Model Selection

The user can use the model selection function to use AIC (or AICc) to identify the best model to fit a specific kinetics.
We start defining the list of the models and the upper and lower bounds of the parameters of each one:

```julia
# Models candidates and their parameter bounds
list_of_models = ["aHPM",   "baranyi_roberts"]
ub_1 =[ 0.1 , 0.1 , 0.1 , 5.001    ]
lb_1 =[ 0.0001 , 0.001,0.0  , 0.001 ]
p1_guess = lb_1 .+ (ub_1.-lb_1)./2

ub_2 =[ 1.2 , 5.1 , 500.0  ,5.0,5.0  ]
lb_2 =[ 0.0001 , 0.1, 0.00 ,0.2,0.2 ]
p2_guess = lb_2 .+ (ub_2.-lb_2)./2


```


The model selection process is runned with the ```ODE_Model_selection``` function. 

```julia
results_ms = ODE_Model_selection(
    data_OD,
    "test", 
    "test_model_selection",
    list_of_models,
    list_guess;
    multistart = true,
    lb_param_array = list_lb,
    ub_param_array = list_ub  
)

Plots.scatter(data_OD[1,:],data_OD[2,:], xlabel="Time", ylabel="Arb. Units", label=["Data " nothing],color=:red,markersize =2 ,size = (300,300))
Plots.plot!(results_ms[4],results_ms[3], xlabel="Time", ylabel="Arb. Units", label=["fit " nothing],color=:red,markersize =2 ,size = (300,300))


```



The results of the model selection process are stored in the ```results_ms``` variable.

### ODE segmentation
For a single kinetics it is possible to run the segmentaion in two different ways. A manual selection of the change points or a using a change points detection algorithm to find them.

First of all we generate a synthetic daset composed by more than one model  (skip this part if you import a real dataset in the data variable).


```julia
#First segment ODE

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


Plots.scatter(times_sim,sol_sim, xlabel="Time", ylabel="Arb. Units", label=["Data " nothing],  color=:blue, size = (300,300))



data_OD = Matrix(transpose(hcat(times_sim,sol_sim)))

# Adding uniform noise to the dataset
noise_uniform = rand(Uniform(-0.01, 0.01), length(sol_sim))
data_OD = Matrix(transpose(hcat(times_sim, sol_sim)))
data_OD[2, :] = data_OD[2, :] .+ noise_uniform


```
Then we generate the list of models that will be used and their parameters's lower/upper bounds:

```julia
# Initializing all the models for selection
# Initializing all the models for selection
ub_exp = [0.1]
lb_exp = [-0.01]
p1_guess = lb_exp .+(ub_exp.-lb_exp)/.2

ub_logistic = [0.9, 5.0]
lb_logistic = [0.0001, 0.001]
p2_guess = lb_logistic .+(ub_logistic.-lb_logistic)/.2

ub_hpm = [0.1, 20.0, 50.001]
lb_hpm = [0.0001, 0.000001, 0.001]
p3_guess = lb_hpm .+(ub_hpm.-lb_hpm)/.2


list_of_models = ["exponential",  "logistic","HPM"]
list_ub_param = [ub_exp,ub_logistic, ub_hpm]
list_lb_param = [lb_exp, lb_logistic,lb_hpm]
list_guess = [p1_guess, p2_guess, p3_guess]

```

First, we fit giving to KinBiont the list of change points:

```julia
@time seg_fitting = selection_ODE_fixed_intervals(
   data_OD, # dataset first row times second row OD
    "test", # name of the well
    "", #label of the experiment
    list_of_models, # ode models to use
    list_guess,
    cdp_list;
)
Plots.scatter(data_OD[1, :], data_OD[2, :], xlabel="Time", ylabel="Arb. Units", label=["Data " nothing], color=:blue, markersize=2, size=(300, 300))
Plots.plot!(seg_fitting[3], seg_fitting[4], xlabel="Time", ylabel="Arb. Units", label=["fit " nothing], color=:red, markersize=2, size=(300, 300))

```
In alternative, if the change points are not known we can run a cpd algorithm and perfom the fitting with :


```julia
n_change_points =2
@time seg_fitting = segmentation_ODE(
   data_OD, # dataset first row times second row OD
    "test", # name of the well
    "", #label of the experiment
    list_of_models, # ode models to use
    list_guess,
    3;
    detect_number_cpd=true,
    fixed_cpd=false,
)

```
### Fitting NL Models
With KinBiont it is possible to fit any non-linear model this can be done by calling the function ```NL_model_selection``` in different ways.


First we declare upper and lower bound and the model (note that in this case we use array of array because the input can be more than one model)

```julia
nl_model = ["NL_Richards"]
p_guess = [[1.0,1.0,0.01,300.0]]
lb_nl =[[0.01,0.01,0.000001,00.01]]
ub_nl =p_guess.*3
```


To perform a single fit on a time series (i.e., ```data_OD```) then we run the following:
```julia
@time nl_fit = NL_model_selection(data_OD, # dataset first row times second row OD
"test", 
"test_model_selection",
nl_model, #  model to use
p_guess;
)


```


The user can specify any parameter of the optimizer, for the bound in this case it is done via:
```julia
@time nl_fit =  NL_model_selection(data_OD, # dataset first row times second row OD
    "test", 
    "test_model_selection",
    nl_model, #  model to use
    p_guess;
    lb_param_array =lb_nl,
    ub_param_array = ub_nl,
    maxiters=2000000,
    abstol=0.00001,
)

```

Alternately, the user can  opt to a Bootstrap approach (```method_of_fitting ="Bootstrap"```): 



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

As for the ODE model, also for the NL fit it is possible to perform a sensitivity analysis with respect to the initial starting guess of the parameters, in this case we just use ```method_of_fitting= "Morris_sensitivity"``` and nrep as the number of Morris steps


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
and the we perform the fit, tuning the AIC parameter (```beta_param```) :
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
Note that in the case of the NL fit, the fit of the segment goes backward and since it is not possible to specificy the bounduary condition of the fit the user can force the optimization problem to mantain the continuty between segment with the parameter  ```penality_CI=8.0```. 

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

Instead fitting a single kinetics the user can supply a ".csv" file (formatted as described in the section), and KinBiont will proceed to perform all the analysis on all the wells of the experiment. Note that the user can supply a annotation .csv in this case becomes possible to subtract the blanks and fit the average of replicates.

To use the following functions the user should input to KinBiont variables that contains the string of the path to the .csv files:

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

If the paths are provided to ```fit_one_file_Log_Lin```, the user will obtain a matrix containing the results for each well:
```julia
fit_one_file_Log_Lin(
   "test", #label of the experiment
    path_to_data; # path to the folder to analyze
    path_to_annotation = path_to_annotation,# path to the annotation of the wells
    display_plots=true,# display plots in julia or not
    do_blank_subtraction="avg_blank", # string on how to use blank (NO,avg_subtraction,time_avg)
    avg_replicate=true, # if true the average between replicates is fitted. If false all replicate are fitted indipendelitly
    correct_negative="remove", # if "thr_correction" it put a thr on the minimum value of the data with blank subracted, if "blank_correction" uses blank distrib to impute negative values
    blank_value = 0.0,
    blank_array = [0.0],)
```


Note that if the user wants to subtract blank but their are non in the file with the data the to optional arguments can be used    ```blank_value``` for one value of blanks,
  ```blank_array``` for an array of  blanks


### ODE fitting
We declare the model and the bounds of it:

```julia


model ="aHPM"
# Upper bounds of the parameters of the ODE
ub_ahpm = [1.2, 1.1, 2.0, 20]

# Lower bounds of the parameters of the ODE
lb_ahpm = [0.0001, 0.00000001, 0.00, 0]
```

We proceed fitting with the ```fit_file_ODE```  function:

```julia

 fit_file_ODE(
    "test", #label of the experiment
    path_to_data, # path to the folder to analyze
    model, # string of the used model
    lb_ahpm,# array of the array of the lower bound of the parameters
    ub_ahpm; # array of the array of the upper bound of the parameters
    path_to_annotation = path_to_annotation,# path to the annotation of the wells
    optmizator=BBO_adaptive_de_rand_1_bin_radiuslimited(), # selection of optimization method 
    integrator=Tsit5(), # selection of sciml integrator
    display_plots=true,# display plots in julia or not
    do_blank_subtraction="avg_blank", # string on how to use blank (NO,avg_subtraction,time_avg)
    correct_negative="removr", # if "thr_correction" it put a thr on the minimum value of the data with blank subracted, if "blank_correction" uses blank distrib to impute negative values
    PopulationSize=300,
    maxiters=2000000,
    abstol=0.00001,
)
```
The user can fit a custom model by declaring the function in advance

```julia
# Custom ODE function
function ODE_custom(du, u, param, t)
    du[1] = u[1] * (1 - u[1]) * param[2] + param[1] * u[1]
end


# Bounds for the custom ODE parameters
custom_ub = [1.2, 1.1]
custom_lb = [0.0001, 0.00000001]
```

Then, we can call the function ```fit_file_custom_ODE```:

```julia

 fit_file_custom_ODE(
    "test", #label of the experiment
    path_to_data, # path to the folder to analyze
    ODE_custom, #  used model
    custom_lb,# array of the array of the lower bound of the parameters
    custom_ub, # array of the array of the upper bound of the parameters
    1;
   path_to_annotation = path_to_annotation,# path to the annotation of the wells
    optmizator=BBO_adaptive_de_rand_1_bin_radiuslimited(), # selection of optimization method 
    integrator=Tsit5(), # selection of sciml integrator
    display_plots=true,# display plots in julia or not
    do_blank_subtraction="avg_blank", # string on how to use blank (NO,avg_subtraction,time_avg)
    correct_negative="removr", # if "thr_correction" it put a thr on the minimum value of the data with blank subracted, if "blank_correction" uses blank distrib to impute negative values
    PopulationSize=300,
    maxiters=2000000,
    abstol=0.00001,
)
```
It is also possible to perform the model selection on an entire file. First we declare the lis of models

```julia
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
After, we run the function ```ODE_model_selection_file```:

```julia

 ODE_model_selection_file(
    "test", 
    path_to_data, # path to the folder to analyze
    list_of_models, # ode model to use 
    list_lb_param::Any, # lower bound param
    list_ub_param::Any; # upper bound param
    path_to_annotation = path_to_annotation,# path to the annotation of the wells
    optmizator=BBO_adaptive_de_rand_1_bin_radiuslimited(), # selection of optimization method 
    integrator=Tsit5(), # selection of sciml integrator
    display_plot_best_model=true, # one wants the results of the best fit to be plotted
    save_plot_best_model=false,
    beta_penality=2.0, # penality for AIC evaluation
    correct_negative="remove", # if "thr_correction" it put a thr on the minimum value of the data with blank subracted, if "blank_correction" uses blank distrib to impute negative values
    PopulationSize=300,
    maxiters=2000000,
    abstol=0.00001,
    thr_lowess=0.05,
    correction_AIC=false,
)
```

### Fitting NL Models


The user can fit any NL model by calling the function ```fit_NL_model_file```.
As usual the user should declare the model and the bounds



```julia
nl_lb = [0.0001 , 0.00000001, 0.00,0.0 ]
nl_ub = [2.0001 , 10.00000001, 5.00,5.0]
model_nl = "NL_Richards"

```

Alternatively, as in the single kinetics case, the model can also be a function, e.g.,:

```julia
function model_nl(p, times)
    model = p[1] .* exp.(times .* p[2])
    return model
end


nl_ub =  [2.0001 , 10.00000001]
nl_lb =  [0.0001 , 0.00000001 ]

```

Then, we proceed fitting, note that is possible to call any of the previous ```method_of_fitting```

```julia

fit_NL_model_file(
    "test", #label of the experiment
    path_to_data, # path to the folder to analyze
    model_nl, # string of the used model
    nl_lb,# array of the array of the lower bound of the parameters
    nl_ub; # array of the array of the upper bound of the parameters
    path_to_annotation = path_to_annotation,# path to the annotation of the wells
    method_of_fitting="MCMC",
    nrep=50,
    optmizator=BBO_adaptive_de_rand_1_bin_radiuslimited(), # selection of optimization method 
    display_plots=true,# display plots in julia or not
    correct_negative="remove", # if "thr_correction" it put a thr on the minimum value of the data with blank subracted, if "blank_correction" uses blank distrib to impute negative values
    PopulationSize=300,
    maxiters=2000000,
    abstol=0.00001,
)
```
The user can perform a NL model selection. To do that we should start declaring the list of models and upper/lower bounds:


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
and the we can call the NL model selection function:

```julia
 fit_NL_model_selection_file(
    "test", #label of the experiment
    path_to_data, # path to the folder to analyze
    list_models_f, # ode model to use
    list_lb, # lower bound param
   list_ub; # upper bound param
    path_to_annotation = path_to_annotation,# path to the annotation of the wells
    method_of_fitting="MCMC",
    nrep=100,
    optmizator=BBO_adaptive_de_rand_1_bin_radiuslimited(), # selection of optimization method 
    display_plots=true,# display plots in julia or not
    do_blank_subtraction="avg_blank", # string on how to use blank (NO,avg_subtraction,time_avg)
    correct_negative="remove", # if "thr_correction" it put a thr on the minimum value of the data with blank subracted, if "blank_correction" uses blank distrib to impute negative values
    thr_negative=0.01,  # used only if correct_negative == "thr_correction"
    PopulationSize=300,
    maxiters=2000000,
    abstol=0.00001,
    thr_lowess=0.05,
    beta_param=2.0,
    correction_AIC=false,
)
```

### ODE segmentation

It is possible to apply the ODE segmentation to the an entire .csv file. Note that in this case all the change points will be detectec using a off-line change point algorithm.

First we declare the models and upper/lower bounds:
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

Finally, we perform the fit:
```julia

 segmentation_ODE_file(
    "test seg ODE",
    path_to_data, # path to the folder to analyze
    list_of_models, # ode model to use 
    lb_param_array, # lower bound param
    ub_param_array,# upper bound param
    n_max_change_point;
    path_to_annotation = path_to_annotation,# path to the annotation of the wells
    detect_number_cpd=false,
    fixed_cpd=true,
    optmizator=BBO_adaptive_de_rand_1_bin_radiuslimited(), # selection of optimization method 
    integrator=Tsit5(), # selection of sciml integrator
    correct_negative="remove",
    thr_negative=0.01,
    smoothing=false, # the smoothing is done or not?
    display_plots=true, # do plots or no
    win_size=8, 
    PopulationSize=300,
    maxiters=2000000,
    abstol=0.00001,
    pt_avgh = 3,
    type_of_smoothing="rolling_avg",
    )
```

### NL segmentation

It is possible to apply the NL segmentation to the an entire .csv file. Note that in this case all the change points will be detectec using a off-line change point algorithm.

We start declaring models and upper/lower bounds:

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

n_change_points =2
```

Then, we call the function to perform the fit:
```julia

fit_NL_segmentation_file(
    "test", #label of the experiment
    path_to_data, # path to the folder to analyze
    list_models_f, # ode model to use
    list_lb, # lower bound param
    list_ub, # upper bound param
    n_change_points;
    detect_number_cpd=false,
    fixed_cpd=true,
    path_to_annotation = path_to_annotation,# path to the annotation of the wells
    method_of_fitting="MCMC",
    nrep=50,
    optmizator=BBO_adaptive_de_rand_1_bin_radiuslimited(), # selection of optimization method 
    loss_type="RE", # string of the type of the used loss
    display_plots=true,# display plots in julia or not
    do_blank_subtraction="avg_blank", # string on how to use blank (NO,avg_subtraction,time_avg)
    correct_negative="remove", # if "thr_correction" it put a thr on the minimum value of the data with blank subracted, if "blank_correction" uses blank distrib to impute negative values
    PopulationSize=200,
    maxiters=2000000,
    abstol=0.00001,
    penality_CI=8.0,
    beta_smoothing_ms=2.0,
    win_size=7, # number of the point of cpd sliding win
    correction_AIC=false,
)
```

## Downstream ML anlyses

# Symbolic regression
# Decision tree/forest regression


