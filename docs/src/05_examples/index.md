# [Examples and Tutorial](@id examples)

This section provides some copy-and-paste examples of Kinbiont.jl

```@contents
Pages = ["index.md"]
Depth = 3
```

## Simulating/Loading single kinetics data
### Loading  data from a .csv
It is possible to load single curves using CSV package in Julia and format them to be compatible with Kinbiont. In this example we suppose that the .csv file is formatted as required in Kinbiont documentation.

```julia
df_data  =CSV.file("your_path/data_examples/plate_data.csv")
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
data_ODsmooth = correction_OD_multiple_scattering(data_ODsmooth, "/your_path/data_examples/cal_curve_examples.csv")
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
   optimizer=LN_PRAXIS(), 
       multistart = true,
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

First, we fit giving to Kinbiont the list of change points:

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
With Kinbiont it is possible to fit any non-linear model this can be done by calling the function ```NL_model_selection``` in different ways.


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



### NL Sensitivity Analysis

As for the ODE model, also for the NL fit it is possible to perform a sensitivity analysis with respect to the initial starting guess of the parameters, in this case we just use ```method_of_fitting= "Morris_sensitivity"``` and nrep as the number of Morris steps


```julia
@time nl_fit =  NL_model_selection(data_OD, # dataset first row times second row OD
"test", 
"test_model_selection",
nl_model, #  model to use
p_guess;
nrep =3,
method_of_fitting ="Morris_sensitivity",
multistart = true,
lb_param_array = lb_nl,
ub_param_array = ub_nl
)

Plots.scatter(data_OD[1,:],data_OD[2,:], xlabel="Time", ylabel="Arb. Units", label=["Data " nothing],color=:red,markersize =2 ,size = (300,300))
Plots.plot!(nl_fit[4],nl_fit[3], xlabel="Time", ylabel="Arb. Units", label=["fit " nothing],color=:red,markersize =2 ,size = (300,300))



```



### NL Model Selection
To perform model selection we just specify an array with the list of all the used models:
```julia

nl_ub_1 =  [2.0001 , 10.00000001, 500.00]
nl_lb_1 =  [0.0001 , 0.00000001, 0.00 ]
p1_guess = nl_lb_1 .+ (nl_ub_1.-nl_lb_1)./2

nl_ub_2 =  [2.0001 , 10.00000001, 5.00,5.0]
nl_lb_2 =  [0.0001 , 0.00000001, 0.00,0.0 ]
p2_guess = nl_lb_2 .+ (nl_ub_2.-nl_lb_2)./2

nl_ub_3 =  [2.0001 , 10.00000001]
nl_lb_3 =  [0.0001 , 0.00000001]
p3_guess = nl_lb_3 .+ (nl_ub_3.-nl_lb_3)./2

nl_ub_4 =  [2.0001 , 10.00000001, 500.00]
nl_lb_4 =  [0.0001 , 0.00000001, 0.00 ]
p4_guess = nl_lb_4 .+ (nl_ub_4.-nl_lb_4)./2

list_models_f = ["NL_Gompertz","NL_Bertalanffy","NL_exponential","NL_Gompertz"]
list_lb =[nl_lb_1,nl_lb_2,nl_lb_3,nl_lb_4]
list_ub = [nl_ub_1,nl_ub_2,nl_ub_3,nl_ub_4]
list_guess = [p1_guess,p2_guess,p3_guess,  p4_guess]

```
and the we perform the fit:
```julia

results_ms = ODE_Model_selection(
    data_OD,
    "test", 
    "test_model_selection",
    list_models_f,
    list_guess;
    lb_param_array = list_lb,
    ub_param_array = list_ub  
)

```

### NL segmentation

As the ODE case, for NL segmentation we can perform a segmente fit with manual selection of the change points or a using a change points detection algorithm to find them.
Note that in the case of the NL fit, the fit of the segment goes backward and since it is not possible to specificy the bounduary condition of the fit the user can force the optimization problem to mantain the continuty between segment with the parameter  ```penality_CI=8.0```. 

we specify the list of used models and parameters bounds:

```julia

ub_1 = [0.3 , 0.1]
lb_1 = [0.01 , -0.01]
p1_guess = lb_1 .+(ub_1.-lb_1)/.2

ub_2 = [1.9, 0.1,500.0]
lb_2 = [0.0001,0.001 ,0.001]
p2_guess = lb_2 .+(ub_2.-lb_2)/.2

ub_3 = [0.1, 1.1, 500.0]
lb_3 = [0.0001, 0.000001, 0.001]
p3_guess = lb_3 .+(ub_3.-lb_3)/.2


list_of_models_nl = ["NL_exponential",  "NL_logistic","NL_Gompertz"]
list_ub_param = [ub_1,ub_2, ub_3]
list_lb_param = [lb_1, lb_2,lb_3]
list_guess = [p1_guess, p2_guess, p3_guess]

```


As before we specify the intervals of the segment and we proceed to fit:


```julia
cdp_list = [100.0, 200.0]

@time seg_fitting = selection_NL_fixed_interval(
   data_OD, # dataset first row times second row OD
    "test", # name of the well
    "", #label of the experiment
    list_of_models_nl, # ode models to use
    list_guess,
    cdp_list;
    pt_smooth_derivative = 3
)
```
Another option could be to  run a cpd algorithm and perfom the fitting:




```julia
n_change_points = 2

@time seg_fitting = segmentation_NL(
   data_OD, # dataset first row times second row OD
    "test", # name of the well
    "", #label of the experiment
    list_of_models_nl, # ode models to use
    list_guess,
    2;
    detect_number_cpd=true,
    fixed_cpd=false,
    beta_smoothing_ms=2.0, #  parameter of the AIC penality
    correction_AIC=true,
)


```
## Fitting a .csv file

Instead fitting a single kinetics the user can supply a ".csv" file (formatted as described in the section), and Kinbiont will proceed to perform all the analysis on all the wells of the experiment. Note that the user can supply a annotation .csv in this case becomes possible to subtract the blanks and fit the average of replicates.

To use the following functions the user should input to Kinbiont variables that contains the string of the path to the .csv files:

```julia
path_to_data = "your_path/data_examples/plate_data.csv"
path_to_annotation ="your_path_to_annotation/annotation.csv"
```
In the following examples we assume that these two variables have a value.



### Log-Lin fitting

If the paths are provided to ```fit_one_file_Log_Lin```, the user will obtain a matrix containing the results for each well:

```julia
fit_log_lin = fit_one_file_Log_Lin(
    " ", #label of the experiment
    path_to_data; # path to the folder to analyze
    path_to_annotation=path_to_annotation,# path to the annotation of the wells
  )
```



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


model = "baranyi_richards"

lb_param = [0.001,0.1,0.0,0.01]
ub_param =[0.1,5.0 ,1000.0,5.01]
param_guess =[0.01,1.0 ,500.0,1.01]

fit_od = fit_file_ODE(
    "test", #label of the experiment
    path_to_data, # path to the folder to analyze
    model, # string of the used model
    param_guess;
    path_to_annotation=path_to_annotation,# path to the annotation of the wells
    lb = lb_param, 
    ub =ub_param,
)

```
It is also possible to perform the model selection on an entire file. First we declare the lis of models

```julia
model_1 = "baranyi_richards"

lb_param_1 = [0.001,    0.1 , 0.0   , 0.01]
ub_param_1 = [0.1  ,    5.0 , 1000.0, 5.01]
param_guess_1 =[0.01,1.0 ,500.0,1.01]


model_2 = "aHPM"

lb_param_2 = [0.001,0.0001,0.01,0.01]
ub_param_2 =[0.1,0.1 ,2.0,5.01]
param_guess_2 =[0.01,0.01 ,1.0,1.01]

list_of_models = [model_1,  model_2]
list_ub_param = [ub_param_1,ub_param_2]
list_lb_param = [lb_param_1, lb_param_2]
list_guess = [param_guess_1, param_guess_2]

```
After, we run the function ```ODE_model_selection_file```:

```julia

@time ms_file = ODE_model_selection_file(
    "test", #label of the experiment
    path_to_data, # path to the folder to analyze
    list_of_models, # ode model to use 
    list_guess;
    lb_param_array=list_lb_param, # lower bound param
    ub_param_array=list_ub_param, # upper bound param
    path_to_annotation=path_to_annotation,# path to the annotation of the wells
)
```

### Fitting NL Models


The user can fit any NL model by calling the function ```fit_NL_model_file```.
As usual the user should declare the model and the bounds



```julia
nl_model = ["NL_Richards"]
p_guess = [[1.0,1.0,0.01,300.0]]
lb_nl =[[0.01,0.01,0.000001,00.01]]
ub_nl =p_guess.*50

```
and then fit:
```julia
@time fit_nl = fit_NL_model_selection_file(
    "TEST", #label of the experiment
    path_to_data    , # path to the folder to analyze
    nl_model, # ode model to use
    p_guess;# initial guess param
    lb_param_array =lb_nl, # lower bound param
    ub_param_array = ub_nl, # upper bound param
    path_to_annotation = path_to_annotation,# path to the annotation of the wells
)

```
The user can perform a NL model selection. To do that we should start declaring the list of models and upper/lower bounds:


```julia


nl_model = ["NL_Richards","NL_Bertalanffy"]
p_guess = [[1.0,1.0,0.01,300.0],[0.08,1.0,0.01,1.0]]
lb_nl =[[0.01,0.01,0.000001,00.01],[0.00,0.01,0.001,00.01]]
ub_nl =p_guess.*50

```
and the we can call the NL model selection function:

```julia
@time fit_nl = fit_NL_model_selection_file(
    "TEST", #label of the experiment
    path_to_data    , # path to the folder to analyze
    nl_model, # ode model to use
    p_guess;# initial guess param
    lb_param_array =lb_nl, # lower bound param
    ub_param_array = ub_nl, # upper bound param
    path_to_annotation = path_to_annotation,# path to the annotation of the wells
  
)
```

### ODE segmentation

It is possible to apply the ODE segmentation to the an entire .csv file. Note that in this case all the change points will be detectec using a off-line change point algorithm.

First we declare the models and upper/lower bounds:
```julia
model_1 = "baranyi_richards"

lb_param_1 = [0.001,    0.1 , 0.0   , 0.01]
ub_param_1 = [0.1  ,    5.0 , 1000.0, 5.01]
param_guess_1 =[0.01,1.0 ,500.0,1.01]


model_2 = "aHPM"

lb_param_2 = [0.001,0.0001,0.01,0.01]
ub_param_2 =[0.1,0.1 ,2.0,5.01]
param_guess_2 =[0.01,0.01 ,1.0,1.01]

list_of_models = [model_1,  model_2]
list_ub_param = [ub_param_1,ub_param_2]
list_lb_param = [lb_param_1, lb_param_2]
list_guess = [param_guess_1, param_guess_2]

```

Finally, we perform the fit:
```julia

@time ms_file =  segmentation_ODE_file(
    " ", #label of the experiment
    path_to_data, # path to the folder to analyze
    list_of_models, # ode model to use 
    list_guess, #  param
    1;
    path_to_annotation=path_to_annotation,# path to the annotation of the wells
    detect_number_cpd=false,
    fixed_cpd=true,
    type_of_curve="deriv",
    win_size=7, # numebr of the point to generate intial condition
     maxiters =5,
    lb_param_array=list_lb_param, # lower bound param
    ub_param_array=list_ub_param, # upper bound param
)
```

### NL segmentation

It is possible to apply the NL segmentation to the an entire .csv file. Note that in this case all the change points will be detectec using a off-line change point algorithm.

We start declaring models and upper/lower bounds:

```julia

nl_model = ["NL_Richards","NL_Bertalanffy"]
p_guess = [[1.0,1.0,0.01,300.0],[0.08,1.0,0.01,1.0]]
lb_nl =[[0.01,0.01,0.000001,00.01],[0.00,0.01,0.001,00.01]]
ub_nl =p_guess.*50
n_change_points =2
```

Then, we call the function to perform the fit:
```julia

ms_segmentation = fit_NL_segmentation_file(
    "test", #label of the experiment
    path_to_data, # path to the folder to analyze
    nl_model, # ode model to use
    p_guess,# initial guess param
    1;
    lb_param_array=lb_nl, # lower bound param
    ub_param_array=ub_nl, # upper bound param
    path_to_annotation = path_to_annotation,# path to the annotation of the wells
)
```

## Downstream ML anlyses

### Symbolic regression

This example demonstrates how to use the `Kinbiont` and `SymbolicRegression` packages to analyze kinetics data.

Set up paths to your data, annotation, calibration curve, and result directories (see examples):

```julia
using Kinbiont
using Plots
using Tables
using SymbolicRegression

path_to_data = "your_path/data_examples/plate_data.csv"
path_to_annotation = "your_path/data_examples/annotation.csv"
path_to_calib = "your_path/data_examples/cal_curve_example.csv"
path_to_results = "your_path//seg_res/"
```

We fit with segmentation and 1 change point, we declare the models

```julia
model1 = "HPM_exp"
lb_param1 = [0.00001, 0.000001]
ub_param1 = [0.5, 1.5]
param_guess1 = [0.01, 0.01]

model2 = "logistic"
lb_param2 = [0.00001, 0.000001]
ub_param2 = [0.5, 2.5]
param_guess2 = [0.01, 1.01]

list_of_models = [model1, model2]
list_guess = [param_guess1, param_guess2]
list_lb = [lb_param1, lb_param2]
list_ub = [ub_param1, ub_param2]
```

We perform  the fitting

```julia
fit_file = segmentation_ODE_file(
    "seg_exp_1_test_2", # Label of the experiment
    path_to_data, # Path to the folder to analyze
    list_of_models, # ODE models to use
    list_guess, # Parameter guesses
    1;
    path_to_annotation=path_to_annotation, # Path to the annotation of the wells
    detect_number_cpd=false,
    fixed_cpd=false,
    path_to_calib=path_to_calib,
    multiple_scattering_correction=false, # If true, use calibration curve for data correction
    type_of_curve="deriv",
    pt_smooth_derivative=10,
    type_of_smoothing="lowess",
    verbose=true,
    write_res=false,
    path_to_results=path_to_results,
    win_size=12, # Number of points to generate initial condition
    smoothing=true,
    lb_param_array=list_lb, # Lower bounds for parameters
    ub_param_array=list_ub, # Upper bounds for parameters
    maxiters=200000
)
```

Loading annotation to find concentrations and selecting results of a specific strain

```julia
annotation_test = CSV.File(path_to_annotation, header=false)
names_of_annotation = propertynames(annotation_test)
feature_matrix = hcat(annotation_test[:Column1], annotation_test[:Column3])

# Selecting strain S5
index_strain = findall(annotation_test[:Column2] .== "S5")
index_cc = findall(annotation_test[:Column3] .> 0.01)
to_keep = intersect(index_strain, index_cc)
feature_matrix = feature_matrix[to_keep, :]
wells_to_use = feature_matrix[:, 1]
index_res = Any

```

Giving same order to the well to analyze and the features (can be skipped)
```julia
for i in wells_to_use
    if i == wells_to_use[1]
        index_res = findfirst(res_first_seg[:, 1:end] .== i)
    else
        index_res = hcat(index_res, findfirst(res_first_seg[:, 1:end] .== i))
    end
end

iii = [index_res[1, i][2] for i in eachindex(index_res[1, :])]
res_first_seg_ML = res_first_seg[:, iii]
res_first_seg_ML = hcat(res_first_seg[:, 1], res_first_seg_ML)
```

Add x = 0.0, y = 0.0 to data for not growing cells

```julia

feature_matrix =vcat(feature_matrix,["zero" 0.0])
res_first_seg_ML=hcat(res_first_seg_ML , reduce(vcat,["zero" ,"zero", "zero", 0.0 ,  0.0 ,0.0 ,0.0 ,0.0 ,0.0 ]))


```

Declaring the options of symbolic regression

```julia

options = SymbolicRegression.Options(
    binary_operators=[+, /, *, -],
    unary_operators=[],
    constraints=nothing,
    elementwise_loss=nothing,
    loss_function=nothing,
    tournament_selection_n=12,
    tournament_selection_p=0.86,
    topn=12,
    complexity_of_operators=nothing,
    complexity_of_constants=nothing,
    complexity_of_variables=nothing,
    parsimony=0.05,
    dimensional_constraint_penalty=nothing,
    alpha=0.100000,
    maxsize=10,
    maxdepth=nothing
)
```

Run the symbolic regression using dependent variable that is the 7th row of the Kinbiont results (i.e., the growth rate)

```julia

gr_sy_reg = downstream_symbolic_regression(res_first_seg_ML, feature_matrix, 7; options=options)

scatter(feature_matrix[:, 2], res_first_seg_ML[7, 2:end], xlabel="Amino Acid concentration μM", ylabel="Growth rate [1/Min]", label=["Data" nothing])
hline!(unique(gr_sy_reg[3][:, 1]), label=["Eq. 1" nothing], line=(3, :green, :dash))
plot!(unique(convert.(Float64, feature_matrix[gr_sy_reg[4], 2])), unique(gr_sy_reg[3][:, 2]), label=["Eq. 2" nothing], line=(3, :red))
plot!(unique(convert.(Float64, feature_matrix[gr_sy_reg[4], 2])), unique(gr_sy_reg[3][:, 3]), label=["Eq. 3" nothing], line=(3, :blue, :dashdot))
plot!(unique(convert.(Float64, feature_matrix[gr_sy_reg[4], 2])), unique(gr_sy_reg[3][:, 4]), label=["Eq. 4" nothing], line=(2, :black))
plot!(unique(convert.(Float64, feature_matrix[gr_sy_reg[4], 2])), unique(gr_sy_reg[3][:, 5]), label=["Eq. 5" nothing], line=(2, :black))

```

### Decision Tree Regression Analysis



We declare the packages

```julia
using Kinbiont
using Plots
using CSV, DataFrames
using Statistics
using DelimitedFiles
using Random
using DecisionTree
using AbstractTrees
using MLJDecisionTreeInterface
using TreeRecipe
```


Read the data from CSV files:

```julia
Kinbiont_res_test = readdlm("your_path/data_examples/Results_for_ML.csv", ',')
annotation_test = readdlm("your_path/data_examples/annotation_for_ML.csv", ',')
```


Define necessary variables for analysis:

```julia
ordered_strain = annotation_test[:, end]
n_folds = 10
feature_names = unique(annotation_test[1, 2:end])[2:(end-1)]

depth = -1 


# Set random seed for reproducibility
seed = Random.seed!(1234)
```

 Perform Decision Tree Regression

Loop through each strain and perform decision tree regression and we plot the tree trying to analyze th 9th row (i.e. growht rate) of the results

```julia


index_strain = findall("N. soli".== ordered_strain)
feature_matrix = annotation_test[index_strain, 2:(end-1)]
Kinbiont_results = Kinbiont_res_test[:, index_strain]

dt_gr = downstream_decision_tree_regression(Kinbiont_results,
        feature_matrix,
        9;
        do_pruning=false,
        pruning_accuracy=1.00,
        verbose=true,
        do_cross_validation=true,
        max_depth=depth,
        n_folds_cv=n_folds,
        seed=seed
    )
    
# Wrap the decision tree model for visualization
wt = DecisionTree.wrap(dt_gr[1], (featurenames = feature_names,))

# Plot the decision tree
p2 = Plots.plot(wt, 0.9, 0.2; size = (900, 400), connect_labels = ["yes", "no"])
```
