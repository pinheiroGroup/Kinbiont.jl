using Kimchi
using DifferentialEquations
using CSV
using Distributions
using StatsBase
using OptimizationBBO
using Optimization
using OptimizationOptimJL
using Plots

# Now i generate data with change point


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


Plots.scatter(times_sim,sol_sim, xlabel="Time", ylabel="Arb. Units", label=["Data " nothing],  color=:blue, size = (300,300))



data_OD = Matrix(transpose(hcat(times_sim,sol_sim)))

# Adding uniform noise to the dataset
noise_uniform = rand(Uniform(-0.01, 0.01), length(sol_sim))
data_OD = Matrix(transpose(hcat(times_sim, sol_sim)))
data_OD[2, :] = data_OD[2, :] .+ noise_uniform



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
# fitting fixed cp
cdp_list = [100.0, 200.0]

@time seg_fitting = selection_ODE_fixed_intervals(
   data_OD, # dataset first row times second row OD
    "test", # name of the well
    "", #label of the experiment
    list_of_models, # ode models to use
    list_guess,
    cdp_list;
    beta_smoothing_ms=2.0, #  parameter of the AIC penality
    correction_AIC=true,
)
Plots.scatter(data_OD[1, :], data_OD[2, :], xlabel="Time", ylabel="Arb. Units", label=["Data " nothing], color=:blue, markersize=2, size=(300, 300))
Plots.plot!(seg_fitting[3], seg_fitting[4], xlabel="Time", ylabel="Arb. Units", label=["fit " nothing], color=:red, markersize=2, size=(300, 300))


seg_fitting = selection_ODE_fixed_intervals(
   data_OD, # dataset first row times second row OD
    "test", # name of the well
    "", #label of the experiment
    list_of_models, # ode models to use
    list_guess,
    cdp_list;
    lb_param_array=list_lb_param, # lower bound param
    ub_param_array=list_ub_param, # upper bound param
    beta_smoothing_ms=2.0, #  parameter of the AIC penality
    correction_AIC=true,
    multistart=true,
  
)
Plots.scatter(data_OD[1, :], data_OD[2, :], xlabel="Time", ylabel="Arb. Units", label=["Data " nothing], color=:blue, markersize=2, size=(300, 300))
Plots.plot!(seg_fitting[3], seg_fitting[4], xlabel="Time", ylabel="Arb. Units", label=["fit " nothing], color=:red, markersize=2, size=(300, 300))
# fitting fixed number of cp

@time seg_fitting = segmentation_ODE(
   data_OD, # dataset first row times second row OD
    "test", # name of the well
    "", #label of the experiment
    list_of_models, # ode models to use
    list_guess,
    2;
    detect_number_cpd=true,
    fixed_cpd=false,
    beta_smoothing_ms=2.0, #  parameter of the AIC penality
    correction_AIC=true,
)

Plots.scatter(data_OD[1, :], data_OD[2, :], xlabel="Time", ylabel="Arb. Units", label=["Data " nothing], color=:blue, markersize=2, size=(300, 300))
Plots.plot!(seg_fitting[4], seg_fitting[3], xlabel="Time", ylabel="Arb. Units", label=["fit " nothing], color=:red, markersize=2, size=(300, 300))


@time seg_fitting = segmentation_ODE(
   data_OD, # dataset first row times second row OD
    "test", # name of the well
    "", #label of the experiment
    list_of_models, # ode models to use
    list_guess,
    2;
    detect_number_cpd=false,
    fixed_cpd=false,
    beta_smoothing_ms=2.0, #  parameter of the AIC penality
    correction_AIC=true,
)

Plots.scatter(data_OD[1, :], data_OD[2, :], xlabel="Time", ylabel="Arb. Units", label=["Data " nothing], color=:blue, markersize=2, size=(300, 300))
Plots.plot!(seg_fitting[4], seg_fitting[3], xlabel="Time", ylabel="Arb. Units", label=["fit " nothing], color=:red, markersize=2, size=(300, 300))



@time seg_fitting = segmentation_ODE(
   data_OD, # dataset first row times second row OD
    "test", # name of the well
    "", #label of the experiment
    list_of_models, # ode models to use
    list_guess,
    2;
    detect_number_cpd=false,
    fixed_cpd=true,
    beta_smoothing_ms=2.0, #  parameter of the AIC penality
    correction_AIC=true,
)

Plots.scatter(data_OD[1, :], data_OD[2, :], xlabel="Time", ylabel="Arb. Units", label=["Data " nothing], color=:blue, markersize=2, size=(300, 300))
Plots.plot!(seg_fitting[4], seg_fitting[3], xlabel="Time", ylabel="Arb. Units", label=["fit " nothing], color=:red, markersize=2, size=(300, 300))

# testing restart

@time seg_fitting = segmentation_ODE(
   data_OD, # dataset first row times second row OD
    "test", # name of the well
    "", #label of the experiment
    list_of_models, # ode models to use
    list_guess,
    2;
    detect_number_cpd=false,
    fixed_cpd=true,
    multistart =true,
    beta_smoothing_ms=2.0, #  parameter of the AIC penality
    correction_AIC=true,
    lb_param_array=list_lb_param, # lower bound param
    ub_param_array=list_ub_param, # upper bound param
)

Plots.scatter(data_OD[1, :], data_OD[2, :], xlabel="Time", ylabel="Arb. Units", label=["Data " nothing], color=:blue, markersize=2, size=(300, 300))
Plots.plot!(seg_fitting[4], seg_fitting[3], xlabel="Time", ylabel="Arb. Units", label=["fit " nothing], color=:red, markersize=2, size=(300, 300))



# NL segmentation fitting 







# Initializing all the models for selection
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
# fitting fixed cp
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

Plots.scatter(data_OD[1, :], data_OD[2, :], xlabel="Time", ylabel="Arb. Units", label=["Data " nothing], color=:blue, markersize=2, size=(300, 300))
Plots.plot!(seg_fitting[3], seg_fitting[2], xlabel="Time", ylabel="Arb. Units", label=["fit " nothing], color=:red, markersize=2, size=(300, 300))


seg_fitting = selection_NL_fixed_interval(
   data_OD, # dataset first row times second row OD
    "test", # name of the well
    "", #label of the experiment
    list_of_models_nl, # ode models to use
    list_guess,
    cdp_list;
    lb_param_array=list_lb_param, # lower bound param
    ub_param_array=list_ub_param, # upper bound param
    beta_smoothing_ms=2.0, #  parameter of the AIC penality
    correction_AIC=true,
    multistart=true,
  
)
Plots.scatter(data_OD[1, :], data_OD[2, :], xlabel="Time", ylabel="Arb. Units", label=["Data " nothing], color=:blue, markersize=2, size=(300, 300))
Plots.plot!(seg_fitting[3], seg_fitting[2], xlabel="Time", ylabel="Arb. Units", label=["fit " nothing], color=:red, markersize=2, size=(300, 300))
# fitting fixed number of cp

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

Plots.scatter(data_OD[1, :], data_OD[2, :], xlabel="Time", ylabel="Arb. Units", label=["Data " nothing], color=:blue, markersize=2, size=(300, 300))
Plots.plot!(seg_fitting[4], seg_fitting[3], xlabel="Time", ylabel="Arb. Units", label=["fit " nothing], color=:red, markersize=2, size=(300, 300))


@time seg_fitting = segmentation_NL(
   data_OD, # dataset first row times second row OD
    "test", # name of the well
    "", #label of the experiment
    list_of_models_nl, # ode models to use
    list_guess,
    2;
    detect_number_cpd=false,
    fixed_cpd=false,
    beta_smoothing_ms=2.0, #  parameter of the AIC penality
    correction_AIC=true,
)

Plots.scatter(data_OD[1, :], data_OD[2, :], xlabel="Time", ylabel="Arb. Units", label=["Data " nothing], color=:blue, markersize=2, size=(300, 300))
Plots.plot!(seg_fitting[4], seg_fitting[3], xlabel="Time", ylabel="Arb. Units", label=["fit " nothing], color=:red, markersize=2, size=(300, 300))



@time seg_fitting = segmentation_NL(
   data_OD, # dataset first row times second row OD
    "test", # name of the well
    "", #label of the experiment
    list_of_models_nl, # ode models to use
    list_guess,
    2;
    detect_number_cpd=false,
    fixed_cpd=true,
    beta_smoothing_ms=2.0, #  parameter of the AIC penality
    correction_AIC=true,
)

Plots.scatter(data_OD[1, :], data_OD[2, :], xlabel="Time", ylabel="Arb. Units", label=["Data " nothing], color=:blue, markersize=2, size=(300, 300))
Plots.plot!(seg_fitting[4], seg_fitting[3], xlabel="Time", ylabel="Arb. Units", label=["fit " nothing], color=:red, markersize=2, size=(300, 300))

# testing restart

@time seg_fitting = segmentation_NL(
   data_OD, # dataset first row times second row OD
    "test", # name of the well
    "", #label of the experiment
    list_of_models_nl, # ode models to use
    list_guess,
    2;
    detect_number_cpd=false,
    fixed_cpd=true,
    multistart =true,
    beta_smoothing_ms=2.0, #  parameter of the AIC penality
    correction_AIC=true,
    lb_param_array=list_lb_param, # lower bound param
    ub_param_array=list_ub_param, # upper bound param
)

Plots.scatter(data_OD[1, :], data_OD[2, :], xlabel="Time", ylabel="Arb. Units", label=["Data " nothing], color=:blue, markersize=2, size=(300, 300))
Plots.plot!(seg_fitting[4], seg_fitting[3], xlabel="Time", ylabel="Arb. Units", label=["fit " nothing], color=:red, markersize=2, size=(300, 300))
