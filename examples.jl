using Kimchi
using DifferentialEquations
using CSV
using Distributions
using Plots
using StatsBase
using Optimization
using OptimizationOptimJL
using OptimizationMOI
using OptimizationNLopt
using OptimizationCMAEvolutionStrategy
using OptimizationBBO
using OptimizationEvolutionary
using OptimizationGCMAES
using OptimizationOptimisers
using OptimizationMultistartOptimization
using OptimizationPRIMA
using OptimizationPolyalgorithms

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

Plots.scatter(data_OD[1,:],data_OD[2,:], xlabel="Time", ylabel="Arb. Units", label=["Data " nothing],color=:red,markersize =2 ,size = (300,300))


res_log_lin = fitting_one_well_Log_Lin(
    data_OD, # dataset first row times second row OD
   "test", # name of the well
    "test log-lin fitting"; #label of the experiment
    type_of_smoothing="rolling_avg", # type of smoothing
    pt_avg=7, # number of the point for rolling avg not used in the other cases
    pt_smoothing_derivative=7, # number of poits to smooth the derivative
    pt_min_size_of_win=7, # minimum size of the exp windows in number of smooted points
)


# Fitting with ODE


model ="aHPM"
# Upper bounds of the parameters of the ODE
ub_ahpm = [1.2, 1.1, 2.0, 20]

# Lower bounds of the parameters of the ODE
lb_ahpm = [0.0001, 0.00000001, 0.00, 0]

# Performing ODE fitting
results_ODE_fit = fitting_one_well_ODE_constrained(
    data_OD, 
    "test",
    "test_ODE",
    model,
    lb_ahpm,
    ub_ahpm;
    optimizer =  NLopt.LN_PRAXIS(),
    smoothing=true, # the smoothing is done or not?
    pt_avg=3, # number of the points to do rolling avg
   lb = lb_ahpm,
   ub = ub_ahpm

)
Plots.plot!(results_ODE_fit[4],results_ODE_fit[3], xlabel="Time", ylabel="Arb. Units", label=["fit " nothing],color=:red,markersize =2 ,size = (300,300))

# Custom ODE function
function ODE_custom(du, u, param, t)
    du[1] = u[1] * (1 - u[1]) * param[2] + param[1] * u[1]
end
custom_ub = [1.2, 1.1]
custom_lb = [0.0001, 0.00000001]



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
    pt_avg=3, # numebr of the point to generate intial condition
    smoothing=true, # the smoothing is done or not?
    PopulationSize=25,
    maxiters=100000,
)



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
    PopulationSize=50,
    maxiters=100000,
    abstol=0.00001,
)


# Models candidates and their parameter bounds
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
    PopulationSize=15,
    maxiters=100000,
    abstol=0.00001,
    correction_AIC=false,
)

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


Plots.scatter(sol_sim, xlabel="Time", ylabel="Arb. Units", label=["Data " nothing],  color=:blue, size = (300,300))

Plots.scatter(sol_sim[40:60], xlabel="Time", ylabel="Arb. Units", label=["Data " nothing],  color=:blue, size = (300,300))


data_OD = Matrix(transpose(hcat(times_sim,sol_sim)))

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

# fitting fixed cp
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
    path_to_plot="NA", # where save plots
    pt_smooth_derivative=0,
    PopulationSize=25,
    maxiters=100000,
    correction_AIC=false)



# fitting cp selectet with algortimh
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
    win_size=10, #  
    pt_smooth_derivative=0,
    PopulationSize=25,
    maxiters=100000,
    abstol=0.00001,
    type_of_smoothing="rolling_average",
    thr_lowess=0.05,
    correction_AIC=true)