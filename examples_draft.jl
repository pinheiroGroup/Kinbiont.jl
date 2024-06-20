
using Kimchi
using DifferentialEquations
using CSV
using Distributions
using Plots
using StatsBase


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

# Lower bounds of the parameters of the ODE
P_GUESS = [0.01, 0.001, 1.00, 1]
ub_ahpm = P_GUESS.*4
lb_ahpm = P_GUESS./4
# Performing ODE fitting WITHOUT RESTART 
@time results_ODE_fit = fitting_one_well_ODE_constrained(
    data_OD, 
    "test",
    "test_ODE",
    model,
    P_GUESS;
)

Plots.plot!(results_ODE_fit[4],results_ODE_fit[3], xlabel="Time", ylabel="Arb. Units", label=["fit " nothing],color=:red,markersize =2 ,size = (300,300))

# Performing ODE fitting WITHOUT BOUNDS AND RESTART 

@time results_ODE_fit = fitting_one_well_ODE_constrained(
    data_OD, 
    "test",
    "test_ODE",
    model,
    P_GUESS;
    multistart = true,
   lb = lb_ahpm,
   ub = ub_ahpm

)

Plots.plot!(results_ODE_fit[4],results_ODE_fit[3], xlabel="Time", ylabel="Arb. Units", label=["fit with restart " nothing],color=:blue,markersize =2 ,size = (300,300))



# performing fitting 
# Custom ODE function
function ODE_custom(du, u, param, t)
    du[1] = u[1] * (1 - u[1]) * param[2] + param[1] * u[1]
end
custom_ub = [1.2, 1.1]
custom_lb = [0.0001, 0.00000001]

param_guess = [0.01,2.0]

# Performing custom ODE fitting without restart
@time results_ODE_fit = fitting_one_well_custom_ODE(
    data_OD, # dataset first row times second row OD
   "test", # name of the well
    "test custom ode", #label of the experiment
   ODE_custom, # ode model to use
   param_guess,
    1; # number ode in the system
  )

Plots.plot!(results_ODE_fit[4],results_ODE_fit[3], xlabel="Time", ylabel="Arb. Units", label=["fit custom " nothing],color=:blue,markersize =2 ,size = (300,300))

# Performing custom ODE fitting with restart

@time results_ODE_fit = fitting_one_well_custom_ODE(
    data_OD, # dataset first row times second row OD
   "test", # name of the well
    "test custom ode", #label of the experiment
   ODE_custom, # ode model to use
   param_guess,
    1; # number ode in the system
    multistart = true,
   lb = custom_lb,
   ub = custom_ub    )

Plots.plot!(results_ODE_fit[4],results_ODE_fit[3], xlabel="Time", ylabel="Arb. Units", label=["fit custom restart " nothing],markersize =2 ,size = (300,300))

# Number of steps for Morris sensitivity analysis
n_step_sensitivity = 2
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

plot(sensitivity_test[2][3,2:end])
plot(sensitivity_test[2][4,2:end])
plot(sensitivity_test[2][5,2:end])
plot(sensitivity_test[2][6,2:end])
plot(sensitivity_test[2][7,2:end])
plot(sensitivity_test[2][8,2:end])
plot(sensitivity_test[2][9,2:end])

@time sensitivity_test = one_well_morris_sensitivity(
    data_OD, 
    "test",
     "test_sensitivity",
      "aHPM", 
      lb_ahpm,
       ub_ahpm,
    N_step_morris=n_step_sensitivity ;
    multistart=true
)

plot(sensitivity_test[2][3,2:end])
plot(sensitivity_test[2][4,2:end])
plot(sensitivity_test[2][5,2:end])
plot(sensitivity_test[2][6,2:end])
plot(sensitivity_test[2][7,2:end])
plot(sensitivity_test[2][8,2:end])
plot(sensitivity_test[2][9,2:end])
# Models candidates and their parameter bounds
list_of_models = ["aHPM",   "baranyi_roberts"]
ub_1 =[ 0.1 , 0.1 , 0.1 , 5.001    ]
lb_1 =[ 0.0001 , 0.001,0.0  , 0.001 ]
p1_guess = lb_1 .+ (ub_1.-lb_1)./2




ub_4 =[ 1.2 , 5.1 , 500.0  ,5.0,5.0  ]
lb_4 =[ 0.0001 , 0.1, 0.00 ,0.2,0.2 ]
p4_guess = lb_4 .+ (ub_4.-lb_4)./2

list_ub = [ub_1, ub_4]
list_lb =[lb_1,  lb_4]
list_guess = [p1_guess,  p4_guess]

# Performing model selection without box
results_ms = ODE_Model_selection(
    data_OD,
    "test", 
    "test_model_selection",
    list_of_models,
    list_guess;
)

# Performing model selection with box and multistart
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
# testing NL fitting

nl_model = ["NL_Richards"]
p_guess = [[1.0,1.0,0.01,300.0]]
lb_nl =[[0.01,0.01,0.000001,00.01]]
ub_nl =[p_guess.*3]




# single fit
@time nl_fit = NL_model_selection(data_OD, # dataset first row times second row OD
"test", 
"test_model_selection",
nl_model, #  model to use
p_guess;
)

Plots.plot!(nl_fit[4],nl_fit[3], xlabel="Time", ylabel="Arb. Units", label=["fit custom restart " nothing],markersize =2 ,size = (300,300))

# single fit with restart
@time nl_fit =  NL_model_selection(data_OD, # dataset first row times second row OD
"test", 
"test_model_selection",
nl_model, #  model to use
p_guess;
multistart = true,
lb_param_array =lb_nl,
ub_param_array = ub_nl

)

Plots.plot!(nl_fit[4],nl_fit[3], xlabel="Time", ylabel="Arb. Units", label=["fit custom restart " nothing],markersize =2 ,size = (300,300))

# Bootstrap

# sensitivity_test

# model selection


