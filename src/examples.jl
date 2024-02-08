
# uploading package this part of the code will be skipped when the github is online 
# you have only to write using JMAKI
using DifferentialEquations, Optimization, Plots, Random
using CSV, DataFrames
using Statistics
using Optim
#using Zygote
using OptimizationBBO
using NaNMath
using CurveFit
using StatsBase
using Tables
using Distributions
using Interpolations
using GaussianProcesses
using Peaks
using ChangePointDetection

# loading the functions  (will be skipped when the package is deployed)
path_to_functions = "/Users/fabrizio.angaroni/Documents/J-MAKi.jl-main.jl/src/";
include(string(path_to_functions, "/functions.jl"))


# simulating data with a ODE
# string of the model see the github for a partial list
model = "triple_piecewise_damped_logistic"
# initial contion of the ODE, note sciML requires vectors
n_start = [0.1]
# starting time of the ODE
tstart = 0.0
#ending time of the ODE
tmax = 600.0
# delta t for numerical integration
delta_t = 10.0
# sciML numerical integrator
integrator = KenCarp4()

# parameters of the ODE 
param_of_ode = [
    0.06, # growth rate
    1.0, # carrying capacity
    200, # lag time
    0.5, # shape factor 
    0.001, # linear factor 
    450, # stationary start
    -0.0002,# linear param of stationary phase]
]


# calling the function simulation

sim = ODE_sim(
    model, #string of the model
    n_start, # starting condition
    tstart, # start time of the sim
    tmax, # final time of the sim
    delta_t, # delta t for poisson approx
    integrator, # which sciml solver of ode
    param_of_ode, # parameters of the ODE model
)


# ploting scatterplot of data
Plots.scatter(
    sim,
    xlabel = "Time",
    ylabel = "Arb. Units",
    label = ["Data " nothing],
    color = :blue,
    size = (300, 300),
)

#adding uniform random noise
noise_unifom = rand(Uniform(-0.05, 0.05), length(sim.t))


data_t = reduce(hcat, sim.t)
data_o = reduce(hcat, sim.u)
data_OD = vcat(data_t, data_o)
data_OD[2, :] = data_OD[2, :] .+ noise_unifom
# ploting scatterplot of data with noise

Plots.scatter(
    data_OD[1, :],
    data_OD[2, :],
    xlabel = "Time",
    ylabel = "Arb. Units",
    label = ["Data " nothing],
    color = :blue,
    markersize = 2,
    size = (300, 300),
)

# smooting of the data with rolling average
# it is possible to do that also with gaussian processes (gaussian_smoothing(data,optimize_gp))
data_OD_smooth = smoothing_data(data_OD, 7)
data_OD_smooth = Matrix(data_OD_smooth)
Plots.scatter(
    data_OD_smooth[1, :],
    data_OD_smooth[2, :],
    xlabel = "Time",
    ylabel = "Arb. Units",
    label = ["Smoothed data " nothing],
    markersize = 2,
    color = :blue,
    size = (300, 300),
)

# multiple scattering correction, since requires an external file you can comment this part
data_OD_smooth2 = correction_OD_multiple_scattering(
    data_OD_smooth,
    "/Users/fabrizio.angaroni/Documents/calib_OD600_1ml/ms_cal_v1._ten.csv",
)
Plots.scatter(
    data_OD_smooth2[1, :],
    data_OD_smooth2[2, :],
    xlabel = "Time",
    ylabel = "Arb. Units",
    label = ["Pre-processes data" nothing],
    markersize = 2,
    color = :blue,
    size = (300, 300),
)

# log lin fitting
# first, as example, I evaluate the specific growth rate (i.e. the dynamic of the slope of a log linear fitting of a slinding window)
# size of the win. If 0 the methods involve interpolation between the data point 
pt_smoothing_derivative = 7
deriv = specific_gr_evaluation(data_OD_smooth2, 7)
# to be correctet if 0
specific_gr_times = [
    (data_OD_smooth2[1, r] + data_OD_smooth2[1, (r+pt_smoothing_derivative)]) / 2 for
    r = 1:1:(length(data_OD_smooth2[2, :])-pt_smoothing_derivative)
]
data_deriv = Matrix(transpose(hcat(specific_gr_times, deriv)));
Plots.plot(
    data_deriv[1, :],
    data_deriv[2, :],
    xlabel = "Time",
    ylabel = "Specific GR [1/Time]",
    label = nothing,
    markersize = 2,
    color = :blue,
    size = (300, 300),
)

# now i fit for real the data with log lin
# path where save the results plots
res_log_lin = fitting_one_well_Log_Lin(
    data_OD_smooth2, # dataset first row times second row OD
    "test", # name of the well
    "test"; #label of the experiment
)

# ODE fitting


# upper bounds of the parameters of the ODE
ub_dhpm = [1.2, 1.1, 2.0, 20]
# lowert bounds of the parameters of the ODE

lb_dhpm = [0.0001, 0.00000001, 0.00, 0]

# starting guess (i start in the middle. The argument will be changed)

test_param = (ub_dhpm .- lb_dhpm) .+ lb_dhpm ./ 2
path_to_plotting = "/Users/fabrizio.angaroni/Documents/test_jmaki/"

results = fitting_one_well_ODE_constrained(
    data_OD_smooth2, # dataset first row times second row OD
    "", # name of the well
    "", #label of the experiment
    "dHPM", # ode model to use 
    lb_dhpm, # lower bound param
    ub_dhpm; # upper bound param
    do_plot = true, # do plots or no
    path_to_plot = path_to_plotting, # where save plots
)

# testing custom ode

# function of the custom ODE
function ODE_custom(du, u, param, t)

    du[1] = u[1] * (1 - u[1]) * (param[2]) + param[1] * u[1]

    du[2] = +u[1] * (param[2]) + param[4] * u[2] * (1 - (u[1] + u[2]) / param[3])


end

custom_ub = [1.2, 1.1, 2.0, 20]

custom_lb = [0.0001, 0.00000001, 0.00, 0]


results = fitting_one_well_custom_ODE(
    data_OD_smooth2, # dataset first row times second row OD
    "test", # name of the well
    "test_model_custom", #label of the experiment
    ODE_custom, # ode model to use 
    custom_lb, # lower bound param
    custom_ub, # upper bound param
    2; # number ode in the system
    do_plot = true, # do plots or no
    path_to_plot = path_to_plotting, # where save plots
)

# testing sensitivity
n_step_sensitivity = 3

sensitivity_test = one_well_morris_sensitivity(
    data_OD_smooth2, # dataset first row times second row OD
    "test", # name of the well
    "test_sensitivity", #label of the experiment
    "dHPM", # ode model to use 
    lb_dhpm, # lower bound param
    ub_dhpm; # upper bound param
    N_step_morris = n_step_sensitivity,
)
# testing model selection



ub_piece_wise_logistic = [0.06, 2.0, 500.0, 10.0, 0.001]
lb_piece_wise_logistic = [0.0001, 0.001, 0.0, 0.001, -0.001]

ub_baranyi_roberts = [0.06, 2.0, 500.0, 10.0, 10]
lb_baranyi_roberts = [0.0001, 0.001, 0.0, 0.01, 0]



ub_triple_exp = [1.2, 0.001, 0.2, 500.0, 2000]
lb_triple_exp = [0.0001, -0.001, 0.0, 00.0, 200.0]



ub_dhpm = [1.2, 1.1, 2.0, 20]
lb_dhpm = [0.0001, 0.00000001, 0.00, 0]



list_of_models =
    ["dHPM", "piecewise_damped_logistic", "triple_piecewise", "baranyi_roberts"]
list_ub = [ub_dhpm, ub_piece_wise_logistic, ub_triple_exp, ub_baranyi_roberts];
list_lb = [lb_dhpm, lb_piece_wise_logistic, lb_triple_exp, lb_baranyi_roberts];





results_ms = ODE_Model_selection(
    data_OD_smooth2, # dataset first row times second row OD
    "test", # name of the well
    "test_model_selection", #label of the experiment
    list_of_models, # ode model to use 
    list_lb, # lower bound param
    list_ub; # upper bound param
    plot_best_model = true, # one wants the results of the best fit to be plotted
    path_to_plot = path_to_plotting,
    verbose = true,
)









# now i show one of the segmentation fitting
# first i generate a dataset with a sequence of ODEs

# first segment ODE


model = "exponential"
n_start = [0.1]
tstart = 0.0

tmax = 010.0
delta_t = 2.0
param_of_ode = [0.00]
sim_1 = ODE_sim(
    model, #string of the model
    n_start, # starting condition
    tstart, # start time of the sim
    tmax, # final time of the sim
    delta_t, # delta t for poisson approx
    integrator, # which sciml solver of ode
    param_of_ode, # parameters of the ODE model
)

sol_1 = reduce(vcat, sim_1)

# second segment ODE


model = "logistic"
n_start = [sol_1[end]]
tstart = 10.0
tmax = 080.0
delta_t = 2.0
param_of_ode = [0.2, 0.4]


sim_2 = ODE_sim(
    model, #string of the model
    n_start, # starting condition
    tstart, # start time of the sim
    tmax, # final time of the sim
    delta_t, # delta t for poisson approx
    integrator, # which sciml solver of ode
    param_of_ode, # parameters of the ODE model
)

sol_2 = reduce(vcat, sim_2)
# third segment ODE


model = "logistic"
n_start = [sol_2[end]]
tstart = 80.0
tmax = 0200.0
delta_t = 2.0
param_of_ode = [0.1, 0.8]
sim_3 = ODE_sim(
    model, #string of the model
    n_start, # starting condition
    tstart, # start time of the sim
    tmax, # final time of the sim
    delta_t, # delta t for poisson approx
    integrator, # which sciml solver of ode
    param_of_ode, # parameters of the ODE model
)

sol_3 = reduce(vcat, sim_3)
times_sim = vcat(sim_1.t, sim_2.t)
times_sim = vcat(times_sim, sim_3.t)

# binding the simulatios
sol_sim = vcat(sol_1, sol_2)
sol_sim = vcat(sol_sim, sol_3)


Plots.scatter(
    sol_sim,
    xlabel = "Time",
    ylabel = "Arb. Units",
    label = ["Data " nothing],
    color = :blue,
    size = (300, 300),
)


data_OD = Matrix(transpose(hcat(times_sim, sol_sim)))

cpd_lsdd_profile(data_OD, 2)

deriv = specific_gr_evaluation(data_OD, 0)
Plots.scatter(
    deriv,
    xlabel = "Time",
    ylabel = "Arb. Units",
    label = ["Data " nothing],
    color = :blue,
    size = (300, 300),
)




profile = ChangePointDetection.lsdd_profile(deriv; window = 2)
Plots.scatter(
    profile,
    xlabel = "Time",
    ylabel = "Arb. Units",
    label = ["Data " nothing],
    color = :blue,
    size = (300, 300),
)


ChangePointDetection.changepoints(deriv; threshold = mean(profile), window = 2)
ChangePointDetection.changepoints(data_OD; threshold = mean(profile), window = 2)

profile = convert.(Float64, profile)
data_dissim = Matrix(transpose(hcat(data_OD[1, 1:length(profile)], profile)))

pks, vals = findmaxima(data_dissim[2, :]; strict = true)
pks, proms = peakproms(pks, data_dissim[2, :])
pks, widths, leftedge, rightedge = peakwidths(pks, data_dissim[2, :], proms)

selected_change_point_index = peaks_detection(data_dissim, 2)

Plots.vline!(
    selected_change_point_index[1],
    xlabel = "Time",
    ylabel = "Arb. Units",
    label = ["Data " nothing],
    color = :red,
    size = (300, 300),
)


#adding noise_unifo
noise_unifom = rand(Uniform(-0.01, 0.01), length(sol_sim))


data_OD = Matrix(transpose(hcat(times_sim, sol_sim)))
data_OD[2, :] = data_OD[2, :] .+ noise_unifom
Plots.scatter(
    data_OD[1, :],
    data_OD[2, :],
    xlabel = "Time",
    ylabel = "Arb. Units",
    label = ["Data " nothing],
    color = :blue,
    markersize = 2,
    size = (300, 300),
)
# testing change point function 

# inizialixating all the models
ub_exp = [0.1]
lb_exp = [-00.01]

ub_logistic = [0.9, 5.0]
lb_logistic = [0.0001, 0.001]

ub_hpm = [0.1, 20.0, 50.001]
lb_hpm = [0.0001, 0.000001, 0.001]

ub_hpm_exp = [0.1, 20.0]
lb_hpm_exp = [0.0001, 0.0000001]

list_of_models = ["exponential", "HPM", "HPM_exp", "logistic"]#,"dHPM"]
list_ub_param = [ub_exp, ub_hpm, ub_hpm_exp, ub_logistic]#,ub_dhpm]# , ub_triple_exp ,ub_baranyi_roberts  ];
list_lb_param = [lb_exp, lb_hpm, lb_hpm_exp, lb_logistic]#,lb_dhpm]




# fitting with a fixed number of change points

test_fixed_cdp = selection_ODE_fixed_change_points(
    data_OD, # dataset first row times second row OD
    "test", # name of the well
    "test", #label of the experiment
    list_of_models, # ode models to use 
    list_lb_param, # lower bound param
    list_ub_param, # upper bound param
    3; # number of change_points
    do_plot = true, # do plots or no
    path_to_plot = path_to_plotting, # where save plots
    pt_smooth_derivative = 0,
)
