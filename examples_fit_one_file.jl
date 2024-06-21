using Kimchi
using DifferentialEquations
using CSV
using Distributions
using StatsBase
using OptimizationBBO
using Optimization
using OptimizationOptimJL
path_to_data = string("E:/Lavoro/Monod_AA_detection/exp_s7/channel_1.csv")
path_to_annotation = string("E:/Lavoro/Monod_AA_detection/exp_s7/annotation.csv")

fit_log_lin = fit_one_file_Log_Lin(
    " ", #label of the experiment
    path_to_data; # path to the folder to analyze
    path_to_annotation=path_to_annotation,# path to the annotation of the wells
    avg_replicate=true, # if true the average between replicates is fitted. If false all replicate are fitted indipendelitly
  )

# fitting ODE
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
    integrator=Tsit5(), # selection of sciml integrator
    lb = lb_param,
    ub =ub_param
)



fit_od = fit_file_ODE(
    "test", #label of the experiment
    path_to_data, # path to the folder to analyze
    model, # string of the used model
    param_guess;
    path_to_annotation=path_to_annotation,# path to the annotation of the wells
    integrator=Tsit5(), # selection of sciml integrator
    lb = lb_param,
    ub =ub_param,
    multistart = true,
    n_restart = 10,
    abstol = 0.0001
)
# fitting MOdel selection
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


@time ms_file = ODE_model_selection_file(
    "", #label of the experiment
    path_to_data, # path to the folder to analyze
    list_of_models, # ode model to use 
    list_guess;
    lb_param_array=list_lb_param, # lower bound param
    ub_param_array=list_ub_param, # upper bound param
    path_to_annotation=path_to_annotation,# path to the annotation of the wells

)


@time ms_file = ODE_model_selection_file(
    "", #label of the experiment
    path_to_data, # path to the folder to analyze
    list_of_models, # ode model to use 
    list_guess;
    lb_param_array=list_lb_param, # lower bound param
    ub_param_array=list_ub_param, # upper bound param
    path_to_annotation=path_to_annotation,# path to the annotation of the wells
    multistart=true,
    n_restart=20,
)


# fitting segmented ODE

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
    multistart=false,
    n_restart=50,
    lb_param_array=list_lb_param, # lower bound param
    ub_param_array=list_ub_param, # upper bound param
)

# fitting NL

# fitting Segmented BL