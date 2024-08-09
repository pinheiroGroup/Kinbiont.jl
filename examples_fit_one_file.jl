using Kinbiont
using DifferentialEquations
using CSV
using Distributions
using StatsBase
using OptimizationBBO
using Optimization
using OptimizationOptimJL
path_to_data = string("/Users/fabrizio.angaroni/Documents/JMAKi_utilities/real_dataset_tests/dataset/Monod_AA_detection/exp_s7/channel_1.csv")
path_to_annotation = string("/Users/fabrizio.angaroni/Documents/JMAKi_utilities/real_dataset_tests/dataset/Monod_AA_detection/exp_s7/annotation.csv")

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
    ub =ub_param,
    maxiters =5
)



fit_od = fit_file_ODE(
    "test", #label of the experiment
    path_to_data, # path to the folder to analyze
    model, # string of the used model
    param_guess;
    path_to_annotation=path_to_annotation,# path to the annotation of the wells
    integrator=Tsit5(), # selection of sciml integrator
    #lb = lb_param,
   # ub =ub_param,
  #  multistart = true,
    n_restart = 10,
    abstol = 0.0001,
    maxiters =5
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
    maxiters =5
)


@time ms_file = ODE_model_selection_file(
    "", #label of the experiment
    path_to_data, # path to the folder to analyze
    list_of_models, # ode model to use 
    list_guess;
    path_to_annotation=path_to_annotation,# path to the annotation of the wells
    maxiters =5
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
     maxiters =5,
    lb_param_array=list_lb_param, # lower bound param
    ub_param_array=list_ub_param, # upper bound param
)

# fitting NL

nl_model = ["NL_Richards"]
p_guess = [[1.0,1.0,0.01,300.0]]
lb_nl =[[0.01,0.01,0.000001,00.01]]
ub_nl =p_guess.*50

@time fit_nl = fit_NL_model_selection_file(
    "TEST", #label of the experiment
    path_to_data    , # path to the folder to analyze
    nl_model, # ode model to use
    p_guess;# initial guess param
    lb_param_array =lb_nl, # lower bound param
    ub_param_array = ub_nl, # upper bound param
    path_to_annotation = path_to_annotation,# path to the annotation of the wells
  
)


@time fit_nl = fit_NL_model_selection_file(
    "TEST", #label of the experiment
    path_to_data    , # path to the folder to analyze
    nl_model, # ode model to use
    p_guess;# initial guess param
    path_to_annotation = path_to_annotation,# path to the annotation of the wells
  
)


# testing model selection


nl_model = ["NL_Richards","NL_Bertalanffy"]
p_guess = [[1.0,1.0,0.01,300.0],[0.08,1.0,0.01,1.0]]
lb_nl =[[0.01,0.01,0.000001,00.01],[0.00,0.01,0.001,00.01]]
ub_nl =p_guess.*50

@time fit_nl = fit_NL_model_selection_file(
    "TEST", #label of the experiment
    path_to_data    , # path to the folder to analyze
    nl_model, # ode model to use
    p_guess;# initial guess param
    lb_param_array =lb_nl, # lower bound param
    ub_param_array = ub_nl, # upper bound param
    path_to_annotation = path_to_annotation,# path to the annotation of the wells
  
)


# fitting Segmented BL
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