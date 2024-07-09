using Kimchi


include("E:/Lavoro/Kimchi.jl-main/KimchiPlot.jl-main/src/function_to_plot.jl")
path_to_data = string("E:/Lavoro/detection_Monod_laws/data/Exp_1/channel_1.csv")
path_to_annotation = string("E:/Lavoro/detection_Monod_laws/data/Exp_1/annotation.csv")


path_to_calib = string("E:/Lavoro/detection_Monod_laws/data/cal_curve.csv")


path_to_plot = string("E:/Lavoro/detection_Monod_laws/seg/")
path_to_results = string("E:/Lavoro/detection_Monod_laws/seg/")


model1 = "HPM_exp"

lb_param1 = [0.00001, 0.000001]
ub_param1 =[0.5,       1.5]
param_guess1 =[0.01, 0.01]
    
model2 = "logistic"

lb_param2 = [0.00001, 0.000001]
ub_param2 =[0.5,       2.5,    ]
param_guess2 =[0.01, 1.01]
    

list_of_models = [model1,model2]
list_guess=  [param_guess1,param_guess2]
list_lb=  [lb_param1,lb_param2]
list_ub=  [ub_param1,ub_param2]


fit_file =  segmentation_ODE_file(
    "seg_exp_1", #label of the experiment
    path_to_data, # path to the folder to analyze
    list_of_models, # ode model to use 
    list_guess, #  param
    1;
    path_to_annotation=path_to_annotation,# path to the annotation of the wells
    detect_number_cpd=false,
    fixed_cpd=true,
    path_to_calib=path_to_calib,
    multiple_scattering_correction=false, # if true uses the given calibration curve to fix the data
    type_of_curve="deriv",
    pt_smooth_derivative= 10,
    pt_avg  = 4,
    verbose = true,
    win_size=12, # numebr of the point to generate intial condition
    smoothing =true,
    lb_param_array=list_lb, # lower bound param
    ub_param_array=list_ub, # upper bound param
    maxiters = 20000000
)


plot_fit_of_file(
        fit_file; # path to the folder to analyze
        path_to_plot=path_to_plot, # path where to save Plots
        display_plots=true,# display plots in julia or not
        save_plots=true, # save the plot or not
        x_size=500,
        y_size =700,
        guidefontsize=18,
        tickfontsize=16,
        legendfontsize=10
)




path_to_data = string("E:/Lavoro/detection_Monod_laws/data/Exp_2/channel_1.csv")
path_to_annotation = string("E:/Lavoro/detection_Monod_laws/data/Exp_2/annotation.csv")

fit_file =  segmentation_ODE_file(
    "seg_exp_2", #label of the experiment
    path_to_data, # path to the folder to analyze
    list_of_models, # ode model to use 
    list_guess, #  param
    1;
    path_to_annotation=path_to_annotation,# path to the annotation of the wells
    detect_number_cpd=false,
    fixed_cpd=true,
    path_to_calib=path_to_calib,
    multiple_scattering_correction=true, # if true uses the given calibration curve to fix the data
    type_of_curve="deriv",
    pt_smooth_derivative= 10,
    verbose = true,
    pt_avg  = 4,
    win_size=12, # numebr of the point to generate intial condition
    smoothing =true,
    lb_param_array=list_lb, # lower bound param
    ub_param_array=list_ub, # upper bound param
    maxiters = 2000000
)


plot_fit_of_file(
        fit_file; # path to the folder to analyze
        path_to_plot=path_to_plot, # path where to save Plots
        display_plots=true,# display plots in julia or not
        save_plots=true, # save the plot or not
        x_size=500,
        y_size =700,
        guidefontsize=18,
        tickfontsize=16,
        legendfontsize=10
)