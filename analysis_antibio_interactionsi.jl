using Kimchi
using CSV
using Plots
using Tables

df =CSV.File("E:/Lavoro/JMAKi_utilities-main/real_dataset_tests/dataset/species-specific_activity_antibio/Supplementary File 9.txt")
names_of_cols = propertynames(df)
path_to_res ="E:/Lavoro/JMAKi_utilities-main/real_dataset_tests/dataset/species-specific_activity_antibio/res/"
mkpath(path_to_res)
position_of_starts = findall(df[names_of_cols[3]].==0)
# loop on the plate
# segm

df[findall(df[:Query] .== "Chloramphenicol 1 rep 1")]
plot(df[findall(df[:Query] .== "Chloramphenicol 1 rep 1")])





for p in position_of_starts
  index = findfirst(position_of_starts .== p)
  start_in =position_of_starts[index]

    results_Log_Lin = [
        "label_exp",
        "well_name",
        "t_start",
        "t_end",
        "t_of_max",
        "empirical_max_Growth_rate",
        "Growth_rate",
        "sigma_gr",
        "dt",
        "sigma_confidence_dt_upper",
        "sigma_confidence_dt_lower",
        "intercept",
        "sigma_intercept",
        "Pearson_correlation",
    ]

  for nn in names_of_cols[4:end]

       if p!= position_of_starts[end]

        end_in  = position_of_starts[index+1]
        times = df[names_of_cols[3]][start_in : (end_in-1)]
        data = df[nn][start_in: (end_in-1)]
       else
        times = df[names_of_cols[3]][start_in : end]
         data = df[nn][start_in:end]
       end

        
    data_OD = Matrix(transpose(hcat(times, data)))

    fit_log_lin = fitting_one_well_Log_Lin(
        data_OD, # dataset first row times second row OD
        string(nn), # name of the well
        string(df[names_of_cols[1]][start_in],"-",df[names_of_cols[2]][start_in]); #label of the experiment
        type_of_smoothing="rolling_avg", # option, NO, gaussian, rolling avg
        pt_avg=3, # number of the point for rolling avg not used in the other cases
        pt_smoothing_derivative=4, # number of poits to smooth the derivative
        pt_min_size_of_win=4, # minimum size of the exp windows in number of smooted points
        type_of_win="maximum", # how the exp. phase win is selected, "maximum" of "global_thr"
        threshold_of_exp=0.9, # threshold of growth rate in quantile to define the exp windows
        start_exp_win_thr=0.05, # minimum value to consider the start of exp window
    )
     
    results_Log_Lin = hcat(results_Log_Lin, fit_log_lin[2])

  end

   
  label_exp = string(df[names_of_cols[1]][start_in],"-",df[names_of_cols[2]][start_in])
  CSV.write(
    string(path_to_res, label_exp, "_results.csv"),
    Tables.table(Matrix(results_Log_Lin)),)


end   

path_to_res ="E:/Lavoro/JMAKi_utilities-main/real_dataset_tests/dataset/species-specific_activity_antibio/res_deriv/"
mkpath(path_to_res)

for p in position_of_starts
  index = findfirst(position_of_starts .== p)
  start_in =position_of_starts[index]
    results_deriv_ = ["label_exp", "name_well", "max_gr", "min_gr", "t_of_max", "od_of_max", "max_deriv", "min_deriv", "end_time", "delta_N", "segment_number"]

  for nn in names_of_cols[4:end]

       if p!= position_of_starts[end]

        end_in  = position_of_starts[index+1]
        times = df[names_of_cols[3]][start_in : (end_in-1)]
        data = df[nn][start_in: (end_in-1)]
       else
        times = df[names_of_cols[3]][start_in : end]
         data = df[nn][start_in:end]
       end

        
    data_OD = Matrix(transpose(hcat(times, data)))

    results_deriv = segment_gr_analysis(
        data_OD, # dataset first row times second row OD
        string(nn), # name of the well
        string(df[names_of_cols[1]][start_in],"-",df[names_of_cols[2]][start_in]); #label of the experiment
        type_of_smoothing="rolling_avg", # option, NO, gaussian, rolling avg
        pt_avg=3, # number of the point for rolling avg not used in the other cases
        pt_smoothing_derivative=4, # number of poits to smooth the derivative

    )
     
    results_deriv_ = hcat(results_deriv_, results_deriv[2])

  end
   
  label_exp = string(df[names_of_cols[1]][start_in],"-",df[names_of_cols[2]][start_in])
  
   CSV.write(
    string(path_to_res, label_exp, "_results.csv"),
    Tables.table(Matrix(results_deriv_)),)


end   



nl_model = ["NL_piecewise_logistic"]

  

lb_nl =[[0.001,0.00001,0.000001,00.0001,-0.001]]


results_fit = Kimchi.initialize_res_ms(lb_nl)


path_to_res ="E:/Lavoro/JMAKi_utilities-main/real_dataset_tests/dataset/species-specific_activity_antibio/res_fit/"
mkpath(path_to_res)
path_to_plot ="E:/Lavoro/JMAKi_utilities-main/real_dataset_tests/dataset/species-specific_activity_antibio/plot_fit/"
mkpath(path_to_plot)

count = 0
for p in position_of_starts
  index = findfirst(position_of_starts .== p)
  start_in =position_of_starts[index]
    results_deriv_ = ["label_exp", "name_well", "max_gr", "min_gr", "t_of_max", "od_of_max", "max_deriv", "min_deriv", "end_time", "delta_N", "segment_number"]

  for nn in names_of_cols[4:end]
    count = count+1

       if p!= position_of_starts[end]

        end_in  = position_of_starts[index+1]
        times = df[names_of_cols[3]][start_in : (end_in-1)]
        data = df[nn][start_in: (end_in-1)]
       else
        times = df[names_of_cols[3]][start_in : end]
         data = df[nn][start_in:end]
       end

        
      data_OD = Matrix(transpose(hcat(times, data)))

      p_guess = [[times[4],median(data[1:4]) ,maximum(data), maximum(Kimchi.deriv_evaluation(data_OD)),0.001]]
    
      lb_nl =[[0.001,0.00001,0.000001,00.0001,-0.001]]
      ub_nl =p_guess.*10
     nl_fit =  NL_model_selection(data_OD, # dataset first row times second row OD
      "test", 
      "test_model_selection",
      nl_model, #  model to use
      p_guess;
      type_of_loss ="RE",
      smoothing = true,
      type_of_smoothing="lowess",
      #multistart = true,
      lb_param_array =lb_nl,
      ub_param_array = ub_nl

    )
     
    results_fit = hcat(results_fit, nl_fit[2])
  if count%100 == 0
    label_exp = string(df[names_of_cols[1]][start_in],"-",df[names_of_cols[2]][start_in])

    display( scatter(times,data,
    xlabel="Time",
    ylabel="Arb. Units",
    size=(700,450),
    legend=:outertopright,
    guidefontsize=15,
    tickfontsize=15,
    legendfontsize=15,))
    display( plot!(nl_fit[4],nl_fit[3],
     xlabel="Time",
    ylabel="Arb. Units",
    size=(700,450),
    legend=:outertopright,
    guidefontsize=15,
    tickfontsize=15,
    legendfontsize=15,))
    png(string(path_to_plot, label_exp, ".png"))

  end
end
  label_exp = string(df[names_of_cols[1]][start_in],"-",df[names_of_cols[2]][start_in])
  
   CSV.write(
    string(path_to_res, label_exp, "_results.csv"),
    Tables.table(Matrix(results_fit)),)


end