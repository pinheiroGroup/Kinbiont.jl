




function cpd_local_detection(data::Matrix{Float64},
    n_max_cp::Int;
    type_of_detection="lsdd",
    type_of_curve="deriv", 
    pt_derivative = 0,
    size_win =2)
     

    """
        perform the analyses with the required change point detection algorithm
        type_of_detection="lsdd" or piecewise linear fitting on the specific growth rate
        pt_derivative number of point to evaluate the derivative/specific gr (if 0 numerical derivative if >1 specific gr with that size of sliding window)
        size_win Int size of the used window in all of the methods
    """


    if type_of_detection == "lsdd" &&  type_of_curve=="deriv"
        list_of_cdps = cpd_lsdd_profile(data,n_max_cp;pt_deriv=pt_derivative,window_size=size_win)

    elseif type_of_detection == "lsdd"  && type_of_curve !="deriv"

        list_of_cdps = cpd_lsdd_profile(data,n_max_cp;pt_deriv=pt_derivative,window_size=size_win, type_of_curve="original")


    else    
        list_of_cdps= detect_list_change_points(data,n_max_cp;win_size=size_win)

    end


    return list_of_cdps
end







function cpd_lsdd_profile(data::Matrix{Float64},n_max::Int; window_size = 2, type_of_curve= "deriv",pt_deriv=0)
 


    """
    evaluate change points using peak detection on a lsdd profile
    type_of_detection="lsdd" or piecewise linear fitting on the specific growth rate
    pt_derivative number of point to evaluate the derivative/specific gr (if 0 numerical derivative if >1 specific gr with that size of sliding window)
    size_win Int size of the used window in all of the methods
    """

    # evaluating the profile of lsdd on the data or on the derivative of the data
    if type_of_curve == "deriv"

        deriv = specific_gr_evaluation(data,pt_deriv )

        profile = ChangePointDetection.lsdd_profile(deriv; window = window_size)

    else

        profile = ChangePointDetection.lsdd_profile(data[2,:]; window = window_size)

    end

    # adding time to profile
    profile = convert.(Float64,profile)
    data_dissim = Matrix(transpose(hcat(data[1,1:length(profile)], profile)))

    selected_change_point_index = peaks_detection(data_dissim,n_max)
   
    return selected_change_point_index
end



function detect_list_change_points( data::Matrix{Float64},n_max::Int;win_size=2)

  
   
        """
         evaluate change points using piecewise linear fitting 
         n_max
         size_win Int size of the used window in all of the methods
        """
   
   
       curve_dissimilitary_deriv = curve_dissimilitary_lin_fitting( data, # dataset x times y OD/fluorescence
       1, # index of start
       win_size, # size sliding window
       )
       data_dissim = Matrix(transpose(hcat(data[1,convert.(Int,curve_dissimilitary_deriv[1,:])], curve_dissimilitary_deriv[2,:])));
   
           
       selected_change_point_index = peaks_detection(data_dissim,n_max)
   
       return selected_change_point_index
         
   
   
       
   end
function detect_list_change_points_derivative( data::Matrix{Float64},n_max::Int,win_size::Int,pt_smoothing_derivative::Int)

 if pt_smoothing_derivative == 0


    derivative_interpolation = specific_gr_interpol_evaluation(data)
    data_deriv =      Matrix(transpose(hcat(data[1,:], derivative_interpolation)));

 else

    derivative_interpolation = specific_gr_evaluation(data,pt_smoothing_derivative)

    specific_gr_times = [(data[1, r] + data[1, (r+pt_smoothing_derivative)]) / 2 for r in 1:1:(length(data[2, :])-pt_smoothing_derivative)]
   
    data_deriv =      Matrix(transpose(hcat(specific_gr_times, derivative_interpolation)));

 end



    curve_dissimilitary_deriv = curve_dissimilitary_lin_fitting( data_deriv, # dataset x times y OD/fluorescence
    1, # index of start
    win_size, # size sliding window
    )
    data_dissim = Matrix(transpose(hcat(data_deriv[1,convert.(Int,curve_dissimilitary_deriv[1,:])], curve_dissimilitary_deriv[2,:])));

        
    selected_change_point_index = peaks_detection(data_dissim,n_max)

    return selected_change_point_index
      


    
end

function specific_gr_interpol_evaluation(data_testing)
    
    itp = interpolate((data_testing[1,:],), data_testing[2,:], Gridded(Linear()));
    specific_gr_interpol = only.(Interpolations.gradient.(Ref(itp), data_testing[1,:]))
    return specific_gr_interpol
end

function peaks_detection(data::Matrix{Float64},
    n_max::Int)

    """
    peaks detection
    n_max maximum number of peaks 
    size_win Int size of the used window in all of the methods
   """

    index_of_peaks =  findmaxima(data[2,:]; strict=true)

    array_prominence = peakproms(index_of_peaks[1],data[2,:])[2]
    index_prominence = peakproms(index_of_peaks[1],data[2,:])[1]
    if length(array_prominence) < n_max
        println("Warning: the max number of peaks is too much")
        top_prominence = sort(array_prominence)
    
    else
        top_prominence = sort(array_prominence)[((end - n_max)+1):end]

    end
    
        index_top_peaks = [ findall(array_prominence .== i)[1] for i in top_prominence]
    times_top_peaks = data[1,    index_prominence[index_top_peaks]    ]
    values_top_peaks = data[2,    index_prominence[index_top_peaks]    ]

    return index_prominence[index_top_peaks] ,times_top_peaks,values_top_peaks
end

function curve_dissimilitary_lin_fitting( data::Matrix{Float64}, # dataset x times y OD/fluorescence
    start_time_Index::Int,
    size_wind::Int, # size sliding window
)
    discrepancy_measure_curve = [start_time_Index,0.0]
    ending = convert(Int, (length(data[2,:]) - floor(size_wind/2)*2))

    for index_t in start_time_Index:ending
        # defining the window
        middle_index= convert(Int,(index_t+floor(size_wind/2)))
        end_index = convert(Int,(index_t+floor(size_wind/2)*2))

        win_1_data = data[2,index_t:middle_index]
        win_2_data  = data[2,middle_index:end_index]
        win_tot_data  =  data[2,index_t:end_index]

        win_1_times = data[1,index_t:middle_index]
        win_2_times   = data[1,middle_index:end_index]
        win_tot_times   =  data[1,index_t:end_index]
    

        #fitting total data
        data_total =      Matrix(transpose(hcat(win_tot_times, win_tot_data)));
        fit_total =  curve_fit(LinearFit, data_total[1,:], data_total[2,:])   
        # residual calculation    
        res_total          = sum([ abs((data_total[2,ll] -    fit_total.coefs[2] *data_total[1,ll] -  fit_total.coefs[1])) for ll in 1 : length(data_total[1,:])   ] )


        
       
    
      #fitting win 1
      data_1 =      Matrix(transpose(hcat(win_1_times, win_1_data)));
      fit_1 = curve_fit(LinearFit, data_1[1,:], data_1[2,:])   
      # residual calculation    
      res_win_1          = sum([ abs((data_1[2,ll] -    fit_1.coefs[2] *data_1[1,ll] -  fit_1.coefs[1])) for ll in 1 : length(data_1[1,:])   ] )


      #fitting win 2
      data_2 =      Matrix(transpose(hcat(win_2_times, win_2_data)));
      fit_2 = curve_fit(LinearFit, data_2[1,:], data_2[2,:])   
      # residual calculation    
      res_win_2          = sum([ abs((data_2[2,ll] -    fit_2.coefs[2] *data_2[1,ll] -  fit_2.coefs[1])) for ll in 1 : length(data_2[1,:])   ] )

        #evaluation of the cost
    
        cost =res_total -res_win_1 -res_win_2
        discrepancy_measure_curve = hcat(discrepancy_measure_curve, [index_t + floor(size_wind/2),cost ])
       # println(discrepancy_measure_curve)    
         # stop when first change point is fitted


       
    end

    return discrepancy_measure_curve




end










function model_selector(model::String,u0,tspan)

    """
    generate sciML OD problem for fitting the ODE
   
   """


    if model == "exponetial"


        ODE_prob = ODEProblem(ODE_exponential, u0, tspan, nothing)


    end
    
    if model == "hyper_gompertz"


        ODE_prob = ODEProblem(hyper_gompertz, u0, tspan, nothing)


    end

    if model == "hyper_logistic"
        ODE_prob = ODEProblem(hyper_logistic, u0, tspan, nothing)

    end

    if model == "bertalanffy_richards"

        ODE_prob = ODEProblem(bertalanffy_richards, u0, tspan, nothing)

    end
    if model == "ode_von_bertalanffy"

        ODE_prob = ODEProblem(ODE_von_bertalanffy, u0, tspan, nothing)

    end
    if model == "triple_piecewise_bertalanffy_richards"

        ODE_prob = ODEProblem(ODE_triple_piecewise_bertalanffy_richards, u0, tspan, nothing)

    end

    if model == "gbsm_piecewise"

        ODE_prob = ODEProblem(ODE_gbsm_piecewise, u0, tspan, nothing)

    end
    if model == "logistic"

        ODE_prob = ODEProblem(logistic, u0, tspan, nothing)

    end 

    if model == "dlogistic"

        ODE_prob = ODEProblem(dlogistic, u0, tspan, nothing)

    end 
    if model == "exponential"

        ODE_prob = ODEProblem(ODE_exponential, u0, tspan, nothing)

    end

    if model == "gompertz"

        ODE_prob = ODEProblem(gompertz, u0, tspan, nothing)

    end

    if model == "baranyi_richards"

        ODE_prob = ODEProblem(baranyi_richards, u0, tspan, nothing)

    end
 

    if model == "baranyi_exp"

        ODE_prob = ODEProblem(baranyi_exp, u0, tspan, nothing)

    end
    if model == "baranyi_roberts"

        ODE_prob = ODEProblem(baranyi_roberts, u0, tspan, nothing)

    end


    if model == "huang"
        ## attention!!!!!
        u0 = [log(u0)]

        ODE_prob = ODEProblem(huang, u0, tspan, nothing)


    end

    if model == "piecewise_damped_logistic"


        ODE_prob = ODEProblem(ODE_piecewise_damped_logistic, u0, tspan, nothing)


    end



    if model == "Diauxic_replicator_1"


        ODE_prob = ODEProblem(ODE_Diauxic_replicator_1, u0, tspan, nothing)


    end
    if model == "Diauxic_replicator_2"


        ODE_prob = ODEProblem(ODE_Diauxic_replicator_2, u0, tspan, nothing)


    end
    if model == "Diauxic_piecewise_damped_logistic"


        ODE_prob = ODEProblem(ODE_Diauxic_piecewise_damped_logistic, u0, tspan, nothing)


    end
    if model == "triple_piecewise_damped_logistic"


        ODE_prob = ODEProblem(ODE_triple_piecewise_damped_logistic, u0, tspan, nothing)


    end


    if model == "triple_piecewise"


        ODE_prob = ODEProblem(ODE_triple_piecewise, u0, tspan, nothing)

    end
    if model == "triple_piecewise_sublinear"


        ODE_prob = ODEProblem(ODE_triple_piecewise_sub_linear, u0, tspan, nothing)

    end



    if model == "four_piecewise"


        ODE_prob = ODEProblem(ODE_four_piecewise, u0, tspan, nothing)

    end
    if model == "HPM_3_death"


        ODE_prob = ODEProblem(HPM_3_death, u0, tspan, nothing)

    end

    if model == "HPM_3_death_resistance"


        ODE_prob = ODEProblem(HPM_3_death_resistance, u0, tspan, nothing)

    end
    
    if model == "dHPM_3_death_resistance"


        ODE_prob = ODEProblem(dHPM_3_death_resistance, u0, tspan, nothing)

    end
    if model == "HPM_3_inhibition"


        ODE_prob = ODEProblem(HPM_3_inhibition, u0, tspan, nothing)

    end
    if model == "HPM"
            
        ODE_prob = ODEProblem(ODEs_McKellar, u0, tspan, nothing)
    end
    if model == "HPM_exp"
            
        ODE_prob = ODEProblem(ODEs_HPM_exp, u0, tspan, nothing)
    end

    if model == "dHPM"
        
        ODE_prob = ODEProblem(ODEs_damped_McKellar, u0, tspan, nothing)
    end
        # defining the loss functions for this case

    if model == "HPM_inhibition"
        
            ODE_prob = ODEProblem(ODEs_HPM_inhibition, u0, tspan, nothing)
    end
    if model == "dHPM_inhibition"
        
           ODE_prob = ODEProblem(ODEs_dHPM_inhibition, u0, tspan, nothing)
    end
    if model == "ODEs_HPM_SR"
        
        ODE_prob = ODEProblem(ODEs_HPM_SR, u0, tspan, nothing)
    end


    return ODE_prob

end   

function model_selector_sim(model::String,u0,tspan,param)
    
    """ 
    generate sciML ODE problem for simulations
   """
    if model == "hyper_gompertz"


        ODE_prob = ODEProblem(hyper_gompertz, u0, tspan, param)


    end

    if model == "hyper_logistic"
        ODE_prob = ODEProblem(hyper_logistic, u0, tspan, param)

    end

    if model == "bertalanffy_richards"

        ODE_prob = ODEProblem(bertalanffy_richards, u0, tspan, param)

    end
    if model == "ode_von_bertalanffy"

        ODE_prob = ODEProblem(ODE_von_bertalanffy, u0, tspan, param)

    end
    if model == "triple_piecewise_bertalanffy_richards"

        ODE_prob = ODEProblem(ODE_triple_piecewise_bertalanffy_richards, u0, tspan, param)

    end

    if model == "gbsm_piecewise"

        ODE_prob = ODEProblem(ODE_gbsm_piecewise, u0, tspan, param)

    end
    if model == "logistic"

        ODE_prob = ODEProblem(logistic, u0, tspan, param)

    end
    if model == "dlogistic"

        ODE_prob = ODEProblem(dlogistic, u0, tspan, param)

    end
    if model == "exponential"

        ODE_prob = ODEProblem(ODE_exponential, u0, tspan, param)

    end

    if model == "gompertz"

        ODE_prob = ODEProblem(gompertz, u0, tspan, param)

    end

    if model == "baranyi_richards"

        ODE_prob = ODEProblem(baranyi_richards, u0, tspan, param)

    end
    if model == "baranyi_exp"

        ODE_prob = ODEProblem(baranyi_exp, u0, tspan, param)

    end

    if model == "baranyi_roberts"

        ODE_prob = ODEProblem(baranyi_roberts, u0, tspan, param)

    end


    if model == "huang"
        ## attention!!!!!
        u0 = [log(u0)]



        function loss_ode_L2_relative_huang(p)

            sol = solve(ODE_prob, KenCarp4(autodiff=false), p=p, saveat=tsteps, verbose=false, abstol=1e-10, reltol=1e-10)
            sol_t = reduce(hcat, sol.u)
            lossa = 0.5 * NaNMath.sum(abs2.(1.0 .- (log.(data[2, :]) ./ sol_t[1, 1:end]))) / length(data[2, :])

            return lossa, sol
        end

        ODE_prob = ODEProblem(huang, u0, tspan, param)


    end

    if model == "piecewise_damped_logistic"


        ODE_prob = ODEProblem(ODE_piecewise_damped_logistic, u0, tspan, param)


    end



    if model == "Diauxic_replicator_1"


        ODE_prob = ODEProblem(ODE_Diauxic_replicator_1, u0, tspan, param)


    end
    if model == "Diauxic_replicator_2"


        ODE_prob = ODEProblem(ODE_Diauxic_replicator_2, u0, tspan, param)


    end
    if model == "Diauxic_piecewise_damped_logistic"


        ODE_prob = ODEProblem(ODE_Diauxic_piecewise_damped_logistic, u0, tspan, param)


    end
    if model == "triple_piecewise_damped_logistic"


        ODE_prob = ODEProblem(ODE_triple_piecewise_damped_logistic, u0, tspan, param)


    end


    if model == "triple_piecewise"


        ODE_prob = ODEProblem(ODE_triple_piecewise, u0, tspan, param)

    end
    if model == "triple_piecewise_sublinear"


        ODE_prob = ODEProblem(ODE_triple_piecewise_sub_linear, u0, tspan, param)

    end



    if model == "four_piecewise"


        ODE_prob = ODEProblem(ODE_four_piecewise, u0, tspan, param)

    end
  
    if model == "HPM_3_death"


        ODE_prob = ODEProblem(HPM_3_death, u0, tspan, param)

    end

    if model == "HPM_3_death_resistance"


        ODE_prob = ODEProblem(HPM_3_death_resistance, u0, tspan, param)

    end

    if model == "dHPM_3_death_resistance"


        ODE_prob = ODEProblem(dHPM_3_death_resistance, u0, tspan, param)

    end
      
    if model == "HPM_3_inhibition"


        ODE_prob = ODEProblem(HPM_3_inhibition, u0, tspan, param)

    end
    if model == "HPM"
            
        ODE_prob = ODEProblem(ODEs_McKellar, u0, tspan, param)
    end

    if model == "HPM_exp"
            
        ODE_prob = ODEProblem(ODEs_HPM_exp, u0, tspan, param)
    end

    if model == "dHPM"
        
        ODE_prob = ODEProblem(ODEs_damped_McKellar, u0, tspan, param)
    end
        # defining the loss functions for this case

    if model == "HPM_inhibition"
        
            ODE_prob = ODEProblem(ODEs_HPM_inhibition, u0, tspan, param)
    end
    if model == "dHPM_inhibition"
        
           ODE_prob = ODEProblem(ODEs_dHPM_inhibition, u0, tspan, param)
    end
    if model == "ODEs_HPM_SR"
        
        ODE_prob = ODEProblem(ODEs_HPM_SR, u0, tspan, param)
    end


    return ODE_prob

end  

function gaussian_smoothing(data::Matrix{Float64};optimize_gp=false)


    """
    Gaussian smoothing 
   """
    xtrain = data[1,:];
     ytrain = data[2,:];

     kernel = SE(4.0,4.0) +   RQ(0.0,0.0,-1.0) + SE(-2.0,-2.0);#GaussianProcesses.Periodic(0.0,1.0,0.0)*SE(4.0,0.0)
         gp = GP(xtrain,ytrain,MeanZero(),kernel,-2.0)   #Fit the GP
    
    if optimize_gp == true
        optimize!(gp) 
    end
    μ, Σ = predict_y(gp,range(data[1,1],stop=data[1,end],length=500));

    μ, Σ = predict_y(gp,xtrain);
    #a = range(data[1,1],stop=data[1,end],length=500)
   # range_vec = [ a[i] for i in 1:length(a)]

  return xtrain ,μ, Σ
end

function specific_gr_evaluation(data_smooted::Matrix{Float64},
     pt_smoothing_derivative::Int)


    """
    specific gr evaluation with slinding window log-lin fitting
    data_smooted = matrix of data
    pt_smoothing_derivative = size of the win, if <2 the the numerical derivative of (log) data is evaluate with interpolation algorithm
    """
    if pt_smoothing_derivative > 1
        specific_gr = [curve_fit(LinearFit, data_smooted[1, r:(r+pt_smoothing_derivative)], log.(data_smooted[2, r:(r+pt_smoothing_derivative)])).coefs[2] for r in 1:1:(length(data_smooted[2, :])-pt_smoothing_derivative)]

    else
        
        itp = interpolate((data_smooted[1,:],), log.(data_smooted[2,:]), Gridded(Linear()));
        specific_gr = only.(Interpolations.gradient.(Ref(itp), data_smooted[1,:]))
    end


    return specific_gr
end







# min max function hard box constrain on the parmeters

function max_min(min_param, max_param, param)
    return max(min_param, min(max_param, param))
end



# sigmoid function  box constrain on the parmeters

function si_ghc(min_param, max_param, param, slope)
    return max_param * (1 / (1 + exp(-slope * param))) + min_param
end




########################## No CONSTRAIN
# heterogeneous population model (HPM) with dormant active and starving population




function HPM_3_inhibition(du, u, param, t)
          
    du[1] = - u[1] * param[2]
          
    du[2] = u[1] * param[2] + param[1] * u[2] -  param[3]*  u[2]
    
    du[3] = param[3]*  u[2]# - param[4]*  u[3]

          
end

function HPM_3_death(du, u, param, t)
          
    du[1] = - u[1] * param[2]
          
    du[2] = u[1] * param[2] + param[1] * u[2] -  param[3]*  u[2]
    
    du[3] = param[3]*  u[2] - param[4]*  u[3]

          
end

function HPM_3_death_resistance_old(du, u, param, t)
          
    du[1] = - u[1] * param[2]
          
    du[2] = u[1] * param[2] + param[1] * u[2] * (1 - (u[3]+u[1]+u[2])/param[6]) -  param[3]*  u[2] 

    
    du[3] = param[3]*  u[2] + param[4]*  u[3]* (1 - (u[3]+u[1]+u[2])/param[5])

          
end
function HPM_3_death_resistance_old_old(du, u, param, t)
          
    du[1] = - u[1] * param[2]
          
    du[2] = u[1] * param[2] + param[1] * u[2] -  param[3]*  u[2] *(1-u[2])
  
    
    du[3] = param[3]* u[2]* (1-u[2]) + param[4]*  u[3]* (1 - (u[3]+u[1]+u[2])/param[5])
  
          
  end
function dHPM_3_death_resistance_old_2(du, u, param, t)
          
    du[1] = - u[1] * param[2]
          
    du[2] = u[1] * param[2] + param[1] * u[2] -  param[3]*  u[2]
    
    du[3] = param[3]*  u[2] + param[4]*  u[3]* (1 - NaNMath.pow((u[3]+u[1]+u[2])/param[5],param[6]))

          
end
function dHPM_3_death_resistance(du, u, param, t)
          
    du[1] = - u[1] * param[2]
          
    du[2] = u[1] * param[2] + param[1] * u[2]  -  param[3]*  u[2]#*(1-u[2])
    
    du[3] = param[3]*  u[2]  + param[4]*  u[3]* (1 - NaNMath.pow((u[3]+u[1]+u[2])/param[5],param[6]))

          
end

#function HPM_3(du, u, param, t)
          
#    du[1] = - u[1] * param[2]
          
#    du[2] = u[1] * param[2] + param[1] * u[2] -  param[3]* (1/(1+exp(-param[5]*t) ) ) * u[2]

#    du[3] = +  param[3]* (1/(1+exp(-param[5]*t) ) ) * u[2] - u[3]* param[4]

          
#end

#function HPM_3(du, u, param, t)
          
#    du[1] = - u[1] * param[2]
          
#    du[2] = u[1] * param[2] + param[1] * u[2]*(1-u[3]/param[5]) - param[3]* u[2]

#    du[3] = + param[3]* u[2] - u[3]* param[4]

          
#end


# McKellar (1997): heterogeneous population model (HPM).

function ODEs_McKellar(du, u, param, t)

    du[1] = - u[1] * (param[2]) 

    du[2] = u[1] * (param[2])  +  param[1] *  u[2] * ( 1 - (u[1] + u[2])/param[3]) 

end
function ODEs_HPM_exp(du, u, param, t)

    du[1] = - u[1] * (param[2]) 

    du[2] = u[1] * (param[2])  +  param[1] *  u[2] 

end
function ODEs_HPM_inhibition(du, u, param, t)

    du[1] =  - u[1] * (param[2])  +  param[1] *  u[1] 

    du[2] = + u[1] * (param[2]) + param[4] * u[2]* ( 1 - (u[1] + u[2])/param[3]) 

end

function ODEs_dHPM_inhibition(du, u, param, t)

    du[1] =  - u[1] * (param[2])  +  param[1] *  u[1] 

    du[2] = + u[1] * (param[2]) + param[4] * u[2]* ( 1 -  NaNMath.pow( (u[1] + u[2])/param[3],param[5])) 

end


function ODEs_HPM_SR(du, u, param, t)



    du[1] =  -  param[4]  / (1+ param[3] * exp(- param[2] * t) )*  u[1] +  param[1] *  u[1] -   param[5]  * u[1] 

    du[2] =  +  param[5]  * u[1] 

end

# damped heterogeneous population model (dHPM).

function ODEs_damped_McKellar(du, u, param, t)

    du[1] = - u[1] * (param[2]) 

    du[2] = u[1] * (param[2])  +  param[1] *  u[2] * ( 1 - NaNMath.pow( (u[1] + u[2])/param[3],param[4])) 

end

# Diauxic  replicator model 1 from "Diauxic behaviour for biological processes at various timescales"


function ODE_Diauxic_replicator_1(du, u, param, t)

    if t <= param[3]
        du[1] = param[5]

    else

        du[1] = u[1] * (param[2] - u[1]) * (NaNMath.pow(u[1] - param[1], 2) + param[4])

    end

end

# Diauxic  replicator model 2    "Diauxic behaviour for biological processes at various timescales"



function ODE_Diauxic_replicator_2(du, u, param, t)

    if t <= param[3]
        du[1] = param[5]

    else


        du[1] = u[1] * (param[2] - u[1]) * (NaNMath.pow(u[1] - param[1], 2) * NaNMath.pow(u[1] - param[6], 2) + param[4])



    end

end

# empirical Diauxic



function ODE_Diauxic_piecewise_damped_logistic(du, u, param, t)

    if (t <=  param[4] && t <= param[6] && t <= param[10])
        du[1] = u[1] * param[5]

    elseif (t > param[4] && t <= param[6] && t <= param[10])

        du[1] = param[1] * u[1] * (1 - NaNMath.pow((u[1] / param[2]), param[3]))

    elseif (t > param[4] && t > param[6] && t <= param[10] )

        du[1] = u[1] * param[11]  

    elseif (t > param[4] && t > param[6] && t> param[10] )
        du[1] = u[1] * param[7] * (1 - NaNMath.pow((u[1] / param[8]), param[9]))

    end


end



function ODE_Diauxic_piecewise_damped_logistic_bk(du, u, param, t)

    if (t <=  param[4] && t <= param[6] && t <= param[10])
        du[1] = u[1] * param[5]

    elseif (t > param[4] && t <= param[6] && t <= param[10])

        du[1] = param[1] * u[1] * (1 - NaNMath.pow((u[1] / param[2]), param[3]))
    elseif (t > param[4] && t > param[6] && t <= param[10] )

        du[1] = u[1] * param[11]

    elseif (t > param[4] && t > param[6] && t> param[10] )
        du[1] = u[1] * param[7] * (1 - NaNMath.pow((u[1] / param[8]), param[9]))

    end


end


# global optimizator piecewise 



# custom piecewise model
function ODE_exponential(du, u, param, t)

        du[1] = param[1]*u[1]


end
function ODE_piecewise_damped_logistic(du, u, param, t)

    if t <= param[3]
        du[1] = param[5]

    else

        du[1] = param[1] * u[1] * (1 - NaNMath.pow((u[1] / param[2]), param[4]))


    end

end
function ODE_triple_piecewise_damped_logistic(du, u, param, t)

    if t <= param[3]
        du[1] = u[1] * param[5]

    elseif ( t <=  param[6] && t > param[3] )


        du[1] = param[1] * u[1] * (1 - NaNMath.pow((u[1] / param[2]), param[4]))


    elseif  ( t >  param[6] && t > param[3] ) 
        du[1] = u[1] * param[7]
    end

end


function ODE_triple_piecewise(du, u, param, t)

    if t <= param[4]

        du[1] = u[1] * param[2]

    elseif (t <= param[5] && t > param[4] )

        du[1] = u[1] * param[1]


    elseif (t > param[5] && t > param[4] )

        du[1] = u[1] * param[3]
    end

end


function ODE_triple_piecewise_sub_linear(du, u, param, t)

    if t <= param[4]

        du[1] = u[1] * param[2]

    elseif t <= param[5]

        du[1] = u[1] * param[1]


    else

        du[1] =(1- NaNMath.log( u[1]/param[6] ))  * param[3]
    end

end


function ODE_gbsm_piecewise(du, u, param, t)

    if t <= param[4]

        du[1] = param[1] * (1 / (1 + NaNMath.pow(abs((t - param[4]) / param[2]), 2 * param[3])))

    else

        du[1] = param[1] * (1 / (1 + NaNMath.pow(abs((t - param[4]) / param[5]), 2 * param[6])))


    end

end

function ODE_four_piecewise(du, u, param, t)

    if t <= param[5]

        du[1] = u[1] * param[2]

    elseif (t <= param[6] && t > param[5])

        du[1] = u[1] * param[1]

    elseif (t <= param[7] && t >  param[6] && t >  param[5])

        du[1] = u[1] * param[3]

    elseif (t > param[7] && t >  param[6] && t >  param[5])

        du[1] = u[1] * param[4]
    end

end

# hyper gompertz curve from "A Theory of Growth" Turner, Brandley and Kirk 1976

function hyper_gompertz(du, u, param, t)
    du[1] = param[1] * u[1] * NaNMath.pow(NaNMath.log(max(param[2], u[1]) / u[1]), 1 + param[3])
end

# hyper logistic curve from "A Theory of Growth" Turner, Brandley and Kirk 1976

function hyper_logistic(du, u, param, t)
    du[1] = (param[1] / param[2]) * NaNMath.pow(u[1], (1 - param[3])) * NaNMath.pow(max(param[2], u[1] + 0.00001) - u[1], 1 + param[3])

end


# Bertalanffy Richards curve from "A Theory of Growth" Turner, Brandley and Kirk 1976

function bertalanffy_richards(du, u, param, t)
    du[1] = NaNMath.pow(param[1] / param[2], param[3]) * u[1] * (NaNMath.pow(max(param[2], u[1] + 0.00001), param[3]) - NaNMath.pow(u[1], param[3]))

end
# Multiplicative modelling of four-phase microbial growth, ODE von Bertalanffy

function ODE_von_bertalanffy(du, u, param, t)
    du[1] = u[1] * ( param[1] * param[2] * t ^(param[2] - 1 ) - param[3] * param[4] * t ^(param[4] - 1 ))

end

# third part from A general model for ontogenetic growth
function ODE_triple_piecewise_bertalanffy_richards(du, u, param, t)

    if t <= param[3]

        du[1] = u[1] * param[2]

    elseif t <= param[4]

        du[1] = u[1] * param[1]


    else

        du[1] = param[5]*  NaNMath.pow(u[1],param[6]) * ( 1 - NaNMath.pow(u[1]/param[7],param[6]) )
    end

end


# Logistic curve 

function logistic(du, u, param, t)
    du[1] = (param[1] / param[2]) * u[1] * (param[2] - u[1])

end
# dLogistic curve 

function dlogistic(du, u, param, t)
    du[1] = (param[1] ) * u[1] * ( 1 - NaNMath.pow(u[1]/param[2],param[3]) )

end

# Gompertz  curve 


function gompertz(du, u, param, t)
    du[1] = (param[1]) * u[1] * log(param[2] / u[1])

end

# Baranyi-Richards model from Baranyi, Roberts, and   McClure. "A non-autonomous differential equation to model bacterial growth" 1993


function baranyi_richards(du, u, param, t)
    du[1] = param[1] * (1 - (u[1] / param[2])) * (t^param[4]) / ((param[3])^(param[4]) + t^(param[4])) * u[1]
end
function baranyi_exp(du, u, param, t)
    du[1] = param[1] *  (t^param[3]) / ((param[2])^(param[3]) + t^(param[3])) * u[1]
end
#  Baranyi-Roberts model from: Baranyi and Roberts. "A dynamic approach to predicting bacterial growth in food" 1994

function baranyi_roberts(du, u, param, t)
    du[1] = param[1] * (1 - NaNMath.pow(u[1] / max(param[2], u[1] + 0.00001), param[5])) * (NaNMath.pow(t, param[4]) / (NaNMath.pow(param[3], param[4]) + NaNMath.pow(t, param[4]))) * u[1]
end

# Huang model from Huang "Optimization of a new mathematical model for bacterial growth" 2013

function huang(du, u, param, t)
    du[1] = param[1] * (1 - exp(u[1] - param[2])) / (1 + exp(4.0 * (t - param[3])))
end





#######################################################################


# function to do the smoothing average of data
function smoothing_data(data::Matrix{Float64},
    pt_avg::Int)
    """
    rolling average smoothing of the data
    data  = matrix of data
    pt_avg = size of the windows of the rolling average
    """
    times = [sum(@view data[1, i:(i+pt_avg-1)]) / pt_avg for i in 1:(length(data[1, :])-(pt_avg-1))]
    values = [sum(@view data[2, i:(i+pt_avg-1)]) / pt_avg for i in 1:(length(data[2, :])-(pt_avg-1))]
    smoothed_data = transpose(hcat(times, values))
    return smoothed_data

end






function vectorize_df_results(well_name::String,
    model::String,
    res::Any,
    th_gr::Any,
    em_gr::Any,
    loss::Float64)
    """
        internal function to reder as vector the results of the hardcoded ODE 
    """


    

    if model == "HPM_3_inhibition"

        res_param = [model,
            well_name,
            res[1], # growth rate
            res[2], # rate exit lag
            res[3], # inhibition rate
            th_gr,
            em_gr,
            loss]  # error function 
    end
    if model == "HPM_3_death"

        res_param = [model,
            well_name,
            res[1], # growth rate
            res[2], # rate exit lag
            res[3], # inhibitionrate
            res[4], # death rate
            th_gr,
            em_gr,
            loss]  # error function 
    end

    if model == "HPM_3_death_resistance"

        res_param = [model,
            well_name,
            res[1], # growth rate
            res[2], # rate exit lag
            res[3], # inhibitionrate
            res[4], # death rate
            res[5], # n_res
            res[6], # n_max
            th_gr,
            em_gr,
            loss]  # error function 
    end

    
    if model == "dHPM_3_death_resistance"

        res_param = [model,
            well_name,
            res[1], # growth rate
            res[2], # rate exit lag
            res[3], # inhibitionrate
            res[4], # death rate
            res[5], # n_res
            res[6], # shape
            th_gr,
            em_gr,
            loss]  # error function 
    end

    if model == "HPM"

        res_param = [model,
            well_name,
            res[1], # growth rate
            res[2], # rate exit lag
            res[3], # N max
            th_gr,
            em_gr,
            loss]  # error function 
    end
    if model == "HPM_exp"

        res_param = [model,
            well_name,
            res[1], # growth rate
            res[2], # rate exit lag
            th_gr,
            em_gr,
            loss]  # error function 
    end

    if model == "HPM_inhibition"

        res_param = [model,
            well_name,
            res[1], # growth rate
            res[2], # rate inhibition
            res[3],# growth inibithed
            res[4], # N max
            th_gr,
            em_gr,
            loss]  # error function 
    end

    if model == "dHPM_inhibition"

        res_param = [model,
            well_name,
            res[1], # growth rate
            res[2], # rate inhibiti
            res[3],# growth inibithed
            res[4], # N max
            res[5], # shape
            th_gr,
            em_gr,
            loss]  # error function 
    end
    if model == "ODEs_HPM_SR"

        res_param = [model,
            well_name,
            res[1], # growth rate
            res[2], # growth of the phages
            res[3], # maximum poin of the growth of phages
            res[4], # maximum death rate
            res[5], # resistance rate
            th_gr,
            em_gr,
            loss]  # error function 
    end

    if model == "dHPM"

        res_param = [model,
            well_name,
            res[1], # growth rate
            res[2], # rate exit lag
            res[3], # N max
            res[4], # shape
            th_gr,
            em_gr,
            loss]  # error function 
    end


    if model == "hyper_gompertz"

        res_param = [model,
            well_name,
            res[1], # growth rate
            res[2], # carrying capacity
            res[3], # shape factor
            th_gr,
            em_gr,
            loss]  # error function 
    end

    if model == "hyper_logistic"

        res_param = [model,
            well_name,
            res[1], # growth rate
            res[2], # carrying capacity
            res[3], # shape factor
            th_gr,
            em_gr,
            loss]  # error function 


    end


    if model == "gbsm_piecewise"

        res_param = [model,
            well_name,
            res[1], # growth rate
            res[2], # a_1
            res[3], # b_1
            res[4], # c
            res[5], # a_1
            res[6], # b_1
            th_gr,
            em_gr,
            loss]  # error function 

    end

    if model == "bertalanffy_richards"

        res_param = [model,
            well_name,
            res[1], # growth rate
            res[2], # carrying capacity
            res[3], # shape factor
            th_gr,
            em_gr,
            loss]  # error function 

    end


    if model == "ode_von_bertalanffy"

        res_param = [model,
            well_name,
            res[1], # alpha
            res[2], # beta
            res[3], # a
            res[4], # b
            th_gr,
            em_gr,
            loss]  # error function 

    end

    if model == "triple_piecewise_bertalanffy_richards"

        res_param = [model,
            well_name,
            res[1], # growth rate
            res[2], # lag growth
            res[3], # lag time
            res[4], # stat time 
            res[5], # growth stat
            res[6], # shape factor
            res[7], # n max
            th_gr,
            em_gr,
            loss]  # error function 

    end
    if model == "exponential"

        res_param = [model,
            well_name,
            res[1], # growth rate
            th_gr,
            em_gr,
            loss]  # error function 

    end


    if model == "logistic"

        res_param = [model,
            well_name,
            res[1], # growth rate
            res[2], # carrying capacity
            th_gr,
            em_gr,
            loss]  # error function 

    end
    if model == "dlogistic"

        res_param = [model,
            well_name,
            res[1], # growth rate
            res[2], # carrying capacity
            res[3], # shape
            th_gr,
            em_gr,
            loss]  # error function 

    end
    if model == "gompertz"

        res_param = [model,
            well_name,
            res[1], # growth rate
            res[2], # carrying capacity
            th_gr,
            em_gr,
            loss]  # error function 
    end




    if model == "piecewise_damped_logistic"

        res_param = [model,
            well_name,
            res[1], # growth rate
            res[2], # carrying capacity
            res[3], # lag time
            res[4], # shape factor 
            res[5], # linear factor 
            th_gr,
            em_gr,
            loss]  # error function 


    end

    if model == "triple_piecewise_damped_logistic"

        res_param = [model,
            well_name,
            res[1], # growth rate
            res[2], # carrying capacity
            res[3], # lag time
            res[4], # shape factor 
            res[5], # linear factor 
            res[6], # stationary start
            res[7], # linear param of stationary phase
            th_gr,
            em_gr,
            loss]  # error function 


    end

    if model == "triple_piecewise"

        res_param = [model,
            well_name,
            res[1], # growth rate 1
            res[2], # growth rate 2
            res[3], # growth rate 3
            res[4], # lag
            res[5], # stationary start      
            th_gr,
            em_gr,
            loss]  # error function 


    end
    if model == "triple_piecewise_sublinear"

        res_param = [model,
            well_name,
            res[1], # growth rate 1
            res[2], # growth rate 2
            res[3], # growth rate 3
            res[4], # lag
            res[5], # stationary start
            res[6], # carrying        
            th_gr,
            em_gr,   
            loss]  # error function 


    end


    if model == "four_piecewise"

        res_param = [model,
            well_name,
            res[1], # growth rate 1
            res[2], # growth rate 2
            res[3], # growth rate 3
            res[4], # growth rate 4
            res[5], # lag
            res[6], # decrease of gr start  
            res[7], # stationary start  
            th_gr,
            em_gr,
            loss]  # error function 


    end

    if model == "triple_piecewise"

        res_param = [model,
            well_name,
            res[1], # growth rate 1
            res[2], # growth rate 2
            res[3], # growth rate 3
            res[4], # lag
            res[5], # stationary start      
            th_gr,
            em_gr,
            loss]  # error function 


    end
    if model == "baranyi_richards"

        res_param = [model,
            well_name,
            res[1], # growth rate
            res[2], # carrying capacity
            res[3], # lag time
            res[4], # shape factor
            th_gr,
            em_gr,
            loss]  # error function 

    end
    if model == "baranyi_exp"

        res_param = [model,
            well_name,
            res[1], # growth rate
            res[2], # lag time
            res[3], # shape factor
            th_gr,
            em_gr,
            loss]  # error function 

    end
    if model == "baranyi_roberts"

        res_param = [model,
            well_name,
            res[1], # growth rate
            res[2], # carrying capacity
            res[3], # lag time
            res[4], # shape factor 1
            res[5], # shape factor 1
            th_gr,
            em_gr,
            loss]  # error function 

    end


    if model == "huang"

        res_param = [model,
            well_name,
            res[1], # growth rate
            res[2], # carrying capacity
            res[3], # lag time
            th_gr,
            em_gr,
            loss]  # error function 

    end



    if model == "Diauxic_replicator_1"

        res_param = [model,
            well_name,
            res[1], # growth rate
            res[2], # carrying capacity
            res[3], # lag time
            res[4], # arbritrary small constrain
            res[5], #growth stationary phase
            th_gr,
            em_gr,
            loss]  # error function 

    end



    if model == "Diauxic_replicator_2"

        res_param = [model,
            well_name,
            res[1], # growth rate
            res[2], # carrying capacity
            res[3], # lag time
            res[4], # arbritrary small constrain
            res[5], #growth stationary phase
            res[6], #growth second phase
            th_gr,
            em_gr,
            loss]  # error function 

    end



    if model == "Diauxic_piecewise_damped_logistic"
        res_param = [model,
            well_name,
            res[1], # growth rate
            res[2], # carrying capacity
            res[3], # shape
            res[4], # lag time 
            res[5], # linear factor 
            res[6], # diaux shift 
            res[7], # growth second part
            res[8], # carrying capacity 2
            res[9],#  shape 2
            res[10],# end second lag
            res[11],# gr second lag
            th_gr,
            em_gr,
            loss]  # error function 

    end

    return res_param

end







function inzialize_df_results(model::String)

    """
        internal function to name the vectors the results of the hardcoded ODE 
    """


    if model == "HPM_3_death_resistance"

        param_names = ["model", "well", "gr", "exit_lag_rate", "inactivation_rate","death_rate","n_res","n_max","th_max_gr","emp_max_gr","loss"]

    end

    if model == "dHPM_3_death_resistance"

        param_names = ["model", "well", "gr", "exit_lag_rate", "inactivation_rate","death_rate","n_res","shape","th_max_gr","emp_max_gr","loss"]

    end
    if model == "HPM_3_inhibition"

        param_names = ["model", "well", "gr", "exit_lag_rate", "inactivation_rate", "th_max_gr","emp_max_gr","loss"]

    end

    if model == "HPM_3_death"

        param_names = ["model", "well", "gr", "exit_lag_rate", "inactivation_rate","death_rate","th_max_gr","emp_max_gr","loss"]
    
    end

    if model == "HPM"

        param_names = ["model", "well", "gr", "exit_lag_rate", "N_max","th_max_gr","emp_max_gr","loss"]

    end
    if model == "HPM_exp"

        param_names = ["model", "well", "gr", "exit_lag_rate","th_max_gr","emp_max_gr","loss"]

    end
    if model == "ODEs_HPM_SR"

        param_names = ["model", "well", "gr", "gr_phage","scale", "death_rate","resistance_rate","th_max_gr","emp_max_gr", "loss"]

    end

    if model == "HPM_inhibition"

        param_names = ["model", "well", "gr", "inhibition_rate", "gr_inhibition","N_max","th_max_gr","emp_max_gr", "loss"]

    end


    if model == "dHPM"

        param_names = ["model", "well", "gr", "exit_lag_rate", "N_max", "shape" ,"th_max_gr","emp_max_gr", "loss"]

    end
    if model == "dHPM_inhibition"

        param_names = ["model", "well", "gr", "inhibition_rate", "gr_inhibition","N_max", "shape","th_max_gr","emp_max_gr" , "loss"]

    end
    if model == "hyper_gompertz"

        param_names = ["model", "well", "gr", "N_max", "shape","th_max_gr","emp_max_gr", "loss"]

    end

    if model == "hyper_logistic"
        param_names = ["model", "well", "doubling_time", "gr", "N_max", "shape","th_max_gr","emp_max_gr", "loss"]


    end


    if model == "gbsm_piecewise"
        param_names = ["model", "well", "gr", "a_1", "b_1", "c", "a_2", "b_2","th_max_gr","emp_max_gr", "loss"]


    end

    if model == "bertalanffy_richards"
        param_names = ["model", "well", "gr", "N_max", "shape","th_max_gr","emp_max_gr", "loss"]


    end

    
    if model == "ode_von_bertalanffy"
        param_names = ["model", "well", "alpha", "beta", "a","b","th_max_gr","emp_max_gr", "loss"]


    end

    if model == "triple_piecewise_bertalanffy_richards"
        param_names = ["model", "well", "gr", "gr_lag","t_lag","t_stationary","gr_stat", "shape", "N_max","th_max_gr","emp_max_gr","loss"]


    end

    if model == "logistic"
        param_names = ["model", "well", "gr", "N_max","th_max_gr","emp_max_gr", "loss"]


    end
    if model == "logistic"
        param_names = ["model", "well", "gr", "N_max","shape","th_max_gr","emp_max_gr", "loss"]


    end

    if model == "exponential"
        param_names = ["model", "well", "gr", "th_max_gr","emp_max_gr", "loss"]


    end

    if model == "gompertz"
        param_names = ["model", "well", "gr", "N_max","th_max_gr","emp_max_gr", "loss"]



    end

    if model == "baranyi_richards"
        param_names = ["model", "well", "gr", "N_max", "lag_time", "shape","th_max_gr","emp_max_gr", "loss"]


    end
    if model == "baranyi_exp"
        param_names = ["model", "well", "gr", "lag_time", "shape","th_max_gr","emp_max_gr", "loss"]


    end
    if model == "baranyi_roberts"
        param_names = ["model", "well", "gr", "N_max", "lag_time", "shape_1", "shape_2","th_max_gr","emp_max_gr", "loss"]

    end


    if model == "huang"
        param_names = ["model", "well", "gr", "N_max", "lag","th_max_gr","emp_max_gr", "loss"]



    end


    if model == "piecewise_damped_logistic"

        param_names = ["model", "well", "gr", "N_max", "lag", "shape", "linear_const","th_max_gr","emp_max_gr", "loss"]

    end

    if model == "triple_piecewise_damped_logistic"

        param_names = ["model", "well", "gr", "N_max", "lag", "shape", "linear_const", "t_stationary", "linear_lag","th_max_gr","emp_max_gr", "loss"]

    end

    if model == "triple_piecewise"

        param_names = ["model", "well", "gr", "gr_2", "gr_3", "lag", "t_stationary","th_max_gr","emp_max_gr", "loss"]

    end

    if model == "triple_piecewise_sublinear"

        param_names = ["model", "well", "gr", "gr_2", "gr_3", "lag", "t_stationary","N_max","th_max_gr","emp_max_gr", "loss"]

    end
    if model == "four_piecewise"

        param_names = ["model", "well", "gr", "gr_2", "gr_3", "gr_4", "lag", "t_decay_gr", "t_stationary","th_max_gr","emp_max_gr", "loss"]

    end

    if model == "Diauxic_replicator_1"

        param_names = ["model", "well", "gr", "N_max", "lag", "arbitrary_const", "linear_const","th_max_gr","emp_max_gr", "loss"]

    end

    if model == "Diauxic_replicator_2"

        param_names = ["model", "well", "gr", "N_max", "lag", "arbitrary_const", "linear_const", "growth_stationary","th_max_gr","emp_max_gr", "loss"]

    end

    if model == "Diauxic_piecewise_damped_logistic"

        param_names = ["model", "well", "gr_1", "N_max", "shape_1", "lag", "linear_const", "t_shift", "gr_2","N_max_2","shape_2","end_second_lag","lag_2_gr","th_max_gr","emp_max_gr", "loss"]

    end




    return param_names

end

function guess_param(
    lb_param::Vector{Float64},
    ub_param::Vector{Float64}
)
"""
internal function to set the start of the optimization problem in the middle of the  box constrains
"""
  param = lb_param .+ (ub_param-lb_param)./2

return param

end

######################################################################

# negative blank correction

function thr_negative_correction(data::Matrix{Float64}, # dataset x times y OD/fluorescence ,
    thr_negative::Float64 # the value at which data are setted if negative
    )
    times = data[1, :]
    values = data[2, :]
    index_neg = findall(x -> x < thr_negative, values)
    values[index_neg] .= thr_negative
    data_corrected = transpose(hcat(times, values))

    return data_corrected

end

function blank_distrib_negative_correction(data::Matrix{Float64}, # dataset x times y OD/fluorescence ,
    blank_array::Vector{Float64}
)

    times = data[1, :]
    values = data[2, :]
    #println(values)
    index_neg = findall(x -> x < 0.0, values)

    number_of_neg = length(index_neg)
    replacement_values = abs.(sample(blank_array, number_of_neg) .- mean(blank_array))
    values[index_neg] .= replacement_values
    data_corrected = transpose(hcat(times, values))

    return data_corrected


end


####################################################################

#######

function fitting_one_well_Log_Lin_Fit(data::Matrix{Float64}, # dataset x times y OD/fluorescence
    name_well::String, # name of the well
    label_exp::String; #label of the experiment
    do_plot=false, # do plots or no
    path_to_plot='', # where save plots
    type_of_smoothing="rolling_avg", # option, NO, gaussian, rolling avg
    pt_avg=7, # number of the point for rolling avg not used in the other cases
    pt_smoothing_derivative=7, # number of poits to smooth the derivative
    pt_min_size_of_win=7, # minimum size of the exp windows in number of smooted points
    type_of_win="maximum", # how the exp. phase win is selected, "maximum" of "global_thr"
    threshold_of_exp=0.9, # threshold of growth rate in quantile to define the exp windows
    multiple_scattering_correction=false, # if true uses the given calibration curve to fix the data
    calibration_OD_curve ="" #  the path to calibration curve to fix the data
    )


    if  multiple_scattering_correction == true

        data = correction_OD_multiple_scattering(data,calibration_OD_curve)
    
    
    end
    if type_of_smoothing == "rolling_avg"

        data_smooted = smoothing_data(data, pt_avg)

    elseif  type_of_smoothing =="gaussian"

         temp =  gaussian_smoothing(data)
        data_smooted = Matrix(transpose(hcat(temp[1], temp[2])))

    else 

        data_smooted = data

    end

    # local fitting generating of specific growth rate (derivative)

    specific_gr = specific_gr_evaluation(data_smooted, pt_smoothing_derivative)

    specific_gr_times = [(data_smooted[1, r] + data_smooted[1, (r+pt_smoothing_derivative)]) / 2 for r in 1:1:(length(data_smooted[2, :])-pt_smoothing_derivative)]
    # selecting the max

    gr_max = maximum(specific_gr)
    index_of_max = findall(x -> x == gr_max, specific_gr)[1]

    # selecting the windows
    # lower threshold of the exp phase using quantile
    lb_of_distib = quantile(specific_gr, threshold_of_exp)

    # searching for t_start_exp 
    t_start = 0.0

    if type_of_win == "maximum"
        for yy in 1:(index_of_max-2)
            if specific_gr[index_of_max-yy] >= lb_of_distib
                t_start = copy(specific_gr_times[index_of_max-yy])
            else

                if specific_gr[(index_of_max-yy-1)] < lb_of_distib
                    break
                end

            end

        end

        # index of start
        index_of_t_start = findfirst(x -> x > t_start, data_smooted[1, :])[1]


        # searching t_end of the exp phase

        t_end = specific_gr_times[end]

        for yy in index_of_max:(length(specific_gr)-1)

            if specific_gr[yy] >= lb_of_distib
                t_end = copy(specific_gr_times[yy])
            else
                if specific_gr[(yy+1)] < lb_of_distib
                    break
                end
            end


        end



        index_of_t_end = findfirst(x -> x > t_end, data_smooted[1, :])[1]
    end

    # selection of exp win with a global thr on the growht rate
    if type_of_win == "global_thr"
        index_of_max = findfirst(x -> x == maximum(specific_gr), specific_gr)[1]
        index_gr_max = index_of_max + findfirst(x -> x < lb_of_distib, specific_gr[index_of_max:end])[1]
        index_gr_min = findlast(x -> x > lb_of_distib, specific_gr[1:index_of_max])[1]

        t_start = specific_gr_times[index_gr_min]

        t_end = specific_gr_times[index_gr_max]


        index_of_t_start = findfirst(x -> x > t_start, data_smooted[1, :])[1]
        index_of_t_end = findall(x -> x > t_end, data_smooted[1, :])[1]

    end
    # checking the minimum size of the window before fitting
    
    if (index_of_t_end - index_of_t_start) < pt_min_size_of_win

        index_of_t_start = convert(Int, index_of_max - floor(pt_min_size_of_win / 2))
        index_of_t_end = convert(Int, index_of_max + floor(pt_min_size_of_win / 2))

        if index_of_t_start < 1
            index_of_t_start = 2
        end   
        if index_of_t_end > length(data_smooted[1,:])
            index_of_t_end = length(data_smooted[1,:]) -1
        end

    end


    # fitting data

    data_to_fit_times = data_smooted[1, index_of_t_start:index_of_t_end]
    data_to_fit_values = log.(data_smooted[2, index_of_t_start:index_of_t_end])

    fitting_results = curve_fit(LinearFit, data_to_fit_times, data_to_fit_values)


    # residual calculation
    residual = [ ((data_to_fit_values[ll] -    fitting_results.coefs[2] *data_to_fit_times[ll] -  fitting_results.coefs[1])^2) for ll in 1 : length(data_to_fit_values)   ] 
    # sum_of_squares calculation
    sum_of_squares =  [ ((data_to_fit_values[ll] -   mean(data_to_fit_values) )^2) for ll in 1 : length(data_to_fit_values)   ]
    # coeff of determination 
    rsquared = 1 - sum(residual) / sum(sum_of_squares)
    # confidence interval growth rate
    a = ((1/(length(residual)-2))*sum(residual))
    std_error = sqrt( a/sum(sum_of_squares  ))
    p95 = ccdf(TDist(length(sum_of_squares)-2),(1-0.05/2))
    confidence_coeff_2 = std_error * p95
    # confidence interval intercept
    a_1 = sqrt( (1/(length(residual))) *  sum(data_to_fit_values.^2))
    std_error_intercept = std_error * a_1
    confidence_coeff_1 = std_error_intercept * p95
      

    ###
    # EVALUATING CONFIDENCE BANDS
    term_1 =  ((1/(length(residual)-2))*sum(residual))
    term_2 = (1/(length(residual))) .+ ( data_to_fit_times .-mean(data_to_fit_values))./(sum(sum_of_squares))
    confidence_band = p95 .* sqrt.( term_1 .*  term_2)
    
    fitted_line = [fitting_results.coefs[2] * data_to_fit_times[ll] + fitting_results.coefs[1] for ll in 1:length(data_to_fit_times)]
    # storing results

    results_lin_log_fit = [label_exp, name_well, data_to_fit_times[1], data_to_fit_times[end], specific_gr_times[index_of_max],gr_max ,fitting_results.coefs[2], confidence_coeff_2, log(2)/( fitting_results.coefs[2]) ,log(2)/( fitting_results.coefs[2] -confidence_coeff_2 ) , log(2)/( fitting_results.coefs[2] + confidence_coeff_2 ) ,fitting_results.coefs[1],confidence_coeff_1,rsquared]

    # plotting if requested
    if do_plot == true

        mkpath(path_to_plot)
        display(Plots.scatter(data_smooted[1, :], log.(data_smooted[2, :]), xlabel="Time", ylabel="Log(Arb. Units)", label=["Data " nothing], markersize=1, color=:black, title=string(label_exp, " ", name_well)))
        display(Plots.plot!(data_to_fit_times, fitted_line,ribbon= confidence_band, xlabel="Time ", ylabel="Log(Arb. Units)", label=[string("Fitting Log-Lin ") nothing], c=:red))
        display(Plots.vline!([data_to_fit_times[1], data_to_fit_times[end]], c=:black, label=[string("Window of exp. phase ") nothing]))
        png(string(path_to_plot, label_exp, "_Log_Lin_Fit_", name_well, ".png"))

        display(Plots.scatter(specific_gr_times, specific_gr, xlabel="Time ", ylabel="1 /time ", label=[string("Dynamics growth rate ") nothing], c=:red))
        display(Plots.vline!([data_to_fit_times[1], data_to_fit_times[end]], c=:black, label=[string("Window of exp. phase ") nothing]))
        png(string(path_to_plot, label_exp, "_dynamics_gr_", name_well, ".png"))



    end


    return results_lin_log_fit

end

######################################################################


function fitting_one_well_Log_Lin_Fit_for_errors(data::Matrix{Float64}, # dataset x times y OD/fluorescence
    pt_avg::Int, # numebr of the point to generate intial condition
    pt_smoothing_derivative::Int, # number of poits to smooth the derivative
    pt_min_size_of_win::Int, # minimum size of the exp windows in number of smooted points
    threshold_of_exp::Float64, # threshold of growth rate in quantile to define the exp windows
    correction::Bool
    )


    if  correction == true

        data = correction_OD_multiple_scattering(data,calibration_OD_curve)
    
    
    end


    data_smooted = smoothing_data(data, pt_avg)

    # local fitting generating of specific growth rate (derivative)

    specific_gr =specific_gr_evaluation(data_smooted, pt_smoothing_derivative)

    specific_gr_times = [(data_smooted[1, r] + data_smooted[1, (r+pt_smoothing_derivative)]) / 2 for r in 1:1:(length(data_smooted[2, :])-pt_smoothing_derivative)]
    # selecting the max

    gr_max = maximum(specific_gr)
    index_of_max = findall(x -> x == gr_max, specific_gr)[1]

    # selecting the windows
    # lower threshold of the exp phase using quantile
    lb_of_distib = quantile(specific_gr, threshold_of_exp)

    # searching for t_start_exp 
    t_start = 0.0

    for yy in 1:(index_of_max-2)

        if specific_gr[index_of_max-yy] >= lb_of_distib
            t_start = copy(specific_gr_times[index_of_max-yy])
        else

            if specific_gr[(index_of_max-yy-1)] < lb_of_distib
                break
            end

        end

    end

    # index of start
    index_of_t_start = findfirst(x -> x > t_start, data_smooted[1, :])[1]


    # searching t_end of the exp phase

    t_end = specific_gr_times[end]

    for yy in index_of_max:(length(specific_gr)-1)

        if specific_gr[yy] >= lb_of_distib
            t_end = copy(specific_gr_times[yy])
        else
            if specific_gr[(yy+1)] < lb_of_distib
                break
            end
        end


    end

    #  

    index_of_t_end = findfirst(x -> x > t_end, data_smooted[1, :])[1]


    # checking the minimum size of the window before fitting


    if (index_of_t_end - index_of_t_start) < pt_min_size_of_win

        index_of_t_start = convert(Int, index_of_max - floor(pt_min_size_of_win / 2))
        index_of_t_end = convert(Int, index_of_max + floor(pt_min_size_of_win / 2))

    end


    # fitting data

    data_to_fit_times = data_smooted[1, index_of_t_start:index_of_t_end]
    data_to_fit_values = log.(data_smooted[2, index_of_t_start:index_of_t_end])

    fitting_results = curve_fit(LinearFit, data_to_fit_times, data_to_fit_values)



    # returning array of results for this well 

    return fitting_results.coefs[2], data_to_fit_times[1], data_to_fit_times[end], specific_gr, specific_gr_times

end



function fitting_one_well_constrained(data::Matrix{Float64}, # dataset x times y OD/fluorescence
    name_well::String, # name of the well
    label_exp::String, #label of the experiment
    model::String, # ode model to use 
    lb_param::Vector{Float64}, # lower bound param
    ub_param::Vector{Float64}; # upper bound param
    param= lb_param .+ (ub_param.-lb_param)./2# initial guess param
    optmizator =   BBO_adaptive_de_rand_1_bin_radiuslimited(), # selection of optimization method 
    do_plot=false, # do plots or no
    path_to_plot="NA", # where save plots
    pt_avg=1, # numebr of the point to generate intial condition
    pt_smooth_derivative::Int,
    smoothing=false, # the smoothing is done or not?
    type_of_loss="RE", # type of used loss 
    blank_array=zeros(100), # data of all blanks
    error_analysis=false, # do error analysis
    write_res=false,
    multiple_scattering_correction=false, # if true uses the given calibration curve to fix the data
    calibration_OD_curve="NA"  #  the path to calibration curve to fix the data
    )


    if  multiple_scattering_correction == true

        data = correction_OD_multiple_scattering(data,calibration_OD_curve)
    
    
    end





    max_t = data[1, end]
    min_t = data[1, 1]
    tspan = (min_t, max_t)
    tsteps = data[1, :]


    # "starting condition  using data if smoothing average is used skip this part

    if smoothing == true

        u0 = [data[2, 1]]

    else
        u0 = [Statistics.mean(data[2, 1:pt_avg])]
    end

    #  definition the possible loss (not optimal here but its working)



    if model ∉ ["HPM", "HPM_exp","dHPM","HPM_inhibition","dHPM_inhibition","ODEs_HPM_SR","HPM_3_death","HPM_3_inhibition","HPM_3_death_resistance","dHPM_3_death_resistance"]
              
        ODE_prob =model_selector(model,u0,tspan)

        function loss_ode_L2(p)
            sol = solve(ODE_prob, KenCarp4(autodiff=true), p=p, saveat=tsteps, verbose=false, abstol=1e-10, reltol=1e-10)
            sol_t = reduce(hcat, sol.u)
            if size(sol_t)[2] == size(data)[2]
                lossa = NaNMath.sum(abs2.((data[2, :] .- sol_t[1, 1:end]))) / length(data[2, :])
            else
                lossa = 10.0^9 * length(data[2, :])

            # lossa = NaNMath.sum(abs2.((data[2, :] .- sol_t[1, 1:end]))) / length(data[2, :]) + 1000 * (size(data)[2] - size(sol_t)[2])
            end

            return lossa, sol
        end


        function loss_ode_L2_derivative(p)
            sol = solve(ODE_prob, KenCarp4(autodiff=true), p=p, saveat=tsteps, verbose=false, abstol=1e-10, reltol=1e-10)
            sol_t = reduce(hcat, sol.u)
            sol_time = reduce(hcat, sol.t)

            derivative_theo = [(sol_t[1, i] - sol_t[1, i-1]) / (sol_time[1, i] - sol_time[1, i-1]) for i in 2:(size(sol_t)[2])]

            derivative_data = [(data[2, i] - data[2, i-1]) / (data[1, i] - data[1, i-1]) for i in 2:(size(data)[2])]

            if size(derivative_theo)[1] == size(derivative_data)[1]
                lossa = NaNMath.sum(abs2.((derivative_theo[2, :] .- derivative_data[1, 1:end]))) / length(derivative_data[2, :])
            else
                lossa = 10.0^9 * length(data[2, :])
               # lossa = NaNMath.sum(abs2.((derivative_theo[2, :] .- derivative_data[1, 1:end]))) / length(derivative_data[2, :]) + 1000
            end

            return lossa, sol
        end


        function loss_ode_RE(p)
            sol = solve(ODE_prob, KenCarp4(autodiff=true), p=p, saveat=tsteps, verbose=false, abstol=1e-10, reltol=1e-10)
            sol_t = reduce(hcat, sol.u)
            if size(sol_t)[2] == size(data)[2]
                lossa = 0.5 * NaNMath.sum(abs2.(1.0 .- (data[2, :] ./ sol_t[1, 1:end]))) / length(data[2, :])
            else
                lossa = 10.0^9 * length(data[2, :])

                #lossa = 0.5 * NaNMath.sum(abs2.(1.0 .- (data[2, 1:size(sol_t)[2]] ./ sol_t[1, 1:end]))) / length(data[2, 1:size(sol_t)[2]]) + 1000 * (size(data)[2] - size(sol_t)[2])
            end

            return lossa, sol
        end





        function loss_ode_blank_weighted_L2(p)
            sol = solve(ODE_prob, KenCarp4(autodiff=true), p=p, saveat=tsteps, verbose=false, abstol=1e-10, reltol=1e-10)
            sol_t = reduce(hcat, sol.u)
            empirical_blank_distrib = (blank_array .- mean(blank_array))
            # generation of the empirica distrib respect the mean of the noise 
            test = StatsBase.fit(Histogram, empirical_blank_distrib)
            binning_distrib = test.edges
            probabiltity_distrib = test.weights ./ sum(test.weights)
            if size(sol_t)[2] == size(data)[2]
                lossa = 0.0



                for ll in 1:size(sol_t)[2]
                    dist = (data[2, ll] - sol_t[1, ll])
                    index = findfirst(x -> (x > dist), binning_distrib[1])
                    if (typeof(index) == Nothing || index > length(probabiltity_distrib))

                        lossa = lossa + abs2.(dist) / length(data[2, :])
                    else

                        prob = probabiltity_distrib[index]


                        lossa = lossa + abs2.((1 - prob) * dist) / length(data[2, :])

                    end

                end
            else

                lossa = 10.0^9 * length(data[2, :])

            end

            return lossa, sol
        end



    end
    
    ## function for model predicifion

   
    # solving models that are not the huang model


    

    if model ∉ ["huang"  ,"HPM_exp" , "HPM", "dHPM","HPM_inhibition","dHPM_inhibition", "ODEs_HPM_SR","HPM_3_death","HPM_3_inhibition","HPM_3_death_resistance","dHPM_3_death_resistance"]

        ODE_prob =model_selector(model,u0,tspan)


# adtype = Optimization.AutoZygote()
# adtype = Optimization.AutoForwardDiff()
        adtype = ADTypes.AutoForwardDiff()

        if type_of_loss == "L2"
            optf = Optimization.OptimizationFunction((x, p) -> loss_ode_L2(x))
        end

        if type_of_loss == "RE"
            optf = Optimization.OptimizationFunction((x, p) -> loss_ode_RE(x))
        end

        if type_of_loss == "blank_weighted_L2"
            optf = Optimization.OptimizationFunction((x, p) -> loss_ode_blank_weighted_L2(x))
        end

        if type_of_loss == "L2_derivative"
            optf = Optimization.OptimizationFunction((x, p) -> loss_ode_L2_derivative(x))
        end

        optprob_const = Optimization.OptimizationProblem(optf, param, u0, lb=lb_param, ub=ub_param)
        res = Optimization.solve(optprob_const, BBO_adaptive_de_rand_1_bin_radiuslimited())


    end
    
    ## evaluation for HPM 
    if model == "HPM_3_death" || model == "HPM_3_inhibition"|| model == "HPM_3_death_resistance"|| model == "dHPM_3_death_resistance"
        if smoothing == true

            u0 = [data[2, 1],0.0,0.0]
    
        else
            u0 = [Statistics.mean(data[2, 1:pt_avg]),0.0,0.0]
        end
        # function selction here
        ODE_prob =model_selector(model,u0,tspan)
        
        function loss_ode_L2_3eq(p)
            sol = solve(ODE_prob, KenCarp4(autodiff=true), p=p, saveat=tsteps, verbose=false, abstol=1e-10, reltol=1e-10)
            sol_t = reduce(hcat, sol.u)
            sol_t_2 = sol_t[2,:]
            sol_t_2 =  reshape(sol_t_2, 1,length(sol_t_2))
            sol_t_3 = sol_t[3,:]
            sol_t_3 =  reshape(sol_t_3, 1,length(sol_t_3))
            sol_t_1 = sol_t[1,:]
            sol_t_1 =  reshape(sol_t_1, 1,length(sol_t_1))

            sol_t = sol_t_1 .+ sol_t_2 .+sol_t_3
 
            if size(sol_t)[2] == size(data)[2]
                lossa = NaNMath.sum(abs2.((data[2, :] .- sol_t[1, 1:end]))) / length(data[2, :])
            else
                lossa = 10.0^9 * length(data[2, :])

                            #                lossa = NaNMath.sum(abs2.((data[2, :] .- sol_t[1, 1:end]))) / length(data[2, :]) + 1000 * (size(data)[2] - size(sol_t)[2])
            end
        
            return lossa, sol
        end



        

        function loss_ode_RE_3eq(p)
            sol = solve(ODE_prob, KenCarp4(autodiff=true), p=p, saveat=tsteps, verbose=false, abstol=1e-10, reltol=1e-10)

            sol_t = reduce(hcat, sol.u)
            sol_t_2 = sol_t[2,:]
            sol_t_2 =  reshape(sol_t_2, 1,length(sol_t_2))

            sol_t_1 = sol_t[1,:]
            sol_t_1 =  reshape(sol_t_1, 1,length(sol_t_1))


            sol_t_3 = sol_t[3,:]
            sol_t_3 =  reshape(sol_t_3, 1,length(sol_t_3))

            sol_t = sol_t_1 .+ sol_t_2 .+sol_t_3
 

            if size(sol_t)[2] == size(data)[2]
                lossa = 0.5 * NaNMath.sum(abs2.(1.0 .- (data[2, :] ./ sol_t[1, 1:end]))) / length(data[2, :])
            else
                lossa = 10.0^9 * length(data[2, :])

#                lossa = 0.5 * NaNMath.sum(abs2.(1.0 .- (data[2, 1:size(sol_t)[2]] ./ sol_t[1, 1:end]))) / length(data[2, 1:size(sol_t)[2]]) + 1000 * (size(data)[2] - size(sol_t)[2]) 
            end

            return lossa, sol
        end        
        


        function loss_ode_blank_weighted_L2_3eq(p)
            sol = solve(ODE_prob, KenCarp4(autodiff=true), p=p, saveat=tsteps, verbose=false, abstol=1e-10, reltol=1e-10)
            sol_t = reduce(hcat, sol.u)
            sol_t_2 = sol_t[2,:]
            sol_t_2 =  reshape(sol_t_2, 1,length(sol_t_2))

            sol_t_1 = sol_t[1,:]
            sol_t_1 =  reshape(sol_t_1, 1,length(sol_t_1))

            sol_t_3 = sol_t[3,:]
            sol_t_3 =  reshape(sol_t_3, 1,length(sol_t_3))

            sol_t = sol_t_1 .+ sol_t_2 .+sol_t_3          
            empirical_blank_distrib = (blank_array .- mean(blank_array))
            # generation of the empirica distrib respect the mean of the noise 
            test = StatsBase.fit(Histogram, empirical_blank_distrib)
            binning_distrib = test.edges
            probabiltity_distrib = test.weights ./ sum(test.weights)
            if size(sol_t)[2] == size(data)[2]
                lossa = 0.0
        
        
        
                for ll in 1:size(sol_t)[2]
                    dist = (data[2, ll] - sol_t[1, ll])
                    index = findfirst(x -> (x > dist), binning_distrib[1])
                    if (typeof(index) == Nothing || index > length(probabiltity_distrib))
       
                        lossa = lossa + abs2.(dist) / length(data[2, :])
                    else
        
                        prob = probabiltity_distrib[index]
        
        
                        lossa = lossa + abs2.((1 - prob) * dist) / length(data[2, :])
        
                    end
        
                end
                
            
            else
        
                lossa = 10.0^9 * length(data[2, :])
        
            end
        
            return lossa, sol
        end      


        function loss_ode_L2_derivative_3eq(p)
          sol = solve(ODE_prob, KenCarp4(autodiff=true), p=p, saveat=tsteps, verbose=false, abstol=1e-10, reltol=1e-10)

          sol_t = reduce(hcat, sol.u)
          sol_t_2 = sol_t[2,:]
          sol_t_2 =  reshape(sol_t_2, 1,length(sol_t_2))

          sol_t_1 = sol_t[1,:]
          sol_t_1 =  reshape(sol_t_1, 1,length(sol_t_1))

          sol_t_3 = sol_t[3,:]
          sol_t_3 =  reshape(sol_t_3, 1,length(sol_t_3))

          sol_t = sol_t_1 .+ sol_t_2 .+sol_t_3        
            if size(sol_t)[2] == size(data)[2]
                lossa = 0.5 * NaNMath.sum(abs2.(1.0 .- (data[2, :] ./ sol_t[1, 1:end]))) / length(data[2, :])
            else
                lossa = 10.0^9 * length(data[2, :])

#                lossa = 0.5 * NaNMath.sum(abs2.(1.0 .- (data[2, 1:size(sol_t)[2]] ./ sol_t[1, 1:end]))) / length(data[2, 1:size(sol_t)[2]]) + 1000 * (size(data)[2] - size(sol_t)[2])
            end
        
            return lossa, sol
        end      
        

        
# adtype = Optimization.AutoZygote()
# adtype = Optimization.AutoForwardDiff()
        adtype = ADTypes.AutoForwardDiff()

        if type_of_loss == "L2"
            optf = Optimization.OptimizationFunction((x, p) -> loss_ode_L2_3eq(x))
        end

        if type_of_loss == "RE"
            optf = Optimization.OptimizationFunction((x, p) -> loss_ode_RE_3eq(x))
        end

        if type_of_loss == "blank_weighted_L2"
            optf = Optimization.OptimizationFunction((x, p) -> loss_ode_blank_weighted_L2_3eq(x))
        end

        if type_of_loss == "L2_derivative"
            optf = Optimization.OptimizationFunction((x, p) -> loss_ode_L2_derivative_3eq(x))
        end
        optprob_const = Optimization.OptimizationProblem(optf, param, u0, lb=lb_param, ub=ub_param)
        res = Optimization.solve(optprob_const, BBO_adaptive_de_rand_1_bin_radiuslimited())
        


    
    
    end

    if model == "HPM" ||  model == "dHPM" ||  model == "HPM_inhibition" ||  model == "dHPM_inhibition" || model == "ODEs_HPM_SR" || model == "HPM_exp"
        # specific initial condition for system of ODEs
        # all the biomass starts as dormient


        
        if smoothing == true

            u0 = [data[2, 1],0.0]
    
        else
            u0 = [Statistics.mean(data[2, 1:pt_avg]),0.0]
        end
         #function selection
         ODE_prob =model_selector(model,u0,tspan)


        function loss_ode_L2_2eq(p)
            sol = solve(ODE_prob, KenCarp4(autodiff=true), p=p, saveat=tsteps, verbose=false, abstol=1e-10, reltol=1e-10)
            sol_t = reduce(hcat, sol.u)
            sol_t_2 = sol_t[2,:]
            sol_t_2 =  reshape(sol_t_2, 1,length(sol_t_2))

            sol_t_1 = sol_t[1,:]
            sol_t_1 =  reshape(sol_t_1, 1,length(sol_t_1))

            sol_t = sol_t_1 .+ sol_t_2
 
            if size(sol_t)[2] == size(data)[2]
                lossa = NaNMath.sum(abs2.((data[2, :] .- sol_t[1, 1:end]))) / length(data[2, :])
            else
                lossa = 10.0^9 * length(data[2, :])

                            #                lossa = NaNMath.sum(abs2.((data[2, :] .- sol_t[1, 1:end]))) / length(data[2, :]) + 1000 * (size(data)[2] - size(sol_t)[2])
            end
        
            return lossa, sol
        end



        

        function loss_ode_RE_2eq(p)
            sol = solve(ODE_prob, KenCarp4(autodiff=true), p=p, saveat=tsteps, verbose=false, abstol=1e-10, reltol=1e-10)

            sol_t = reduce(hcat, sol.u)
            sol_t_2 = sol_t[2,:]
            sol_t_2 =  reshape(sol_t_2, 1,length(sol_t_2))

            sol_t_1 = sol_t[1,:]
            sol_t_1 =  reshape(sol_t_1, 1,length(sol_t_1))

            sol_t = sol_t_1 .+ sol_t_2
 

            if size(sol_t)[2] == size(data)[2]
                lossa = 0.5 * NaNMath.sum(abs2.(1.0 .- (data[2, :] ./ sol_t[1, 1:end]))) / length(data[2, :])
            else
                lossa = 10.0^9 * length(data[2, :])

#                lossa = 0.5 * NaNMath.sum(abs2.(1.0 .- (data[2, 1:size(sol_t)[2]] ./ sol_t[1, 1:end]))) / length(data[2, 1:size(sol_t)[2]]) + 1000 * (size(data)[2] - size(sol_t)[2]) 
            end

            return lossa, sol
        end        
        


        function loss_ode_blank_weighted_L2_2eq(p)
            sol = solve(ODE_prob, KenCarp4(autodiff=true), p=p, saveat=tsteps, verbose=false, abstol=1e-10, reltol=1e-10)
            sol_t = reduce(hcat, sol.u)
            sol_t_2 = sol_t[2,:]
            sol_t_2 =  reshape(sol_t_2, 1,length(sol_t_2))

            sol_t_1 = sol_t[1,:]
            sol_t_1 =  reshape(sol_t_1, 1,length(sol_t_1))

            sol_t = sol_t_1 .+ sol_t_2
          
            empirical_blank_distrib = (blank_array .- mean(blank_array))
            # generation of the empirica distrib respect the mean of the noise 
            test = StatsBase.fit(Histogram, empirical_blank_distrib)
            binning_distrib = test.edges
            probabiltity_distrib = test.weights ./ sum(test.weights)
            if size(sol_t)[2] == size(data)[2]
                lossa = 0.0
        
        
        
                for ll in 1:size(sol_t)[2]
                    dist = (data[2, ll] - sol_t[1, ll])
                    index = findfirst(x -> (x > dist), binning_distrib[1])
                    if (typeof(index) == Nothing || index > length(probabiltity_distrib))
       
                        lossa = lossa + abs2.(dist) / length(data[2, :])
                    else
        
                        prob = probabiltity_distrib[index]
        
        
                        lossa = lossa + abs2.((1 - prob) * dist) / length(data[2, :])
        
                    end
        
                end
                
            
            else
        
                lossa = 10.0^9 * length(data[2, :])
        
            end
        
            return lossa, sol
        end      


        function loss_ode_L2_derivative_2eq(p)
          sol = solve(ODE_prob, KenCarp4(autodiff=true), p=p, saveat=tsteps, verbose=false, abstol=1e-10, reltol=1e-10)

          sol_t = reduce(hcat, sol.u)
          sol_t_2 = sol_t[2,:]
          sol_t_2 =  reshape(sol_t_2, 1,length(sol_t_2))

          sol_t_1 = sol_t[1,:]
          sol_t_1 =  reshape(sol_t_1, 1,length(sol_t_1))

          sol_t = sol_t_1 .+ sol_t_2         
        
            if size(sol_t)[2] == size(data)[2]
                lossa = 0.5 * NaNMath.sum(abs2.(1.0 .- (data[2, :] ./ sol_t[1, 1:end]))) / length(data[2, :])
            else
                lossa = 10.0^9 * length(data[2, :])

#                lossa = 0.5 * NaNMath.sum(abs2.(1.0 .- (data[2, 1:size(sol_t)[2]] ./ sol_t[1, 1:end]))) / length(data[2, 1:size(sol_t)[2]]) + 1000 * (size(data)[2] - size(sol_t)[2])
            end
        
            return lossa, sol
        end      
        

        
# adtype = Optimization.AutoZygote()
# adtype = Optimization.AutoForwardDiff()
        adtype = ADTypes.AutoForwardDiff()

        if type_of_loss == "L2"
            optf = Optimization.OptimizationFunction((x, p) -> loss_ode_L2_2eq(x))
        end

        if type_of_loss == "RE"
            optf = Optimization.OptimizationFunction((x, p) -> loss_ode_RE_2eq(x))
        end

        if type_of_loss == "blank_weighted_L2"
            optf = Optimization.OptimizationFunction((x, p) -> loss_ode_blank_weighted_L2_2eq(x))
        end

        if type_of_loss == "L2_derivative"
            optf = Optimization.OptimizationFunction((x, p) -> loss_ode_L2_derivative_2eq(x))
        end
        optprob_const = Optimization.OptimizationProblem(optf, param, u0, lb=lb_param, ub=ub_param)
        res = Optimization.solve(optprob_const, BBO_adaptive_de_rand_1_bin_radiuslimited())
        


    end




 
    #revalution of solution for plot an loss evaluation 


        remade_solution = solve(remake(ODE_prob, p=res.u), KenCarp4(), saveat=tsteps)
        sol_fin = reduce(hcat, remade_solution.u)


    if model == "HPM_3_death" || model == "HPM_3_inhibition"|| model == "HPM_3_death_resistance"|| model == "dHPM_3_death_resistance"

        sol_t_2 = sol_fin[2,:]
        sol_t_2 =  reshape(sol_t_2, 1,length(sol_t_2))

        sol_t_1 = sol_fin[1,:]
        sol_t_1 =  reshape(sol_t_1, 1,length(sol_t_1))


        sol_t_3 = sol_fin[3,:]
        sol_t_3 =  reshape(sol_t_3, 1,length(sol_t_3))

        sol_fin = sol_t_3 .+ sol_t_2.+ sol_t_1



    end 


    if (model == "HPM" || model == "HPM_exp" )

        sol_t_2 = sol_fin[2,:]
        sol_t_2 =  reshape(sol_t_2, 1,length(sol_t_2))

        sol_t_1 = sol_fin[1,:]
        sol_t_1 =  reshape(sol_t_1, 1,length(sol_t_1))

        sol_fin = sol_t_1 .+ sol_t_2



    end 
    if model == "ODEs_HPM_SR"

        sol_t_2 = sol_fin[2,:]
        sol_t_2 =  reshape(sol_t_2, 1,length(sol_t_2))

        sol_t_1 = sol_fin[1,:]
        sol_t_1 =  reshape(sol_t_1, 1,length(sol_t_1))

        sol_fin = sol_t_1 .+ sol_t_2
    



    end 

    if model == "dHPM"

        sol_t_2 = sol_fin[2,:]
        sol_t_2 =  reshape(sol_t_2, 1,length(sol_t_2))

        sol_t_1 = sol_fin[1,:]
        sol_t_1 =  reshape(sol_t_1, 1,length(sol_t_1))

        sol_fin = sol_t_1 .+ sol_t_2

    end 


    if model == "HPM_inhibition"

        sol_t_2 = sol_fin[2,:]
        sol_t_2 =  reshape(sol_t_2, 1,length(sol_t_2))

        sol_t_1 = sol_fin[1,:]
        sol_t_1 =  reshape(sol_t_1, 1,length(sol_t_1))

        sol_fin = sol_t_1 .+ sol_t_2

    end 

    if model == "dHPM_inhibition"

        sol_t_2 = sol_fin[2,:]
        sol_t_2 =  reshape(sol_t_2, 1,length(sol_t_2))

        sol_t_1 = sol_fin[1,:]
        sol_t_1 =  reshape(sol_t_1, 1,length(sol_t_1))

        sol_fin = sol_t_1 .+ sol_t_2

    end 


  



    # plotting if required
    if do_plot == true
        mkpath(path_to_plot)
        display(Plots.scatter(data[1, :], data[2, :], xlabel="Time", ylabel="Arb. Units", label=["Data " nothing], markersize=2, color=:black, title=string(label_exp, " ", name_well)))
        display(Plots.plot!(remade_solution.t, sol_fin[1, 1:end], xlabel="Time", ylabel="Arb. Units", label=[string("Fitting ", model) nothing], c=:red))
        png(string(path_to_plot, label_exp, "_", model, "_", name_well, ".png"))
    end


    if type_of_loss == "L2"

        if size(sol_fin)[2] == size(data)[2]
            loss = NaNMath.sum(abs2.((data[2, :] - sol_fin[1, 1:end]))) / length(data[2, :])
        else
            loss = NaNMath.sum(abs2.((data[2, 1:size(sol_fin)[2]] - sol_fin[1, 1:end]))) / length(data[2, 1:size(sol_fin)[2]]) + 1000 * (size(data)[2] - size(sol_fin)[2])
        end
    end


    if type_of_loss == "RE"

        if size(sol_fin)[2] == size(data)[2]
            loss = 0.5 * NaNMath.sum(abs2.(1.0 .- (data[2, :] ./ sol_fin[1, 1:end]))) / length(data[2, :])
        else
            loss = 0.5 * NaNMath.sum(abs2.(1.0 .- (data[2, 1:size(sol_fin)[2]] ./ sol_fin[1, 1:end]))) / length(data[2, 1:size(sol_fin)[2]]) + 1000 * (size(data)[2] - size(sol_fin)[2])
        end
    end

    if type_of_loss == "L2_derivative"

        derivative_theo = [(sol_fin[1, i] - sol_fin[1, i-1]) / (data[1, i] - data[1, i-1]) for i in 2:(size(sol_fin)[2])]

        derivative_data = [(data[2, i] - data[2, i-1]) / (data[1, i] - data[1, i-1]) for i in 2:(size(data)[2])]

        if size(derivative_theo)[1] == size(derivative_data)[1]
            loss = NaNMath.sum(abs2.((derivative_theo[2, :] .- derivative_data[1, 1:end]))) / length(derivative_data[2, :])
        else
            loss = NaNMath.sum(abs2.((derivative_theo[2, :] .- derivative_data[1, 1:end]))) / length(derivative_data[2, :]) + 1000
        end

    end

    if type_of_loss == "blank_weighted_L2"


        empirical_blank_distrib = (blank_array .- mean(blank_array))
        # generation of the empirica distrib respect the mean of the noise 
        test = StatsBase.fit(Histogram, empirical_blank_distrib)
        binning_distrib = test.edges
        probabiltity_distrib = test.weights ./ sum(test.weights)

        if size(sol_fin)[2] == size(data)[2]
            loss = 0.0



            for ll in 1:size(sol_fin)[2]
                dist = (data[2, ll] - sol_fin[1, ll])
                index = findfirst(x -> (x > dist), binning_distrib[1])
                if (typeof(index) == Nothing || index > length(probabiltity_distrib))

                    loss = loss + (abs2.(dist)) / length(data[2, :])
                else
                    prob = probabiltity_distrib[index]

                    loss = loss + (abs2.((1 - prob) * dist)) / length(data[2, :])

                end

            end
        else

            loss = 10.0^9 * length(data[2, :])

        end


    end


    # parameters have different form in function of the used model
   
    #max_theoretical gr
    max_th_gr = maximum( specific_gr_evaluation(  Matrix(transpose(hcat(data[1,:],  sol_fin[1, 1:end]))), pt_smooth_derivative))
    # max empirical gr
    max_em_gr = maximum( specific_gr_evaluation(data, pt_smooth_derivative))

    res_temp = res.u

    res_param = vectorize_df_results(name_well,
        model,
        res_temp,
        max_th_gr,
        max_em_gr,
        loss
    )

    gr_param = res_temp[1]


    if error_analysis == true
        error_metrics(
            label_exp, #label of the experiment
            data,# name of the file to analyze
            path_to_plot, # path where save plots
            sol_fin, # 
            name_well,
            gr_param,
            model,
            true, # 1 do and visulaze the plots of data
            write_res,
            correction
        )
    end



    return res_param

end
##########################################################################
##########################################################################



function fit_one_file_ODE(
    label_exp::String, #label of the experiment
    path_to_data::String, # path to the folder to analyze
    path_to_annotation::String,# path to the annotation of the wells
    path_to_results::String, # path where save results
    path_to_plot::String, # path where to save Plots
    model::String, # string of the used model
    loss_type::String, # string of the type of the used loss
    smoothing::Bool, # 1 do smoothing of data with rolling average
    do_plot::Bool, # 1 do and visulaze the plots of data
    verbose::Bool, # 1 true verbose
    write_res::Bool, # write results
    pt_avg::Int, # number of points to do smoothing average
    pt_smooth_derivative::Int, # number of points to do ssmooth_derivative
    lb_param::Vector{Float64},# array of the array of the lower bound of the parameters
    ub_param::Vector{Float64}, # array of the array of the upper bound of the parameters
    blank_subtraction::String, # string on how to use blank (NO,avg_subctraction,time_avg)
    fit_replicate::Bool, # if true the average between replicates is fitted. If false all replicate are fitted indipendelitly
    error_analysis::Bool, # do error analysis
    correct_negative::String, # if "thr_correction" it put a thr on the minimum value of the data with blank subracted, if "blank_correction" uses blank distrib to impute negative values
    thr_negative::Float64,  # used only if correct_negative == "thr_correction"
    correction::Bool, # if true uses the given calibration curve to fix the data
    calibration_OD_curve ::String  #  the path to calibration curve to fix the data
    )




    mkpath(path_to_results)

    if do_plot == true
        mkpath(path_to_plot)
    end

    parameter_of_optimization = inzialize_df_results(model)

 

    ## reading annotation here
    annotation = CSV.File(string(path_to_annotation))
    names_of_annotated_df = [annotation[l][1] for l in 1:length(annotation)]
    # selcting blank wells
    properties_of_annotation = [annotation[l][2] for l in 1:length(annotation)]
    list_of_blank = names_of_annotated_df[findall(x -> x == "b", properties_of_annotation)]
    list_of_discarded = names_of_annotated_df[findall(x -> x == "X", properties_of_annotation)]
    list_of_blank = Symbol.(list_of_blank)
    list_of_discarded = Symbol.(list_of_discarded)



    # reading files
    dfs_data = CSV.File(path_to_data)


    # shaping df for the inference

    names_of_cols = propertynames(dfs_data)


    # excluding blank data and discarded wells 
    names_of_cols = filter!(e -> !(e in list_of_blank), names_of_cols)
    if length(list_of_discarded)>0
        names_of_cols = filter!(e -> !(e in list_of_discarded), names_of_cols)
    end


    times_data = dfs_data[names_of_cols[1]]

    blank_array = reduce(vcat, [(dfs_data[k]) for k in list_of_blank])
    blank_array = convert(Vector{Float64}, blank_array)

    ## BLANK ANALYSIS HERE 
    if blank_subtraction == "avg_blank"

        blank_value = mean([mean(dfs_data[k]) for k in list_of_blank])

    elseif blank_subtraction == "time_blank"

        blank_value = [mean([dfs_data[k][j] for k in list_of_blank]) for j in 1:length(times_data)]

    else

        blank_value = zeros(length(times_data))

    end




    ## considering replicates
    list_replicate = unique(properties_of_annotation)
    list_replicate = filter!(e -> e != "b", list_replicate)

    if fit_replicate == true
        new_data = times_data

        list_replicate = unique(properties_of_annotation)
        list_replicate = filter!(e -> e != "b", list_replicate)


        for replicate_temp in list_replicate

            names_of_replicate_temp = Symbol.(names_of_annotated_df[findall(x -> x == replicate_temp, properties_of_annotation)])
            replicate_mean = [mean([dfs_data[k][j] for k in names_of_replicate_temp]) for j in 1:length(times_data)]

            new_data = hcat(new_data, replicate_mean)

        end
        new_data = DataFrame(new_data, :auto)
        rename!(new_data, vcat(:Time, reduce(vcat, Symbol.(list_replicate))))
        names_of_cols = propertynames(new_data)
        dfs_data = new_data


    end



    # for on the columns to analyze

    for well_name in names_of_cols[2:end]




        if fit_replicate == true

            data_values = copy(dfs_data[!, well_name])

        else
            data_values = copy(dfs_data[well_name])
        end

        # blank subtraction 
        data_values = data_values .- blank_value

        data = Matrix(transpose(hcat(times_data, data_values)))

        if correct_negative == "thr_correction"

            data = thr_negative_correction(data, thr_negative)

        end

        if correct_negative == "blank_correction"

            data = blank_distrib_negative_correction(data, blank_array)

        end


        if smoothing == true

            data = smoothing_data(data, pt_avg)
        end
        # defining time steps of the inference

        max_t = data[1, end]
        min_t = data[1, 1]
        # attention, this is necessart if data are not float
        #println(typeof(data[2,:]) )
        #data[2,:] = parse(Float64,data[2,:] )
        # guessing  parameters

        data = Matrix(data)

        param = guess_param(model, data, min_t, max_t,lb_param,ub_param)


        # inference

        temp_results_1 = fitting_one_well_constrained(data, # dataset x times y OD/fluorescence
            string(well_name), # name of the well
            label_exp, #label of the experiment
            lb_param, # lower bound param
            ub_param, # upper bound param
            param, # initial guess param
            model, # ode model to use 
            do_plot, # do plots or no
            path_to_plot, # where save plots
            pt_avg, # number of the point to generate intial condition
            pt_smooth_derivative,
            smoothing, # the smoothing is done or not?
            loss_type, # type of used loss 
            blank_array, # data of all blanks
            error_analysis,# do error analysis
            write_res,
            correction, # if true uses the given calibration curve to fix the data
            calibration_OD_curve   #  the path to calibration curve to fix the data
            )


        if verbose == true
            println("the results are:")
            println(temp_results_1)
        end

        parameter_of_optimization = hcat(parameter_of_optimization, temp_results_1)

    end


    if write_res == true

        CSV.write(string(path_to_results, label_exp, "_parameters_", model, ".csv"), Tables.table(Matrix(parameter_of_optimization)))


    end
    return parameter_of_optimization




end








function fit_one_file_Log_Lin(
    label_exp::String, #label of the experiment
    path_to_data::String, # path to the folder to analyze
    path_to_annotation::String,# path to the annotation of the wells
    path_to_results::String, # path where save results
    path_to_plot::String, # path where to save Plots
    do_plot::Bool, # 1 do and visulaze the plots of data
    verbose::Bool, # 1 true verbose
    write_res::Bool, # write results
    type_of_smoothing::String, # option, NO, gaussian, rolling avg
    pt_avg::Int, # number of points to do smoothing average
    pt_smoothing_derivative::Int, # number of poits to smooth the derivative
    pt_min_size_of_win::Int, # minimum size of the exp windows in number of smooted points
    type_of_win::String, # how the exp. phase win is selected, "maximum" of "global_thr"
    threshold_of_exp::Float64, # threshold of growth rate in quantile to define the exp windows
    blank_subtraction::String, # string on how to use blank (NO,avg_subctraction,time_avg)
    fit_replicate::Bool, # if true the average between replicates is fitted. If false all replicate are fitted indipendelitly
    correct_negative::String, # if "thr_correction" it put a thr on the minimum value of the data with blank subracted, if "blank_correction" uses blank distrib to impute negative values
    thr_negative::Float64, # used only if correct_negative == "thr_correction"
    correction::Bool, # if true uses the given calibration curve to fix the data
    calibration_OD_curve ::String  #  the path to calibration curve to fix the data
    )



    mkpath(path_to_results)
    if do_plot == true
        mkpath(path_to_plot)
    end





    # reading files
  # dfs_data = CSV.File(path_to_data,header=true,sep=",")
   dfs_data = CSV.File(path_to_data)

    # TEMPORARY results df
    results_Log_Lin = ["label_exp", "well_name", "t_start", "t_end", "t_of_max",  "empirical_max_Growth_rate" ,"Growth_rate", "2sigma_gr","dt","95_confidence_dt_upper","95_confidence_dt_lower","intercept","2sigma_intercept","R^2"]


    # shaping df for the inference

    names_of_cols = propertynames(dfs_data)
    times_data = dfs_data[names_of_cols[1]]
    annotation = CSV.File(string(path_to_annotation),header=false)
    names_of_annotated_df = [annotation[l][1] for l in 1:length(annotation)]
    # selcting blank wells
    properties_of_annotation = [annotation[l][2] for l in 1:length(annotation)]
    list_of_blank = names_of_annotated_df[findall(x -> x == "b", properties_of_annotation)]
    list_of_discarded = names_of_annotated_df[findall(x -> x == "X", properties_of_annotation)]
    list_of_blank = Symbol.(list_of_blank)
    list_of_discarded = Symbol.(list_of_discarded)






    # excluding blank data and discarded wells 
    names_of_cols = filter!(e -> !(e in list_of_blank), names_of_cols)
    if length(list_of_discarded)>0
        names_of_cols = filter!(e -> !(e in list_of_discarded), names_of_cols)
    end


    times_data = dfs_data[names_of_cols[1]]
    ## BLANK ANALYSIS HERE 


    blank_array = reduce(vcat, [(dfs_data[k]) for k in list_of_blank])
    blank_array = convert(Vector{Float64}, blank_array)


    if blank_subtraction == "avg_blank"

        blank_value = mean([mean(dfs_data[k]) for k in list_of_blank])

    elseif blank_subtraction == "time_blank"

        blank_value = [mean([dfs_data[k][j] for k in list_of_blank]) for j in 1:length(times_data)]

    else

        blank_value = zeros(length(times_data))

    end




    ## considering replicates
    list_replicate = unique(properties_of_annotation)
    list_replicate = filter!(e -> e != "b", list_replicate)

    if fit_replicate == true
        new_data = times_data

        list_replicate = unique(properties_of_annotation)
        list_replicate = filter!(e -> e != "b", list_replicate)


        for replicate_temp in list_replicate

            names_of_replicate_temp = Symbol.(names_of_annotated_df[findall(x -> x == replicate_temp, properties_of_annotation)])
            replicate_mean = [mean([dfs_data[k][j] for k in names_of_replicate_temp]) for j in 1:length(times_data)]

            new_data = hcat(new_data, replicate_mean)

        end
        new_data = DataFrame(new_data, :auto)
        rename!(new_data, vcat(:Time, reduce(vcat, Symbol.(list_replicate))))
        names_of_cols = propertynames(new_data)
        dfs_data = new_data


    end

    # for on the columns to analyze

    for well_name in names_of_cols[2:end]



        if fit_replicate == true

            data_values = copy(dfs_data[!, well_name])

        else
            data_values = copy(dfs_data[well_name])
        end

        # blank subtraction 
        data_values = data_values .- blank_value


        data = Matrix(transpose(hcat(times_data, data_values)))


        if correct_negative == "thr_correction"

            data = thr_negative_correction(data, thr_negative)

        end

        if correct_negative == "blank_correction"

            data = blank_distrib_negative_correction(data, blank_array)

        end

        data = Matrix(transpose(hcat(data[1, :], data[2, :])))

        # inference

        temp_results_1 = fitting_one_well_Log_Lin_Fit(data, # dataset x times y OD/fluorescence
            string(well_name), # name of the well
            label_exp, #label of the experiment
            do_plot, # do plots or no
            path_to_plot, # where save plots
            type_of_smoothing,# option, NO, gaussian, rolling avg
            pt_avg, # numebr of the point to generate intial condition
            pt_smoothing_derivative, # number of poits to smooth the derivative
            pt_min_size_of_win, # minimum size of the exp windows in number of smooted points
           type_of_win, # how the exp. phase win is selected, "maximum" of "global_thr"
            threshold_of_exp, # threshold of growth rate in quantile to define the exp windows
            correction, # if true uses the given calibration curve to fix the data
            calibration_OD_curve  #  the path to calibration curve to fix the data
            )

        if verbose == true
            println("the results are:")
            println(temp_results_1)
        end

        results_Log_Lin = hcat(results_Log_Lin, temp_results_1)




        if write_res == true

            CSV.write(string(path_to_results, label_exp, "_results.csv"), Tables.table(Matrix(results_Log_Lin)))


        end

    end

    return results_Log_Lin




end










######################## function for the  plots and  metric of errors ############################################## 







function error_metrics(
    label_exp::String, #label of the experiment
    data::Matrix{Float64},# name of the file to analyze
    path_to_plot::String, # path where save plots
    sol_fin::Any, # path where to save Plots
    name_well::String,
    gr_param::Float64,
    model::String,
    do_plot::Bool, # 1 do and visulaze the plots of data
    write_res::Bool,
     correction::Bool)

    ## absolute error dynamics

    abs_error_dynamics = abs.(data[2, :] .- sol_fin[1, 1:end])


    if write_res == true

        CSV.write(string(path_to_results, label_exp, "_abs_error_dynamics_", model, ".csv"), Tables.table((abs_error_dynamics)))


    end


    if do_plot == true

        mkpath(path_to_plot)
        display(Plots.plot(data[1, :], abs_error_dynamics, xlabel="Time", ylabel="abs error Arb. Units", label=[string("Fitting ", model) nothing], c=:black))
        png(string(path_to_plot, label_exp, "_abs_error_dynamics_", model, "_", name_well, ".png"))
        display(Plots.plot(data[1, :], cumsum(abs_error_dynamics), xlabel="Time", ylabel="cumulative abs error Arb. Units", label=[string("Fitting ", model) nothing], c=:black))
        png(string(path_to_plot, label_exp, "_cumulative_error_dynamics_", model, "_", name_well, ".png"))


    end



    ##  residual distrib plot
    residuals_error = data[2, :] .- sol_fin[1, 1:end]


    if write_res == true


        CSV.write(string(path_to_results, label_exp, "_distrib_residual_", model, ".csv"), Tables.table((residuals_error)))


    end


    # plots of residuals
    if do_plot == true

        mkpath(path_to_plot)
        b_range = range(minimum(residuals_error), maximum(residuals_error), length=100)
        display(Plots.histogram(residuals_error, xlabel="errors ", ylabel="counts", bins=b_range, normalize=:pdf, label=[string("Histogram residuals ") nothing]))
        png(string(path_to_plot, label_exp, "_residual_distrib_", model, "_", name_well, ".png"))

    end
    # comparison of dynamics of growth rate and max growth
    #pt_avg = 7
   # pt_smoothing_derivative = 7
    pt_min_size_of_win = 4
    result_data = fitting_one_well_Log_Lin_Fit_for_errors(data, # dataset x times y OD/fluorescence
        pt_avg, # numebr of the point to generate intial condition
        pt_smoothing_derivative, # number of poits to smooth the derivative
        pt_min_size_of_win, # minimum size of the exp windows in number of smooted points
        0.8, # threshold of growth rate in quantile to define the exp windows
    correction)

    synth_data = Matrix(transpose(hcat(reduce(vcat, data[1, :]), reduce(vcat, sol_fin[1, 1:end]))))


    result_synth = fitting_one_well_Log_Lin_Fit_for_errors(synth_data, # dataset x times y OD/fluorescence
        pt_avg, # numebr of the point to generate intial condition
        pt_smoothing_derivative, # number of poits to smooth the derivative
        pt_min_size_of_win, # minimum size of the exp windows in number of smooted points
        0.98,    # threshold of growth rate in quantile to define the exp windows
        correction)


    # plotting if requested
    if do_plot == true
        mkpath(path_to_plot)
        display(Plots.scatter(result_data[5], result_data[4], xlabel="Time", ylabel="GR 1/T", label=["Data " nothing], markersize=1, color=:black, title=string(label_exp, " ", name_well)))
        display(Plots.plot!(result_synth[5], result_synth[4], xlabel="Time ", ylabel="GR 1/T", label=[string("Fitting ") nothing], c=:red))
        display(Plots.vline!([result_data[2], result_data[3]], c=:black, label=[string("Exp. phase ") nothing]))
     #   display(Plots.hline!([result_data[1]], c=:black, label=nothing))
        display(Plots.hline!([gr_param[1]], c=:red, label=nothing))




        png(string(path_to_plot, label_exp, "gr_comparison", name_well, ".png"))

    end







end
########################function for the simulations############################################## 

##### poisson approx simulation ###############
function stochastic_sim(model::String, #string of the model
    n_start::Int, # number of starting cells
    n_mol_start::Float64, # starting concentration of the limiting nutrients
    tstart::Float64, # start time of the sim
    tmax::Float64, # final time of the sim
    delta_t::Float64, # delta t for poisson approx
    k_1_val::Float64,
    k_2_val::Float64, # monod constant
    alpha_val::Float64, # massimum possible growth rate
    lambda::Float64, # lag time
    n_mol_per_birth::Float64,# nutrient consumed per division (conc)
    volume::Float64
)


    #inizialization of times 
    tot_pop = [copy(n_start)]
    times = [copy(tstart)]

    conc_of_nutriens = [copy(n_mol_start / volume)]
    n_times = floor((tmax - tstart) / delta_t)




    for i in 2:n_times

        # if for defining the lag phase

        if i * delta_t < lambda
            rate_per_cell = 0.0
        else

            if model == "Monod"
                factor = conc_of_nutriens[end] / (k_1_val + conc_of_nutriens[end])
            end

            if model == "Haldane"
                factor = conc_of_nutriens[end] / (k_1_val + conc_of_nutriens[end] + (conc_of_nutriens[end])^2 / k_2_val)

            end

            if model == "Blackman"
                factor = conc_of_nutriens[end] / (k_1_val)
            end



            if model == "Tesseir"
                factor = 1 - exp(conc_of_nutriens[end] * k_1_val)
            end

            if model == "Moser"
                factor = conc_of_nutriens[end]^k_2_val / (conc_of_nutriens[end]^k_2_val + k_1_val)
            end



            if model == "Aiba-Edwards"
                factor = exp(-k_2_val / conc_of_nutriens[end]) * conc_of_nutriens[end] / (conc_of_nutriens[end] + k_1_val)
            end

            if model == "Verhulst"
                factor = 1 - tot_pop[end] / k_1_val
            end

            rate_per_cell = alpha_val * factor
        end


        # evaluating the number of birth events with poisson approx

        total_birth_rate = rate_per_cell .* tot_pop[end] .* delta_t


        n_birth = rand(Poisson(total_birth_rate), 1)




        # updating states 
        new_conc = max(0, conc_of_nutriens[end] - n_birth[1] * n_mol_per_birth / volume)
        net_pop_variation = n_birth[1]

        new_pop = max(0, tot_pop[end] + net_pop_variation)
        conc_of_nutriens = push!(conc_of_nutriens, new_conc)
        tot_pop = push!(tot_pop, new_pop)
        times = push!(times, times[end] + delta_t)




    end

    return tot_pop, conc_of_nutriens, times


end





########## ode simulation ################



function ODE_sim(model::String, #string of the model
    n_start::Vector{Float64}, # starting condition
    tstart::Float64, # start time of the sim
    tmax::Float64, # final time of the sim
    delta_t::Float64, # delta t for poisson approx
    integrator::Any, # which sciml solver of ode
    param_of_ode::Vector{Float64} # parameters of the ODE model
)

    # defining time stepping
    t_steps = tstart:delta_t:tmax
    tspan = (tstart, tmax)

    u0 = n_start

    ODE_prob =model_selector_sim(model,u0,tspan,param_of_ode)



    sim = solve(ODE_prob, integrator, saveat=t_steps)



    return sim
end




function ODE_sim_for_iterate(model::String, #string of the model
    n_start::Vector{Float64}, # starting condition
    array_time::Vector{Float64},
    integrator::Any, # which sciml solver of ode
    param_of_ode::Any # parameters of the ODE model
)

    # defining time stepping
    t_steps = array_time
    tspan = (array_time[1],array_time[end])

    u0 = n_start

    ODE_prob =model_selector_sim(model,u0,tspan,param_of_ode)



    sim = solve(ODE_prob, integrator, saveat=t_steps)

    return sim
end












function generation_of_combination_of_IC_morris(lb_param::Vector{Float64},
    ub_param::Vector{Float64},
    N_step_morris::Int
)

    starting_guess = copy(lb_param)


    delta_vec = [(ub_param[i] - lb_param[i]) / (N_step_morris + 1) for i in 1:length(ub_param)]

    # generating combinations of all possible parameters values

    combinations_par = copy(starting_guess)
    combinations_tot = copy(combinations_par)

    for k in 1:N_step_morris
        for ll in 1:(size(delta_vec)[1])


            combinations_par[ll] = combinations_par[ll] + delta_vec[ll]
            combinations_tot = hcat(combinations_tot, combinations_par)

        end
    end




    println(string("the  number combinations of inizializaton to test is ", size(combinations_tot)[2]))

    return combinations_tot


end


function Sensitivity_ODE_Morris(data::Matrix{Float64}, # dataset x times y OD/fluorescence
    name_well::String, # name of the well
    label_exp::String, #label of the experiment
    lb_param::Vector{Float64}, # lower bound param
    ub_param::Vector{Float64}, # upper bound param
    model::String, # ode model to use 
    do_plot::Bool, # do plots or no
    path_to_plot::String, # where save plots
    path_to_results::String, # where save data frame
    pt_avg::Int, # numebr of the point to generate intial condition
    smoothing::Bool, # the smoothing is done or not?
    type_of_loss::String, # type of used loss 
    blank_array::Vector{Float64}, # data of all blanks
    N_step_morris::Int,
    write_res::Bool,
    correction::Bool, # if true uses the given calibration curve to fix the data
    calibration_OD_curve ::String  #  the path to calibration curve to fix the data
    
)

    if  correction == true

     data = correction_OD_multiple_scattering(data,calibration_OD_curve)


    end

    results_sensitivity =inzialize_df_results(model)

    if write_res == true
        mkpath(path_to_results)
    end



    max_t = data[1, end]
    min_t = data[1, 1]
    tspan = (min_t, max_t)
    tsteps = data[1, :]


    # "starting condition  using data if smoothing average is used skip this part

    if smoothing == true

        u0 = [data[2, 1]]

    else
        u0 = [Statistics.mean(data[2, 1:pt_avg])]
    end

    #  definition the possible loss (not optimal here but its working)





    function loss_ode_L2(p)
        sol = solve(ODE_prob, KenCarp4(autodiff=true), p=p, saveat=tsteps, verbose=false, abstol=1e-10, reltol=1e-10)
        sol_t = reduce(hcat, sol.u)
        if size(sol_t)[2] == size(data)[2]
            lossa = NaNMath.sum(abs2.((data[2, :] .- sol_t[1, 1:end]))) / length(data[2, :])
        else
            lossa =10.0^9

#            lossa = NaNMath.sum(abs2.((data[2, :] .- sol_t[1, 1:end]))) / length(data[2, :]) + 1000 * (size(data)[2] - size(sol_t)[2])
        end

        return lossa, sol
    end


    function loss_ode_L2_derivative(p)
        sol = solve(ODE_prob, KenCarp4(autodiff=true), p=p, saveat=tsteps, verbose=false, abstol=1e-10, reltol=1e-10)
        sol_t = reduce(hcat, sol.u)
        sol_time = reduce(hcat, sol.t)

        derivative_theo = [(sol_t[1, i] - sol_t[1, i-1]) / (sol_time[1, i] - sol_time[1, i-1]) for i in 2:(size(sol_t)[2])]

        derivative_data = [(data[2, i] - data[2, i-1]) / (data[1, i] - data[1, i-1]) for i in 2:(size(data)[2])]

        if size(derivative_theo)[1] == size(derivative_data)[1]
            lossa = NaNMath.sum(abs2.((derivative_theo[2, :] .- derivative_data[1, 1:end]))) / length(derivative_data[2, :])
        else
            lossa = 10.0^9 * length(data[2, :])

#            lossa =10.0^9

            #   lossa = NaNMath.sum(abs2.((derivative_theo[2, :] .- derivative_data[1, 1:end]))) / length(derivative_data[2, :]) + 1000
        end

        return lossa, sol
    end


    function loss_ode_RE(p)
        sol = solve(ODE_prob, KenCarp4(autodiff=true), p=p, saveat=tsteps, verbose=false, abstol=1e-10, reltol=1e-10)
        sol_t = reduce(hcat, sol.u)
        if size(sol_t)[2] == size(data)[2]
            lossa = 0.5 * NaNMath.sum(abs2.(1.0 .- (data[2, :] ./ sol_t[1, 1:end]))) / length(data[2, :])
        else
            lossa = 10.0^9 * length(data[2, :])

#            lossa = 0.5 * NaNMath.sum(abs2.(1.0 .- (data[2, 1:size(sol_t)[2]] ./ sol_t[1, 1:end]))) / length(data[2, 1:size(sol_t)[2]]) + 1000 * (size(data)[2] - size(sol_t)[2])
        end

        return lossa, sol
    end





    function loss_ode_blank_weighted_L2(p)
        sol = solve(ODE_prob, KenCarp4(autodiff=true), p=p, saveat=tsteps, verbose=false, abstol=1e-10, reltol=1e-10)
        sol_t = reduce(hcat, sol.u)
        empirical_blank_distrib = (blank_array .- mean(blank_array))
        # generation of the empirica distrib respect the mean of the noise 
        test = StatsBase.fit(Histogram, empirical_blank_distrib)
        binning_distrib = test.edges
        probabiltity_distrib = test.weights ./ sum(test.weights)
        if size(sol_t)[2] == size(data)[2]
            lossa = 0.0



            for ll in 1:size(sol_t)[2]
                dist = (data[2, ll] - sol_t[1, ll])
                index = findfirst(x -> (x > dist), binning_distrib[1])
                if (typeof(index) == Nothing || index > length(probabiltity_distrib))

                    lossa = lossa + abs2.(dist) / length(data[2, :])
                else

                    prob = probabiltity_distrib[index]


                    lossa = lossa + abs2.((1 - prob) * dist) / length(data[2, :])

                end

            end
        else

            lossa = 10.0^9 * length(data[2, :])

        end

        return lossa, sol
    end






    # selection of the model 

    if model == "hyper_gompertz"


        ODE_prob = ODEProblem(hyper_gompertz, u0, tspan, nothing)


    end

    if model == "hyper_logistic"
        ODE_prob = ODEProblem(hyper_logistic, u0, tspan, nothing)

    end
    
    if model == "ode_von_bertalanffy"

        ODE_prob = ODEProblem(ODE_von_bertalanffy, u0, tspan, param_of_ode)

    end
    if model == "bertalanffy_richards"

        ODE_prob = ODEProblem(bertalanffy_richards, u0, tspan, nothing)

    end
    if model == "triple_piecewise_bertalanffy_richards"

        ODE_prob = ODEProblem(ODE_triple_piecewise_bertalanffy_richards, u0, tspan, nothing)

    end

    if model == "gbsm_piecewise"

        ODE_prob = ODEProblem(ODE_gbsm_piecewise, u0, tspan, nothing)

    end
    if model == "logistic"

        ODE_prob = ODEProblem(logistic, u0, tspan, nothing)

    end
    if model == "dlogistic"

        ODE_prob = ODEProblem(dlogistic, u0, tspan, nothing)

    end
    if model == "exponential"

        ODE_prob = ODEProblem(ODE_exponential, u0, tspan, nothing)

    end
    if model == "gompertz"

        ODE_prob = ODEProblem(gompertz, u0, tspan, nothing)

    end

    if model == "baranyi_richards"

        ODE_prob = ODEProblem(baranyi_richards, u0, tspan, nothing)

    end
    if model == "baranyi_exp"

        ODE_prob = ODEProblem(baranyi_exp, u0, tspan, nothing)

    end
    if model == "baranyi_roberts"

        ODE_prob = ODEProblem(baranyi_roberts, u0, tspan, nothing)

    end


    if model == "huang"
        ## attention!!!!!
        u0 = [log(Statistics.mean(data[2, 1:5]))]



        function loss_ode_L2_relative_huang(p)

            sol = solve(ODE_prob, KenCarp4(autodiff=false), p=p, saveat=tsteps, verbose=false, abstol=1e-10, reltol=1e-10)
            sol_t = reduce(hcat, sol.u)
            lossa = 0.5 * NaNMath.sum(abs2.(1.0 .- (log.(data[2, :]) ./ sol_t[1, 1:end]))) / length(data[2, :])

            return lossa, sol
        end

        ODE_prob = ODEProblem(huang, u0, tspan, nothing)


    end

    if model == "piecewise_damped_logistic"


        ODE_prob = ODEProblem(ODE_piecewise_damped_logistic, u0, tspan, nothing)


    end



    if model == "Diauxic_replicator_1"


        ODE_prob = ODEProblem(ODE_Diauxic_replicator_1, u0, tspan, nothing)


    end
    if model == "Diauxic_replicator_2"


        ODE_prob = ODEProblem(ODE_Diauxic_replicator_2, u0, tspan, nothing)


    end
    if model == "Diauxic_piecewise_damped_logistic"


        ODE_prob = ODEProblem(ODE_Diauxic_piecewise_damped_logistic, u0, tspan, nothing)


    end
    if model == "triple_piecewise_damped_logistic"


        ODE_prob = ODEProblem(ODE_triple_piecewise_damped_logistic, u0, tspan, nothing)


    end


    if model == "triple_piecewise"


        ODE_prob = ODEProblem(ODE_triple_piecewise, u0, tspan, nothing)

    end
    if model == "triple_piecewise_sublinear"


        ODE_prob = ODEProblem(ODE_triple_piecewise_sub_linear, u0, tspan, nothing)

    end


    if model == "four_piecewise"


        ODE_prob = ODEProblem(ODE_four_piecewise, u0, tspan, nothing)

    end

    # solving models that are not the huang model
    param_combination =  generation_of_combination_of_IC_morris(lb_param, ub_param, N_step_morris)

    for i in 1:size(param_combination)[2]
        param = param_combination[:,i]
        if model  ∉ ["huang","HPM", "dHPM","HPM_inhibition","dHPM_inhibition","ODEs_HPM_SR", "HPM_exp"]

    # adtype = Optimization.AutoZygote()
# adtype = Optimization.AutoForwardDiff()
        adtype = ADTypes.AutoForwardDiff()

            if type_of_loss == "L2"
                optf = Optimization.OptimizationFunction((x, p) -> loss_ode_L2(x))
            end

            if type_of_loss == "RE"
                optf = Optimization.OptimizationFunction((x, p) -> loss_ode_RE(x))
            end

            if type_of_loss == "blank_weighted_L2"
                optf = Optimization.OptimizationFunction((x, p) -> loss_ode_blank_weighted_L2(x))
            end

            if type_of_loss == "L2_derivative"
                optf = Optimization.OptimizationFunction((x, p) -> loss_ode_L2_derivative(x))
            end

            optprob_const = Optimization.OptimizationProblem(optf, param, u0, lb=lb_param, ub=ub_param)
            res = Optimization.solve(optprob_const, BBO_adaptive_de_rand_1_bin_radiuslimited())


            remade_solution = solve(remake(ODE_prob, p=res.u), KenCarp4(), saveat=tsteps)
            sol_fin = reduce(hcat, remade_solution.u)

            ## error metrics plot 





            if type_of_loss == "L2"

                if size(sol_fin)[2] == size(data)[2]
                    loss = NaNMath.sum(abs2.((data[2, :] - sol_fin[1, 1:end]))) / length(data[2, :])
                else
                    loss = NaNMath.sum(abs2.((data[2, 1:size(sol_fin)[2]] - sol_fin[1, 1:end]))) / length(data[2, 1:size(sol_fin)[2]]) + 1000 * (size(data)[2] - size(sol_fin)[2])
                end
            end


            if type_of_loss == "RE"

                if size(sol_fin)[2] == size(data)[2]
                    loss = 0.5 * NaNMath.sum(abs2.(1.0 .- (data[2, :] ./ sol_fin[1, 1:end]))) / length(data[2, :])
                else
                    loss = 0.5 * NaNMath.sum(abs2.(1.0 .- (data[2, 1:size(sol_fin)[2]] ./ sol_fin[1, 1:end]))) / length(data[2, 1:size(sol_fin)[2]]) + 1000 * (size(data)[2] - size(sol_fin)[2])
                end
            end

            if type_of_loss == "L2_derivative"

                derivative_theo = [(sol_fin[1, i] - sol_fin[1, i-1]) / (data[1, i] - data[1, i-1]) for i in 2:(size(sol_fin)[2])]

                derivative_data = [(data[2, i] - data[2, i-1]) / (data[1, i] - data[1, i-1]) for i in 2:(size(data)[2])]

                if size(derivative_theo)[1] == size(derivative_data)[1]
                    loss = NaNMath.sum(abs2.((derivative_theo[2, :] .- derivative_data[1, 1:end]))) / length(derivative_data[2, :])
                else
                    loss = NaNMath.sum(abs2.((derivative_theo[2, :] .- derivative_data[1, 1:end]))) / length(derivative_data[2, :]) + 1000
                end

            end

            if type_of_loss == "blank_weighted_L2"


                empirical_blank_distrib = (blank_array .- mean(blank_array))
                # generation of the empirica distrib respect the mean of the noise 
                test = StatsBase.fit(Histogram, empirical_blank_distrib)
                binning_distrib = test.edges
                probabiltity_distrib = test.weights ./ sum(test.weights)

                if size(sol_fin)[2] == size(data)[2]
                    loss = 0.0



                    for ll in 1:size(sol_fin)[2]
                        dist = (data[2, ll] - sol_fin[1, ll])
                        index = findfirst(x -> (x > dist), binning_distrib[1])
                        if (typeof(index) == Nothing || index > length(probabiltity_distrib))

                            loss = loss + (abs2.(dist)) / length(data[2, :])
                        else
                            prob = probabiltity_distrib[index]

                            loss = loss + (abs2.((1 - prob) * dist)) / length(data[2, :])

                        end

                    end
                else

                    loss = 10.0^9 * length(data[2, :])

                end


            end

            res_temp = res.u
            max_th_gr = maximum( specific_gr_evaluation(  Matrix(transpose(hcat(data[1,:],  sol_fin[1, 1:end]))), pt_smoothing_derivative))
            # max empirical gr
            max_em_gr = maximum( specific_gr_evaluation(data, pt_smoothing_derivative))
        
            res_param = vectorize_df_results(name_well,
                model,
                res_temp,
                max_th_gr,
                max_em_gr,
                loss
            )
    
            results_sensitivity =hcat(results_sensitivity,res_param)
    


        end



        if model == "HPM_3_death" || model == "HPM_3_inhibition"|| model == "HPM_3_death_resistance"|| model == "dHPM_3_death_resistance"
            if smoothing == true
    
                u0 = [data[2, 1],0.0,0.0]
        
            else
                u0 = [Statistics.mean(data[2, 1:pt_avg]),0.0,0.0]
            end

            
            ODE_prob = model_selector(model,u0,tspan)
            
            function loss_ode_L2_3eq(p)
                sol = solve(ODE_prob, KenCarp4(autodiff=true), p=p, saveat=tsteps, verbose=false, abstol=1e-10, reltol=1e-10)
                sol_t = reduce(hcat, sol.u)
                sol_t_2 = sol_t[2,:]
                sol_t_2 =  reshape(sol_t_2, 1,length(sol_t_2))
                sol_t_3 = sol_t[3,:]
                sol_t_3 =  reshape(sol_t_3, 1,length(sol_t_3))
                sol_t_1 = sol_t[1,:]
                sol_t_1 =  reshape(sol_t_1, 1,length(sol_t_1))
    
                sol_t = sol_t_1 .+ sol_t_2 .+sol_t_3
     
                if size(sol_t)[2] == size(data)[2]
                    lossa = NaNMath.sum(abs2.((data[2, :] .- sol_t[1, 1:end]))) / length(data[2, :])
                else
                    lossa = 10.0^9 * length(data[2, :])
    
                                #                lossa = NaNMath.sum(abs2.((data[2, :] .- sol_t[1, 1:end]))) / length(data[2, :]) + 1000 * (size(data)[2] - size(sol_t)[2])
                end
            
                return lossa, sol
            end
    
    
    
            
    
            function loss_ode_RE_3eq(p)
                sol = solve(ODE_prob, KenCarp4(autodiff=true), p=p, saveat=tsteps, verbose=false, abstol=1e-10, reltol=1e-10)
    
                sol_t = reduce(hcat, sol.u)
                sol_t_2 = sol_t[2,:]
                sol_t_2 =  reshape(sol_t_2, 1,length(sol_t_2))
    
                sol_t_1 = sol_t[1,:]
                sol_t_1 =  reshape(sol_t_1, 1,length(sol_t_1))
    
    
                sol_t_3 = sol_t[3,:]
                sol_t_3 =  reshape(sol_t_3, 1,length(sol_t_3))
    
                sol_t = sol_t_1 .+ sol_t_2 .+sol_t_3
     
    
                if size(sol_t)[2] == size(data)[2]
                    lossa = 0.5 * NaNMath.sum(abs2.(1.0 .- (data[2, :] ./ sol_t[1, 1:end]))) / length(data[2, :])
                else
                    lossa = 10.0^9 * length(data[2, :])
    
    #                lossa = 0.5 * NaNMath.sum(abs2.(1.0 .- (data[2, 1:size(sol_t)[2]] ./ sol_t[1, 1:end]))) / length(data[2, 1:size(sol_t)[2]]) + 1000 * (size(data)[2] - size(sol_t)[2]) 
                end
    
                return lossa, sol
            end        
            
    
    
            function loss_ode_blank_weighted_L2_3eq(p)
                sol = solve(ODE_prob, KenCarp4(autodiff=true), p=p, saveat=tsteps, verbose=false, abstol=1e-10, reltol=1e-10)
                sol_t = reduce(hcat, sol.u)
                sol_t_2 = sol_t[2,:]
                sol_t_2 =  reshape(sol_t_2, 1,length(sol_t_2))
    
                sol_t_1 = sol_t[1,:]
                sol_t_1 =  reshape(sol_t_1, 1,length(sol_t_1))
    
                sol_t_3 = sol_t[3,:]
                sol_t_3 =  reshape(sol_t_3, 1,length(sol_t_3))
    
                sol_t = sol_t_1 .+ sol_t_2 .+sol_t_3          
                empirical_blank_distrib = (blank_array .- mean(blank_array))
                # generation of the empirica distrib respect the mean of the noise 
                test = StatsBase.fit(Histogram, empirical_blank_distrib)
                binning_distrib = test.edges
                probabiltity_distrib = test.weights ./ sum(test.weights)
                if size(sol_t)[2] == size(data)[2]
                    lossa = 0.0
            
            
            
                    for ll in 1:size(sol_t)[2]
                        dist = (data[2, ll] - sol_t[1, ll])
                        index = findfirst(x -> (x > dist), binning_distrib[1])
                        if (typeof(index) == Nothing || index > length(probabiltity_distrib))
           
                            lossa = lossa + abs2.(dist) / length(data[2, :])
                        else
            
                            prob = probabiltity_distrib[index]
            
            
                            lossa = lossa + abs2.((1 - prob) * dist) / length(data[2, :])
            
                        end
            
                    end
                    
                
                else
            
                    lossa = 10.0^9 * length(data[2, :])
            
                end
            
                return lossa, sol
            end      
    
    
            function loss_ode_L2_derivative_3eq(p)
              sol = solve(ODE_prob, KenCarp4(autodiff=true), p=p, saveat=tsteps, verbose=false, abstol=1e-10, reltol=1e-10)
    
              sol_t = reduce(hcat, sol.u)
              sol_t_2 = sol_t[2,:]
              sol_t_2 =  reshape(sol_t_2, 1,length(sol_t_2))
    
              sol_t_1 = sol_t[1,:]
              sol_t_1 =  reshape(sol_t_1, 1,length(sol_t_1))
    
              sol_t_3 = sol_t[3,:]
              sol_t_3 =  reshape(sol_t_3, 1,length(sol_t_3))
    
              sol_t = sol_t_1 .+ sol_t_2 .+sol_t_3        
                if size(sol_t)[2] == size(data)[2]
                    lossa = 0.5 * NaNMath.sum(abs2.(1.0 .- (data[2, :] ./ sol_t[1, 1:end]))) / length(data[2, :])
                else
                    lossa = 10.0^9 * length(data[2, :])
    
    #                lossa = 0.5 * NaNMath.sum(abs2.(1.0 .- (data[2, 1:size(sol_t)[2]] ./ sol_t[1, 1:end]))) / length(data[2, 1:size(sol_t)[2]]) + 1000 * (size(data)[2] - size(sol_t)[2])
                end
            
                return lossa, sol
            end      
            
    
            
    # adtype = Optimization.AutoZygote()
# adtype = Optimization.AutoForwardDiff()
        adtype = ADTypes.AutoForwardDiff()

            if type_of_loss == "L2"
                optf = Optimization.OptimizationFunction((x, p) -> loss_ode_L2_3eq(x))
            end
    
            if type_of_loss == "RE"
                optf = Optimization.OptimizationFunction((x, p) -> loss_ode_RE_3eq(x))
            end
    
            if type_of_loss == "blank_weighted_L2"
                optf = Optimization.OptimizationFunction((x, p) -> loss_ode_blank_weighted_L2_3eq(x))
            end
    
            if type_of_loss == "L2_derivative"
                optf = Optimization.OptimizationFunction((x, p) -> loss_ode_L2_derivative_3eq(x))
            end
            optprob_const = Optimization.OptimizationProblem(optf, param, u0, lb=lb_param, ub=ub_param)
            res = Optimization.solve(optprob_const, BBO_adaptive_de_rand_1_bin_radiuslimited())
            
    
    
        
        
        end



        if model == "HPM" || model == "dHPM" || model == "HPM_inhibition" || model == "dHPM_inhibition"  || model=="ODEs_HPM_SR" || model == "HPM_exp"





      
            # specific initial condition for system of ODEs
            # all the biomass starts as dormient
        
        
                
                if smoothing == true
        
                    u0 = [data[2, 1],0.0]
            
                else
                    u0 = [Statistics.mean(data[2, 1:pt_avg]),0.0]
                end
        
                if model == "HPM"
                    
                    ODE_prob = ODEProblem(ODEs_McKellar, u0, tspan, nothing)
                end
                if model == "HPM_exp"
                    
                    ODE_prob = ODEProblem(ODEs_HPM_exp, u0, tspan, nothing)
                end
                if model == "dHPM"
                    
                    ODE_prob = ODEProblem(ODEs_damped_McKellar, u0, tspan, nothing)
                end
                if model == "dHPM_inhibition"
                    
                    ODE_prob = ODEProblem(ODEs_dHPM_inhibition, u0, tspan, nothing)
                end

                if model == "HPM_inhibition"
                    
                    ODE_prob = ODEProblem(ODEs_HPM_inhibition, u0, tspan, nothing)
                end
                if model == "ODEs_HPM_SR"
                    
                    ODE_prob = ODEProblem(ODEs_HPM_SR, u0, tspan, nothing)
                end
                    # defining the loss functions for this case
        
        
        
        
                function loss_ode_L2_2eq(p)
                    sol = solve(ODE_prob, KenCarp4(autodiff=true), p=p, saveat=tsteps, verbose=false, abstol=1e-10, reltol=1e-10)
                    sol_t = reduce(hcat, sol.u)
                    sol_t_2 = sol_t[2,:]
                    sol_t_2 =  reshape(sol_t_2, 1,length(sol_t_2))
        
                    sol_t_1 = sol_t[1,:]
                    sol_t_1 =  reshape(sol_t_1, 1,length(sol_t_1))
        
                    sol_t = sol_t_1 .+ sol_t_2
         
                    if size(sol_t)[2] == size(data)[2]
                        lossa = NaNMath.sum(abs2.((data[2, :] .- sol_t[1, 1:end]))) / length(data[2, :])
                    else
                        lossa = 10.0^9 * length(data[2, :])

#                        lossa = NaNMath.sum(abs2.((data[2, :] .- sol_t[1, 1:end]))) / length(data[2, :]) + 1000 * (size(data)[2] - size(sol_t)[2])
                    end
                
                    return lossa, sol
                end
        
        
        
                
        
                function loss_ode_RE_2eq(p)
                    sol = solve(ODE_prob, KenCarp4(autodiff=true), p=p, saveat=tsteps, verbose=false, abstol=1e-10, reltol=1e-10)
        
                    sol_t = reduce(hcat, sol.u)
                    sol_t_2 = sol_t[2,:]
                    sol_t_2 =  reshape(sol_t_2, 1,length(sol_t_2))
        
                    sol_t_1 = sol_t[1,:]
                    sol_t_1 =  reshape(sol_t_1, 1,length(sol_t_1))
        
                    sol_t = sol_t_1 .+ sol_t_2
         
        
                    if size(sol_t)[2] == size(data)[2]
                        lossa = 0.5 * NaNMath.sum(abs2.(1.0 .- (data[2, :] ./ sol_t[1, 1:end]))) / length(data[2, :])
                    else
                        lossa = 10.0^9 * length(data[2, :])

#                        lossa = 0.5 * NaNMath.sum(abs2.(1.0 .- (data[2, 1:size(sol_t)[2]] ./ sol_t[1, 1:end]))) / length(data[2, 1:size(sol_t)[2]]) + 1000 * (size(data)[2] - size(sol_t)[2]) 
                    end
        
                    return lossa, sol
                end        
                
        
        
                function loss_ode_blank_weighted_L2_2eq(p)
                    sol = solve(ODE_prob, KenCarp4(autodiff=true), p=p, saveat=tsteps, verbose=false, abstol=1e-10, reltol=1e-10)
                    sol_t = reduce(hcat, sol.u)
                    sol_t_2 = sol_t[2,:]
                    sol_t_2 =  reshape(sol_t_2, 1,length(sol_t_2))
        
                    sol_t_1 = sol_t[1,:]
                    sol_t_1 =  reshape(sol_t_1, 1,length(sol_t_1))
        
                    sol_t = sol_t_1 .+ sol_t_2
                  
                    empirical_blank_distrib = (blank_array .- mean(blank_array))
                    # generation of the empirica distrib respect the mean of the noise 
                    test = StatsBase.fit(Histogram, empirical_blank_distrib)
                    binning_distrib = test.edges
                    probabiltity_distrib = test.weights ./ sum(test.weights)
                    if size(sol_t)[2] == size(data)[2]
                        lossa = 0.0
                
                
                
                        for ll in 1:size(sol_t)[2]
                            dist = (data[2, ll] - sol_t[1, ll])
                            index = findfirst(x -> (x > dist), binning_distrib[1])
                            if (typeof(index) == Nothing || index > length(probabiltity_distrib))
               
                                lossa = lossa + abs2.(dist) / length(data[2, :])
                            else
                
                                prob = probabiltity_distrib[index]
                
                
                                lossa = lossa + abs2.((1 - prob) * dist) / length(data[2, :])
                
                            end
                
                        end
                        
                    
                    else
                
                        lossa = 10.0^9 * length(data[2, :])
                
                    end
                
                    return lossa, sol
                end      
        
        
                function loss_ode_L2_derivative_2eq(p)
                  sol = solve(ODE_prob, KenCarp4(autodiff=true), p=p, saveat=tsteps, verbose=false, abstol=1e-10, reltol=1e-10)
        
                  sol_t = reduce(hcat, sol.u)
                  sol_t_2 = sol_t[2,:]
                  sol_t_2 =  reshape(sol_t_2, 1,length(sol_t_2))
        
                  sol_t_1 = sol_t[1,:]
                  sol_t_1 =  reshape(sol_t_1, 1,length(sol_t_1))
        
                  sol_t = sol_t_1 .+ sol_t_2         
                
                    if size(sol_t)[2] == size(data)[2]
                        lossa = 0.5 * NaNMath.sum(abs2.(1.0 .- (data[2, :] ./ sol_t[1, 1:end]))) / length(data[2, :])
                    else
                        lossa = 10.0^9 * length(data[2, :])

                        #lossa = 0.5 * NaNMath.sum(abs2.(1.0 .- (data[2, 1:size(sol_t)[2]] ./ sol_t[1, 1:end]))) / length(data[2, 1:size(sol_t)[2]]) + 1000 * (size(data)[2] - size(sol_t)[2])
                    end
                
                    return lossa, sol
                end      
                
        
                
        # adtype = Optimization.AutoZygote()
# adtype = Optimization.AutoForwardDiff()
        adtype = ADTypes.AutoForwardDiff()

                if type_of_loss == "L2"
                    optf = Optimization.OptimizationFunction((x, p) -> loss_ode_L2_2eq(x))
                end
        
                if type_of_loss == "RE"
                    optf = Optimization.OptimizationFunction((x, p) -> loss_ode_RE_2eq(x))
                end
        
                if type_of_loss == "blank_weighted_L2"
                    optf = Optimization.OptimizationFunction((x, p) -> loss_ode_blank_weighted_L2_2eq(x))
                end
        
                if type_of_loss == "L2_derivative"
                    optf = Optimization.OptimizationFunction((x, p) -> loss_ode_L2_derivative_2eq(x))
                end
                optprob_const = Optimization.OptimizationProblem(optf, param, u0, lb=lb_param, ub=ub_param)
                res = Optimization.solve(optprob_const, BBO_adaptive_de_rand_1_bin_radiuslimited())
                
        
        
            
            res_temp = res.u


            remade_solution = solve(remake(ODE_prob, p=res.u), KenCarp4(), saveat=tsteps)
            sol_fin = reduce(hcat, remade_solution.u)
        
        
            if (model == "HPM" || model == "HPM_exp")
        
                sol_t_2 = sol_fin[2,:]
                sol_t_2 =  reshape(sol_t_2, 1,length(sol_t_2))
        
                sol_t_1 = sol_fin[1,:]
                sol_t_1 =  reshape(sol_t_1, 1,length(sol_t_1))
        
                sol_fin = sol_t_1 .+ sol_t_2
        
        
        
            end 
            if model == "ODEs_HPM_SR"
        
                sol_t_2 = sol_fin[2,:]
                sol_t_2 =  reshape(sol_t_2, 1,length(sol_t_2))
        
                sol_t_1 = sol_fin[1,:]
                sol_t_1 =  reshape(sol_t_1, 1,length(sol_t_1))
        
                sol_fin = sol_t_1 .+ sol_t_2
        
        
        
            end 
        
            if model == "dHPM"
        
                sol_t_2 = sol_fin[2,:]
                sol_t_2 =  reshape(sol_t_2, 1,length(sol_t_2))
        
                sol_t_1 = sol_fin[1,:]
                sol_t_1 =  reshape(sol_t_1, 1,length(sol_t_1))
        
                sol_fin = sol_t_1 .+ sol_t_2
        
            end 
            if model == "dHPM_inhibition"
        
                sol_t_2 = sol_fin[2,:]
                sol_t_2 =  reshape(sol_t_2, 1,length(sol_t_2))
        
                sol_t_1 = sol_fin[1,:]
                sol_t_1 =  reshape(sol_t_1, 1,length(sol_t_1))
        
                sol_fin = sol_t_1 .+ sol_t_2
        
            end 
            if model == "HPM_inhibition"
        
                sol_t_2 = sol_fin[2,:]
                sol_t_2 =  reshape(sol_t_2, 1,length(sol_t_2))
        
                sol_t_1 = sol_fin[1,:]
                sol_t_1 =  reshape(sol_t_1, 1,length(sol_t_1))
        
                sol_fin = sol_t_1 .+ sol_t_2
        
            end 
        
            if model == "HPM_3_death" || model == "HPM_3_inhibition" || model == "HPM_3_death_resistance"|| model == "dHPM_3_death_resistance"

                sol_t_2 = sol_fin[2,:]
                sol_t_2 =  reshape(sol_t_2, 1,length(sol_t_2))
        
                sol_t_1 = sol_fin[1,:]
                sol_t_1 =  reshape(sol_t_1, 1,length(sol_t_1))
        
        
                sol_t_3 = sol_fin[3,:]
                sol_t_3 =  reshape(sol_t_3, 1,length(sol_t_3))
        
                sol_fin = sol_t_3 .+ sol_t_2.+ sol_t_1
        
        
        
            end 
        
        
        
            if type_of_loss == "L2"
        
                if size(sol_fin)[2] == size(data)[2]
                    loss = NaNMath.sum(abs2.((data[2, :] - sol_fin[1, 1:end]))) / length(data[2, :])
                else
                    loss = NaNMath.sum(abs2.((data[2, 1:size(sol_fin)[2]] - sol_fin[1, 1:end]))) / length(data[2, 1:size(sol_fin)[2]]) + 1000 * (size(data)[2] - size(sol_fin)[2])
                end
            end
        
        
            if type_of_loss == "RE"
        
                if size(sol_fin)[2] == size(data)[2]
                    loss = 0.5 * NaNMath.sum(abs2.(1.0 .- (data[2, :] ./ sol_fin[1, 1:end]))) / length(data[2, :])
                else
                    loss = 0.5 * NaNMath.sum(abs2.(1.0 .- (data[2, 1:size(sol_fin)[2]] ./ sol_fin[1, 1:end]))) / length(data[2, 1:size(sol_fin)[2]]) + 1000 * (size(data)[2] - size(sol_fin)[2])
                end
            end
        
            if type_of_loss == "L2_derivative"
        
                derivative_theo = [(sol_fin[1, i] - sol_fin[1, i-1]) / (data[1, i] - data[1, i-1]) for i in 2:(size(sol_fin)[2])]
        
                derivative_data = [(data[2, i] - data[2, i-1]) / (data[1, i] - data[1, i-1]) for i in 2:(size(data)[2])]
        
                if size(derivative_theo)[1] == size(derivative_data)[1]
                    loss = NaNMath.sum(abs2.((derivative_theo[2, :] .- derivative_data[1, 1:end]))) / length(derivative_data[2, :])
                else
                    loss = NaNMath.sum(abs2.((derivative_theo[2, :] .- derivative_data[1, 1:end]))) / length(derivative_data[2, :]) + 1000
                end
        
            end
        
            if type_of_loss == "blank_weighted_L2"
        
        
                empirical_blank_distrib = (blank_array .- mean(blank_array))
                # generation of the empirica distrib respect the mean of the noise 
                test = StatsBase.fit(Histogram, empirical_blank_distrib)
                binning_distrib = test.edges
                probabiltity_distrib = test.weights ./ sum(test.weights)
        
                if size(sol_fin)[2] == size(data)[2]
                    loss = 0.0
        
        
        
                    for ll in 1:size(sol_fin)[2]
                        dist = (data[2, ll] - sol_fin[1, ll])
                        index = findfirst(x -> (x > dist), binning_distrib[1])
                        if (typeof(index) == Nothing || index > length(probabiltity_distrib))
        
                            loss = loss + (abs2.(dist)) / length(data[2, :])
                        else
                            prob = probabiltity_distrib[index]
        
                            loss = loss + (abs2.((1 - prob) * dist)) / length(data[2, :])
        
                        end
        
                    end
                else
        
                    loss = 10.0^9 * length(data[2, :])
        
                end
        
        
            end
            pt_smoothing_derivative = 7
            max_th_gr = maximum( specific_gr_evaluation(  Matrix(transpose(hcat(data[1,:],  sol_fin[1, 1:end]))), pt_smoothing_derivative))
            # max empirical gr
            max_em_gr = maximum( specific_gr_evaluation(data, pt_smoothing_derivative))
        

            res_param = vectorize_df_results(name_well,
                model,
                res_temp,
                max_th_gr,
                max_em_gr,
                loss
            )
    
            results_sensitivity =hcat(results_sensitivity,res_param)
    
        
        end





    end
                # plotting if required
    if do_plot == true
        mkpath(path_to_plot)
        display(Plots.scatter(data[1, :], data[2, :], xlabel="Time", ylabel="Arb. Units", label=["Data " nothing], markersize=2, color=:black, title=string(label_exp, " ", name_well)))
        display(Plots.plot!(remade_solution.t, sol_fin[1, 1:end], xlabel="Time", ylabel="Arb. Units", label=[string("Fitting ", model) nothing], c=:red))
        png(string(path_to_plot, label_exp, "_", model, "_", name_well, ".png"))
    end


    if write_res == true

        CSV.write(string(path_to_results, label_exp, "_results_sensitivity.csv"), Tables.table(Matrix(results_sensitivity)))
        CSV.write(string(path_to_results, label_exp, "_configuration_tested.csv"), Tables.table(Matrix(param_combination)))


    end

    return results_sensitivity

end







function stochastic_sim_with_activation(model::String, #string of the model|
    n_start::Int, # number of starting cells
    n_mol_start::Float64, # starting concentration of the limiting nutrients
    tstart::Float64, # start time of the sim
    tmax::Float64, # final time of the sim
    delta_t::Float64, # delta t for poisson approx
    k_1_val::Float64,
    k_2_val::Float64, # monod constant
    alpha_val::Float64, # massimum possible growth rate
    rate_activation::Float64, 
    lambda::Float64, # lag time
    n_mol_per_birth::Float64,# nutrient consumed per division (conc)
    volume::Float64
)


    #inizialization of times 
    tot_pop_dorm = [copy(n_start)]
    tot_pop_active = [0]
    times = [copy(tstart)]

    conc_of_nutriens = [copy(n_mol_start / volume)]
    n_times = floor((tmax - tstart) / delta_t)




    for i in 2:n_times

        # if for defining the lag phase

        if i * delta_t < lambda
            rate_per_cell = 0.0
        else

            if model == "Monod"
                factor = conc_of_nutriens[end] / (k_1_val + conc_of_nutriens[end])
            end

            if model == "Haldane"
                factor = conc_of_nutriens[end] / (k_1_val + conc_of_nutriens[end] + (conc_of_nutriens[end])^2 / k_2_val)

            end

            if model == "Blackman"
                factor = conc_of_nutriens[end] / (k_1_val)
            end



            if model == "Tesseir"
                factor = 1 - exp(conc_of_nutriens[end] * k_1_val)
            end

            if model == "Moser"
                factor = conc_of_nutriens[end]^k_2_val / (conc_of_nutriens[end]^k_2_val + k_1_val)
            end



            if model == "Aiba-Edwards"
                factor = exp(-k_2_val / conc_of_nutriens[end]) * conc_of_nutriens[end] / (conc_of_nutriens[end] + k_1_val)
            end

            if model == "Verhulst"
                factor = 1 - tot_pop[end] / k_1_val
            end

            rate_per_cell = alpha_val * factor
        end


        # evaluating the number of birth events with poisson approx
        total_activation_rate =  tot_pop_dorm[end]  * rate_activation
        total_birth_rate = rate_per_cell .* tot_pop_active[end] .* delta_t

        n_activation = rand(Distributions.Poisson(total_activation_rate), 1)
        n_birth = rand(Distributions.Poisson(total_birth_rate), 1)




        # updating states 
        new_conc = max(0, conc_of_nutriens[end] - n_birth[1] * n_mol_per_birth / volume)
        net_pop_variation = n_birth[1]
        net_activation = n_activation[1]
        new_pop_dormient = max(0, tot_pop_dorm[end] - net_activation)

        new_pop = max(0, tot_pop_active[end] + net_pop_variation + net_activation)
        conc_of_nutriens = push!(conc_of_nutriens, new_conc)
        tot_pop_active = push!(tot_pop_active, new_pop)
        tot_pop_dorm = push!(tot_pop_dorm, new_pop_dormient)

        times = push!(times, times[end] + delta_t)




    end

    return tot_pop_active,tot_pop_dorm , conc_of_nutriens, times


end

function fitting_one_well_constrained_ODE_selection(data::Matrix{Float64}, # dataset x times y OD/fluorescence
    name_well::String, # name of the well
    label_exp::String, #label of the experiment
    lb_param::Vector{Float64}, # lower bound param
    ub_param::Vector{Float64}, # upper bound param
    param::Vector{Float64}, # initial guess param
    model::String, # ode model to use 
    do_plot::Bool, # do plots or no
    path_to_plot::String, # where save plots
    pt_avg::Int, # numebr of the point to generate intial condition
    smoothing::Bool, # the smoothing is done or not?
    type_of_loss::String, # type of used loss 
    blank_array::Vector{Float64}, # data of all blanks
    error_analysis::Bool, # do error analysis
    write_res::Bool
)


    times_data = data[1, :]


    max_t = data[1, end]
    min_t = data[1, 1]
    tspan = (min_t, max_t)
    tsteps = data[1, :]


    # "starting condition  using data if smoothing average is used skip this part

    if smoothing == true

        u0 = [data[2, 1]]

    else
        u0 = [Statistics.mean(data[2, 1:pt_avg])]
    end

    #  definition the possible loss (not optimal here but its working)



    if model ∉ ["HPM","HPM_exp", "dHPM","HPM_inhibition","dHPM_inhibition","ODEs_HPM_SR","HPM_3_death","HPM_3_inhibition","HPM_3_death_resistance","dHPM_3_death_resistance"]


        function loss_ode_L2(p)
            sol = solve(ODE_prob, KenCarp4(autodiff=true), p=p, saveat=tsteps, verbose=false, abstol=1e-10, reltol=1e-10)
            sol_t = reduce(hcat, sol.u)
            if size(sol_t)[2] == size(data)[2]
                lossa = NaNMath.sum(abs2.((data[2, :] .- sol_t[1, 1:end]))) / length(data[2, :])
            else
                lossa = 10.0^9 * length(data[2, :])

#                lossa = NaNMath.sum(abs2.((data[2, :] .- sol_t[1, 1:end]))) / length(data[2, :]) + 1000 * (size(data)[2] - size(sol_t)[2])
            end

            return lossa, sol
        end


        function loss_ode_L2_derivative(p)
            sol = solve(ODE_prob, KenCarp4(autodiff=true), p=p, saveat=tsteps, verbose=false, abstol=1e-10, reltol=1e-10)
            sol_t = reduce(hcat, sol.u)
            sol_time = reduce(hcat, sol.t)

            derivative_theo = [(sol_t[1, i] - sol_t[1, i-1]) / (sol_time[1, i] - sol_time[1, i-1]) for i in 2:(size(sol_t)[2])]

            derivative_data = [(data[2, i] - data[2, i-1]) / (data[1, i] - data[1, i-1]) for i in 2:(size(data)[2])]

            if size(derivative_theo)[1] == size(derivative_data)[1]
                lossa = NaNMath.sum(abs2.((derivative_theo[2, :] .- derivative_data[1, 1:end]))) / length(derivative_data[2, :])
            else
                lossa = 10.0^9 * length(data[2, :])

                #lossa = NaNMath.sum(abs2.((derivative_theo[2, :] .- derivative_data[1, 1:end]))) / length(derivative_data[2, :]) + 1000
            end

            return lossa, sol
        end


        function loss_ode_RE(p)
            sol = solve(ODE_prob, KenCarp4(autodiff=true), p=p, saveat=tsteps, verbose=false, abstol=1e-10, reltol=1e-10)
            sol_t = reduce(hcat, sol.u)
            if size(sol_t)[2] == size(data)[2]
                lossa = 0.5 * NaNMath.sum(abs2.(1.0 .- (data[2, :] ./ sol_t[1, 1:end]))) / length(data[2, :])
            else
                lossa = 10.0^9 * length(data[2, :])

#                lossa = 0.5 * NaNMath.sum(abs2.(1.0 .- (data[2, 1:size(sol_t)[2]] ./ sol_t[1, 1:end]))) / length(data[2, 1:size(sol_t)[2]]) + 1000 * (size(data)[2] - size(sol_t)[2])
            end

            return lossa, sol
        end





        function loss_ode_blank_weighted_L2(p)
            sol = solve(ODE_prob, KenCarp4(autodiff=true), p=p, saveat=tsteps, verbose=false, abstol=1e-10, reltol=1e-10)
            sol_t = reduce(hcat, sol.u)
            empirical_blank_distrib = (blank_array .- mean(blank_array))
            # generation of the empirica distrib respect the mean of the noise 
            test = StatsBase.fit(Histogram, empirical_blank_distrib)
            binning_distrib = test.edges
            probabiltity_distrib = test.weights ./ sum(test.weights)
            if size(sol_t)[2] == size(data)[2]
                lossa = 0.0



                for ll in 1:size(sol_t)[2]
                    dist = (data[2, ll] - sol_t[1, ll])
                    index = findfirst(x -> (x > dist), binning_distrib[1])
                    if (typeof(index) == Nothing || index > length(probabiltity_distrib))

                        lossa = lossa + abs2.(dist) / length(data[2, :])
                    else

                        prob = probabiltity_distrib[index]


                        lossa = lossa + abs2.((1 - prob) * dist) / length(data[2, :])

                    end

                end
            else

                lossa = 10.0^9 * length(data[2, :])

            end

            return lossa, sol
        end



    end


    # selection of the model 

    if model == "hyper_gompertz"


        ODE_prob = ODEProblem(hyper_gompertz, u0, tspan, nothing)


    end

    if model == "hyper_logistic"
        ODE_prob = ODEProblem(hyper_logistic, u0, tspan, nothing)

    end

    if model == "bertalanffy_richards"

        ODE_prob = ODEProblem(bertalanffy_richards, u0, tspan, nothing)

    end
        
    if model == "ode_von_bertalanffy"

        ODE_prob = ODEProblem(ODE_von_bertalanffy, u0, tspan, param_of_ode)

    end
    if model == "triple_piecewise_bertalanffy_richards"

        ODE_prob = ODEProblem(ODE_triple_piecewise_bertalanffy_richards, u0, tspan, nothing)

    end

    if model == "gbsm_piecewise"

        ODE_prob = ODEProblem(ODE_gbsm_piecewise, u0, tspan, nothing)

    end
    if model == "logistic"

        ODE_prob = ODEProblem(logistic, u0, tspan, nothing)

    end
    if model == "dlogistic"

        ODE_prob = ODEProblem(dlogistic, u0, tspan, nothing)

    end
    if model == "exponential"

        ODE_prob = ODEProblem(ODE_exponential, u0, tspan, nothing)

    end
    if model == "gompertz"

        ODE_prob = ODEProblem(gompertz, u0, tspan, nothing)

    end

    if model == "baranyi_richards"

        ODE_prob = ODEProblem(baranyi_richards, u0, tspan, nothing)

    end
    if model == "baranyi_exp"

        ODE_prob = ODEProblem(baranyi_exp, u0, tspan, nothing)

    end
    if model == "baranyi_roberts"

        ODE_prob = ODEProblem(baranyi_roberts, u0, tspan, nothing)

    end


    if model == "huang"
        ## attention!!!!!
        u0 = [log(Statistics.mean(data[2, 1:pt_avg]))]



        function loss_ode_L2_relative_huang(p)

            sol = solve(ODE_prob, KenCarp4(autodiff=false), p=p, saveat=tsteps, verbose=false, abstol=1e-10, reltol=1e-10)
            sol_t = reduce(hcat, sol.u)
            lossa = 0.5 * NaNMath.sum(abs2.(1.0 .- (log.(data[2, :]) ./ sol_t[1, 1:end]))) / length(data[2, :])

            return lossa, sol
        end

        ODE_prob = ODEProblem(huang, u0, tspan, nothing)


    end

    if model == "piecewise_damped_logistic"


        ODE_prob = ODEProblem(ODE_piecewise_damped_logistic, u0, tspan, nothing)


    end



    if model == "Diauxic_replicator_1"


        ODE_prob = ODEProblem(ODE_Diauxic_replicator_1, u0, tspan, nothing)


    end
    if model == "Diauxic_replicator_2"


        ODE_prob = ODEProblem(ODE_Diauxic_replicator_2, u0, tspan, nothing)


    end
    if model == "Diauxic_piecewise_damped_logistic"


        ODE_prob = ODEProblem(ODE_Diauxic_piecewise_damped_logistic, u0, tspan, nothing)


    end
    if model == "triple_piecewise_damped_logistic"


        ODE_prob = ODEProblem(ODE_triple_piecewise_damped_logistic, u0, tspan, nothing)


    end


    if model == "triple_piecewise"


        ODE_prob = ODEProblem(ODE_triple_piecewise, u0, tspan, nothing)

    end

    if model == "triple_piecewise_sublinear"


        ODE_prob = ODEProblem(ODE_triple_piecewise_sub_linear, u0, tspan, nothing)

    end

    if model == "four_piecewise"


        ODE_prob = ODEProblem(ODE_four_piecewise, u0, tspan, nothing)

    end

    # solving models that are not the huang model


    

    if model ∉ ["huang"   , "HPM_exp","HPM", "dHPM", "HPM_inhibition", "dHPM_inhibition" ,"ODEs_HPM_SR","HPM_3_death","HPM_3_inhibition","HPM_3_death_resistance","dHPM_3_death_resistance"]


# adtype = Optimization.AutoZygote()
# adtype = Optimization.AutoForwardDiff()
        adtype = ADTypes.AutoForwardDiff()

        

        if type_of_loss == "L2"
            optf = Optimization.OptimizationFunction((x, p) -> loss_ode_L2(x))
        end

        if type_of_loss == "RE"
            optf = Optimization.OptimizationFunction((x, p) -> loss_ode_RE(x))
        end

        if type_of_loss == "blank_weighted_L2"
            optf = Optimization.OptimizationFunction((x, p) -> loss_ode_blank_weighted_L2(x))
        end

        if type_of_loss == "L2_derivative"
            optf = Optimization.OptimizationFunction((x, p) -> loss_ode_L2_derivative(x))
        end

        optprob_const = Optimization.OptimizationProblem(optf, param, u0, lb=lb_param, ub=ub_param)
        res = Optimization.solve(optprob_const, BBO_adaptive_de_rand_1_bin_radiuslimited())


    end
    
    ## evaluation for HPM 
    if model == "HPM_3_death" || model == "HPM_3_inhibition" || model == "HPM_3_death_resistance"|| model == "dHPM_3_death_resistance"
        if smoothing == true

            u0 = [data[2, 1],0.0,0.0]
    
        else
            u0 = [Statistics.mean(data[2, 1:pt_avg]),0.0,0.0]
        end
  
        if model == "HPM_3_death"


            ODE_prob = ODEProblem(HPM_3_death, u0, tspan, nothing)
    
        end
          
        if model == "HPM_3_inhibition"
    
    
            ODE_prob = ODEProblem(HPM_3_inhibition, u0, tspan, nothing)
        end


        if model == "HPM_3_death_resistance"
    
    
            ODE_prob = ODEProblem(HPM_3_death_resistance, u0, tspan, nothing)
        end


        if model == "dHPM_3_death_resistance"
    
    
            ODE_prob = ODEProblem(dHPM_3_death_resistance, u0, tspan, nothing)
        end
        function loss_ode_L2_3eq(p)
            sol = solve(ODE_prob, KenCarp4(autodiff=true), p=p, saveat=tsteps, verbose=false, abstol=1e-10, reltol=1e-10)
            sol_t = reduce(hcat, sol.u)
            sol_t_2 = sol_t[2,:]
            sol_t_2 =  reshape(sol_t_2, 1,length(sol_t_2))
            sol_t_3 = sol_t[3,:]
            sol_t_3 =  reshape(sol_t_3, 1,length(sol_t_3))
            sol_t_1 = sol_t[1,:]
            sol_t_1 =  reshape(sol_t_1, 1,length(sol_t_1))

            sol_t = sol_t_1 .+ sol_t_2 .+sol_t_3
 
            if size(sol_t)[2] == size(data)[2]
                lossa = NaNMath.sum(abs2.((data[2, :] .- sol_t[1, 1:end]))) / length(data[2, :])
            else
                lossa = 10.0^9 * length(data[2, :])

                            #                lossa = NaNMath.sum(abs2.((data[2, :] .- sol_t[1, 1:end]))) / length(data[2, :]) + 1000 * (size(data)[2] - size(sol_t)[2])
            end
        
            return lossa, sol
        end



        

        function loss_ode_RE_3eq(p)
            sol = solve(ODE_prob, KenCarp4(autodiff=true), p=p, saveat=tsteps, verbose=false, abstol=1e-10, reltol=1e-10)

            sol_t = reduce(hcat, sol.u)
            sol_t_2 = sol_t[2,:]
            sol_t_2 =  reshape(sol_t_2, 1,length(sol_t_2))

            sol_t_1 = sol_t[1,:]
            sol_t_1 =  reshape(sol_t_1, 1,length(sol_t_1))


            sol_t_3 = sol_t[3,:]
            sol_t_3 =  reshape(sol_t_3, 1,length(sol_t_3))

            sol_t = sol_t_1 .+ sol_t_2 .+sol_t_3
 

            if size(sol_t)[2] == size(data)[2]
                lossa = 0.5 * NaNMath.sum(abs2.(1.0 .- (data[2, :] ./ sol_t[1, 1:end]))) / length(data[2, :])
            else
                lossa = 10.0^9 * length(data[2, :])

#                lossa = 0.5 * NaNMath.sum(abs2.(1.0 .- (data[2, 1:size(sol_t)[2]] ./ sol_t[1, 1:end]))) / length(data[2, 1:size(sol_t)[2]]) + 1000 * (size(data)[2] - size(sol_t)[2]) 
            end

            return lossa, sol
        end        
        


        function loss_ode_blank_weighted_L2_3eq(p)
            sol = solve(ODE_prob, KenCarp4(autodiff=true), p=p, saveat=tsteps, verbose=false, abstol=1e-10, reltol=1e-10)
            sol_t = reduce(hcat, sol.u)
            sol_t_2 = sol_t[2,:]
            sol_t_2 =  reshape(sol_t_2, 1,length(sol_t_2))

            sol_t_1 = sol_t[1,:]
            sol_t_1 =  reshape(sol_t_1, 1,length(sol_t_1))

            sol_t_3 = sol_t[3,:]
            sol_t_3 =  reshape(sol_t_3, 1,length(sol_t_3))

            sol_t = sol_t_1 .+ sol_t_2 .+sol_t_3          
            empirical_blank_distrib = (blank_array .- mean(blank_array))
            # generation of the empirica distrib respect the mean of the noise 
            test = StatsBase.fit(Histogram, empirical_blank_distrib)
            binning_distrib = test.edges
            probabiltity_distrib = test.weights ./ sum(test.weights)
            if size(sol_t)[2] == size(data)[2]
                lossa = 0.0
        
        
        
                for ll in 1:size(sol_t)[2]
                    dist = (data[2, ll] - sol_t[1, ll])
                    index = findfirst(x -> (x > dist), binning_distrib[1])
                    if (typeof(index) == Nothing || index > length(probabiltity_distrib))
       
                        lossa = lossa + abs2.(dist) / length(data[2, :])
                    else
        
                        prob = probabiltity_distrib[index]
        
        
                        lossa = lossa + abs2.((1 - prob) * dist) / length(data[2, :])
        
                    end
        
                end
                
            
            else
        
                lossa = 10.0^9 * length(data[2, :])
        
            end
        
            return lossa, sol
        end      


        function loss_ode_L2_derivative_3eq(p)
          sol = solve(ODE_prob, KenCarp4(autodiff=true), p=p, saveat=tsteps, verbose=false, abstol=1e-10, reltol=1e-10)

          sol_t = reduce(hcat, sol.u)
          sol_t_2 = sol_t[2,:]
          sol_t_2 =  reshape(sol_t_2, 1,length(sol_t_2))

          sol_t_1 = sol_t[1,:]
          sol_t_1 =  reshape(sol_t_1, 1,length(sol_t_1))

          sol_t_3 = sol_t[3,:]
          sol_t_3 =  reshape(sol_t_3, 1,length(sol_t_3))

          sol_t = sol_t_1 .+ sol_t_2 .+sol_t_3        
            if size(sol_t)[2] == size(data)[2]
                lossa = 0.5 * NaNMath.sum(abs2.(1.0 .- (data[2, :] ./ sol_t[1, 1:end]))) / length(data[2, :])
            else
                lossa = 10.0^9 * length(data[2, :])

#                lossa = 0.5 * NaNMath.sum(abs2.(1.0 .- (data[2, 1:size(sol_t)[2]] ./ sol_t[1, 1:end]))) / length(data[2, 1:size(sol_t)[2]]) + 1000 * (size(data)[2] - size(sol_t)[2])
            end
        
            return lossa, sol
        end      
        

        
# adtype = Optimization.AutoZygote()
# adtype = Optimization.AutoForwardDiff()
        adtype = ADTypes.AutoForwardDiff()

        if type_of_loss == "L2"
            optf = Optimization.OptimizationFunction((x, p) -> loss_ode_L2_3eq(x))
        end

        if type_of_loss == "RE"
            optf = Optimization.OptimizationFunction((x, p) -> loss_ode_RE_3eq(x))
        end

        if type_of_loss == "blank_weighted_L2"
            optf = Optimization.OptimizationFunction((x, p) -> loss_ode_blank_weighted_L2_3eq(x))
        end

        if type_of_loss == "L2_derivative"
            optf = Optimization.OptimizationFunction((x, p) -> loss_ode_L2_derivative_3eq(x))
        end

        optprob_const = Optimization.OptimizationProblem(optf, param, u0, lb=lb_param, ub=ub_param)
        res = Optimization.solve(optprob_const, BBO_adaptive_de_rand_1_bin_radiuslimited())
        


    
    
    end


    if model == "HPM" ||model == "HPM_exp" ||    model == "dHPM" || model == "HPM_inhibition" ||  model == "dHPM_inhibition" || model == "ODEs_HPM_SR"
        # specific initial condition for system of ODEs
        # all the biomass starts as dormient


        
        if smoothing == true

            u0 = [data[2, 1],0.0]
    
        else
            u0 = [Statistics.mean(data[2, 1:pt_avg]),0.0]
        end

        if model == "HPM"
            
            ODE_prob = ODEProblem(ODEs_McKellar, u0, tspan, nothing)
        end
        
        if model == "HPM_exp"
            
            ODE_prob = ODEProblem(ODEs_HPM_exp, u0, tspan, nothing)
        end
        if model == "ODEs_HPM_SR"
            
            ODE_prob = ODEProblem(ODEs_HPM_SR, u0, tspan, nothing)
        end

        if model == "dHPM"
            
            ODE_prob = ODEProblem(ODEs_damped_McKellar, u0, tspan, nothing)
        end
        if model == "HPM_inhibition"
            
            ODE_prob = ODEProblem(ODEs_HPM_inhibition, u0, tspan, nothing)
        end

        if model == "dHPM_inhibition"
            
            ODE_prob = ODEProblem(ODEs_dHPM_inhibition, u0, tspan, nothing)
        end
            # defining the loss functions for this case




        function loss_ode_L2_2eq(p)
            sol = solve(ODE_prob, KenCarp4(autodiff=true), p=p, saveat=tsteps, verbose=false, abstol=1e-10, reltol=1e-10)
            sol_t = reduce(hcat, sol.u)
            sol_t_2 = sol_t[2,:]
            sol_t_2 =  reshape(sol_t_2, 1,length(sol_t_2))

            sol_t_1 = sol_t[1,:]
            sol_t_1 =  reshape(sol_t_1, 1,length(sol_t_1))

            sol_t = sol_t_1 .+ sol_t_2
 
            if size(sol_t)[2] == size(data)[2]
                lossa = NaNMath.sum(abs2.((data[2, :] .- sol_t[1, 1:end]))) / length(data[2, :])
            else
                lossa = 10.0^9 * length(data[2, :])

                #lossa = NaNMath.sum(abs2.((data[2, :] .- sol_t[1, 1:end]))) / length(data[2, :]) + 1000 * (size(data)[2] - size(sol_t)[2])
            end
        
            return lossa, sol
        end



        

        function loss_ode_RE_2eq(p)
            sol = solve(ODE_prob, KenCarp4(autodiff=true), p=p, saveat=tsteps, verbose=false, abstol=1e-10, reltol=1e-10)

            sol_t = reduce(hcat, sol.u)
            sol_t_2 = sol_t[2,:]
            sol_t_2 =  reshape(sol_t_2, 1,length(sol_t_2))

            sol_t_1 = sol_t[1,:]
            sol_t_1 =  reshape(sol_t_1, 1,length(sol_t_1))

            sol_t = sol_t_1 .+ sol_t_2
 

            if size(sol_t)[2] == size(data)[2]
                lossa = 0.5 * NaNMath.sum(abs2.(1.0 .- (data[2, :] ./ sol_t[1, 1:end]))) / length(data[2, :])
            else
                lossa = 10.0^9 * length(data[2, :])

                # lossa = 0.5 * NaNMath.sum(abs2.(1.0 .- (data[2, 1:size(sol_t)[2]] ./ sol_t[1, 1:end]))) / length(data[2, 1:size(sol_t)[2]]) + 1000 * (size(data)[2] - size(sol_t)[2]) 
            end

            return lossa, sol
        end        
        


        function loss_ode_blank_weighted_L2_2eq(p)
            sol = solve(ODE_prob, KenCarp4(autodiff=true), p=p, saveat=tsteps, verbose=false, abstol=1e-10, reltol=1e-10)
            sol_t = reduce(hcat, sol.u)
            sol_t_2 = sol_t[2,:]
            sol_t_2 =  reshape(sol_t_2, 1,length(sol_t_2))

            sol_t_1 = sol_t[1,:]
            sol_t_1 =  reshape(sol_t_1, 1,length(sol_t_1))

            sol_t = sol_t_1 .+ sol_t_2
          
            empirical_blank_distrib = (blank_array .- mean(blank_array))
            # generation of the empirica distrib respect the mean of the noise 
            test = StatsBase.fit(Histogram, empirical_blank_distrib)
            binning_distrib = test.edges
            probabiltity_distrib = test.weights ./ sum(test.weights)
            if size(sol_t)[2] == size(data)[2]
                lossa = 0.0
        
        
        
                for ll in 1:size(sol_t)[2]
                    dist = (data[2, ll] - sol_t[1, ll])
                    index = findfirst(x -> (x > dist), binning_distrib[1])
                    if (typeof(index) == Nothing || index > length(probabiltity_distrib))
       
                        lossa = lossa + abs2.(dist) / length(data[2, :])
                    else
        
                        prob = probabiltity_distrib[index]
        
        
                        lossa = lossa + abs2.((1 - prob) * dist) / length(data[2, :])
        
                    end
        
                end
                
            
            else
        
                lossa = 10.0^9 * length(data[2, :])
        
            end
        
            return lossa, sol
        end      


        function loss_ode_L2_derivative_2eq(p)
          sol = solve(ODE_prob, KenCarp4(autodiff=true), p=p, saveat=tsteps, verbose=false, abstol=1e-10, reltol=1e-10)

          sol_t = reduce(hcat, sol.u)
          sol_t_2 = sol_t[2,:]
          sol_t_2 =  reshape(sol_t_2, 1,length(sol_t_2))

          sol_t_1 = sol_t[1,:]
          sol_t_1 =  reshape(sol_t_1, 1,length(sol_t_1))

          sol_t = sol_t_1 .+ sol_t_2         
        
            if size(sol_t)[2] == size(data)[2]
                lossa = 0.5 * NaNMath.sum(abs2.(1.0 .- (data[2, :] ./ sol_t[1, 1:end]))) / length(data[2, :])
            else
                lossa = 10.0^9 * length(data[2, :])

#                lossa = 0.5 * NaNMath.sum(abs2.(1.0 .- (data[2, 1:size(sol_t)[2]] ./ sol_t[1, 1:end]))) / length(data[2, 1:size(sol_t)[2]]) + 1000 * (size(data)[2] - size(sol_t)[2])
            end
        
            return lossa, sol
        end      
        

        
# adtype = Optimization.AutoZygote()
# adtype = Optimization.AutoForwardDiff()
        adtype = ADTypes.AutoForwardDiff()

        if type_of_loss == "L2"
            optf = Optimization.OptimizationFunction((x, p) -> loss_ode_L2_2eq(x))
        end

        if type_of_loss == "RE"
            optf = Optimization.OptimizationFunction((x, p) -> loss_ode_RE_2eq(x))
        end

        if type_of_loss == "blank_weighted_L2"
            optf = Optimization.OptimizationFunction((x, p) -> loss_ode_blank_weighted_L2_2eq(x))
        end

        if type_of_loss == "L2_derivative"
            optf = Optimization.OptimizationFunction((x, p) -> loss_ode_L2_derivative_2eq(x))
        end
        optprob_const = Optimization.OptimizationProblem(optf, param, u0, lb=lb_param, ub=ub_param)
        res = Optimization.solve(optprob_const, BBO_adaptive_de_rand_1_bin_radiuslimited())
        


    end




 
    #revalution of solution for plot an loss evaluation 


    remade_solution = solve(remake(ODE_prob, p=res.u), KenCarp4(), saveat=tsteps)
    sol_fin = reduce(hcat, remade_solution.u)

    if model == "HPM_3_death" || model == "HPM_3_inhibition" || model == "HPM_3_death_resistance"  || model == "dHPM_3_death_resistance" 

        sol_t_2 = sol_fin[2,:]
        sol_t_2 =  reshape(sol_t_2, 1,length(sol_t_2))

        sol_t_1 = sol_fin[1,:]
        sol_t_1 =  reshape(sol_t_1, 1,length(sol_t_1))


        sol_t_3 = sol_fin[3,:]
        sol_t_3 =  reshape(sol_t_3, 1,length(sol_t_3))

        sol_fin = sol_t_3 .+ sol_t_2.+ sol_t_1



    end 
    

    if (model == "HPM" ||model == "HPM_exp")

        sol_t_2 = sol_fin[2,:]
        sol_t_2 =  reshape(sol_t_2, 1,length(sol_t_2))

        sol_t_1 = sol_fin[1,:]
        sol_t_1 =  reshape(sol_t_1, 1,length(sol_t_1))

        sol_fin = sol_t_1 .+ sol_t_2



    end 
    if model == "ODEs_HPM_SR"

        sol_t_2 = sol_fin[2,:]
        sol_t_2 =  reshape(sol_t_2, 1,length(sol_t_2))

        sol_t_1 = sol_fin[1,:]
        sol_t_1 =  reshape(sol_t_1, 1,length(sol_t_1))

        sol_fin = sol_t_1 .+ sol_t_2



    end 

    if model == "dHPM"

        sol_t_2 = sol_fin[2,:]
        sol_t_2 =  reshape(sol_t_2, 1,length(sol_t_2))

        sol_t_1 = sol_fin[1,:]
        sol_t_1 =  reshape(sol_t_1, 1,length(sol_t_1))

        sol_fin = sol_t_1 .+ sol_t_2

    end 
    if model == "dHPM_inhibition"

        sol_t_2 = sol_fin[2,:]
        sol_t_2 =  reshape(sol_t_2, 1,length(sol_t_2))

        sol_t_1 = sol_fin[1,:]
        sol_t_1 =  reshape(sol_t_1, 1,length(sol_t_1))

        sol_fin = sol_t_1 .+ sol_t_2

    end 
    if model == "HPM_inhibition"

        sol_t_2 = sol_fin[2,:]
        sol_t_2 =  reshape(sol_t_2, 1,length(sol_t_2))

        sol_t_1 = sol_fin[1,:]
        sol_t_1 =  reshape(sol_t_1, 1,length(sol_t_1))

        sol_fin = sol_t_1 .+ sol_t_2

    end 
    ## error metrics plot 


    # plotting if required



    if type_of_loss == "L2"

        if size(sol_fin)[2] == size(data)[2]
            loss = NaNMath.sum(abs2.((data[2, :] - sol_fin[1, 1:end]))) / length(data[2, :])
        else
            loss = NaNMath.sum(abs2.((data[2, 1:size(sol_fin)[2]] - sol_fin[1, 1:end]))) / length(data[2, 1:size(sol_fin)[2]]) + 1000 * (size(data)[2] - size(sol_fin)[2])
        end
    end


    if type_of_loss == "RE"

        if size(sol_fin)[2] == size(data)[2]
            loss = 0.5 * NaNMath.sum(abs2.(1.0 .- (data[2, :] ./ sol_fin[1, 1:end]))) / length(data[2, :])
        else
            loss = 0.5 * NaNMath.sum(abs2.(1.0 .- (data[2, 1:size(sol_fin)[2]] ./ sol_fin[1, 1:end]))) / length(data[2, 1:size(sol_fin)[2]]) + 1000 * (size(data)[2] - size(sol_fin)[2])
        end
    end

    if type_of_loss == "L2_derivative"

        derivative_theo = [(sol_fin[1, i] - sol_fin[1, i-1]) / (data[1, i] - data[1, i-1]) for i in 2:(size(sol_fin)[2])]

        derivative_data = [(data[2, i] - data[2, i-1]) / (data[1, i] - data[1, i-1]) for i in 2:(size(data)[2])]

        if size(derivative_theo)[1] == size(derivative_data)[1]
            loss = NaNMath.sum(abs2.((derivative_theo[2, :] .- derivative_data[1, 1:end]))) / length(derivative_data[2, :])
        else
            loss = NaNMath.sum(abs2.((derivative_theo[2, :] .- derivative_data[1, 1:end]))) / length(derivative_data[2, :]) + 1000
        end

    end

    if type_of_loss == "blank_weighted_L2"


        empirical_blank_distrib = (blank_array .- mean(blank_array))
        # generation of the empirica distrib respect the mean of the noise 
        test = StatsBase.fit(Histogram, empirical_blank_distrib)
        binning_distrib = test.edges
        probabiltity_distrib = test.weights ./ sum(test.weights)

        if size(sol_fin)[2] == size(data)[2]
            loss = 0.0



            for ll in 1:size(sol_fin)[2]
                dist = (data[2, ll] - sol_fin[1, ll])
                index = findfirst(x -> (x > dist), binning_distrib[1])
                if (typeof(index) == Nothing || index > length(probabiltity_distrib))

                    loss = loss + (abs2.(dist)) / length(data[2, :])
                else
                    prob = probabiltity_distrib[index]

                    loss = loss + (abs2.((1 - prob) * dist)) / length(data[2, :])

                end

            end
        else

            loss = 10.0^9 * length(data[2, :])

        end


    end


    # parameters have different form in function of the used model



    res_temp = res.u

        pt_smoothing_derivative = 7
    max_th_gr = maximum( specific_gr_evaluation(  Matrix(transpose(hcat(data[1,:],  sol_fin[1, 1:end]))), pt_smoothing_derivative))
    # max empirical gr
    max_em_gr = maximum( specific_gr_evaluation(data, pt_smoothing_derivative))


    res_param = vectorize_df_results(name_well,
        model,
        res_temp,
        max_th_gr,
        max_em_gr,
        loss
    )

    gr_param = res_temp[1]


    if error_analysis == true
        error_metrics(
            label_exp, #label of the experiment
            data,# name of the file to analyze
            path_to_plot, # path where save plots
            sol_fin, # 
            name_well,
            gr_param,
            model,
            true, # 1 do and visulaze the plots of data
            write_res,
            correction
        )
    end
    if do_plot == true
        mkpath(path_to_plot)
        display(Plots.scatter(data[1, :], data[2, :], xlabel="Time", ylabel="Arb. Units", label=["Data " nothing], markersize=2, color=:black, title=string(label_exp, " ", name_well)))
        display(Plots.plot!(remade_solution.t, sol_fin[1, 1:end], xlabel="Time", ylabel="Arb. Units", label=[string("Fitting ", model) nothing], c=:red))
        png(string(path_to_plot, label_exp, "_", model, "_", name_well, ".png"))
    end
   # evaluation of rss

        if size(sol_fin)[2] == size(data)[2]
            rss = NaNMath.sum(abs2.((data[2, :] - sol_fin[1, 1:end]))) 
        else
            rss = NaNMath.sum(abs2.((data[2, 1:size(sol_fin)[2]] - sol_fin[1, 1:end]))) + 1000 * (size(data)[2] - size(sol_fin)[2])
        end

    return res_param, rss,res_temp

end







# model selection using aic from :"Model selection for dynamical systems via sparse regression and information criteria"

function  ODE_Model_selection(data::Matrix{Float64}, # dataset x times y OD/fluorescence
    name_well::String, # name of the well
    label_exp::String, #label of the experiment
    lb_param_array::Any, # lower bound param
    ub_param_array::Any, # upper bound param
    models_list::Vector{String}, # ode model to use 
    pt_avg::Int, # numebr of the point to generate intial condition
    smoothing::Bool, # the smoothing is done or not?
    type_of_loss::String, # type of used loss 
    blank_array::Vector{Float64}, # data of all blanks
    plot_best_model::Bool, # one wants the results of the best fit to be plotted
    path_to_plot::String,
    verbose::Bool,
    correction::Bool, # if true uses the given calibration curve to fix the data
    calibration_OD_curve ::String  #  the path to calibration curve to fix the data
    
)

    if  correction == true

     data = correction_OD_multiple_scattering(data,calibration_OD_curve)


    end

    if plot_best_model == true
        mkpath(path_to_plot)
    end



    df_res_optimization = Array{Any}(nothing, length(models_list))
    # inizialization of array of results
    rss_array = ["model", "RSS", "AIC_standard"]

    # subtraction of the array of the blank
    if smoothing == true

        data = smoothing_data(data, pt_avg)
    end
    # defining time steps of the inference

    max_t = data[1, end]
    min_t = data[1, 1]

    data = Matrix(data)

    # fitting all required models

    for  kk in 1:length(models_list)

        if verbose == true
          println(string("fitting ", models_list[kk]))
        end
        
        param = guess_param(models_list[kk], data, min_t, max_t, lb_param_array[kk],ub_param_array[kk])

        results = fitting_one_well_constrained_ODE_selection(data, # dataset x times y OD/fluorescence
        name_well, # name of the well
        label_exp, #label of the experiment
        lb_param_array[kk], # lower bound param
        ub_param_array[kk], # upper bound param
        param, # initial guess param
        models_list[kk], # ode model to use 
        false, # do plots or no
        "path_to_plot", # where save plots
        pt_avg, # numebr of the point to generate intial condition
        smoothing, # the smoothing is done or not?
        type_of_loss, # type of used loss 
        blank_array, # data of all blanks
        false, # do error analysis
        false
    )

    df_res_optimization[kk] = results[3]
    data_size = size(data)[2]
    param_number = length(param)
    results_to_be_pushed = [  models_list[kk], results[2], data_size * log(  results[2] /data_size) + 2* param_number  ]
    rss_array = hcat(rss_array, results_to_be_pushed )

    end


    AIC_array = rss_array[3,2:end]
    min_AIC = minimum(AIC_array)
    

    index_minimal_AIC_model = findfirst(item -> item == min_AIC, AIC_array) + 1


    # string of the model choosen
    model = rss_array[1,index_minimal_AIC_model]
   # param of the best model
    param_min = df_res_optimization[index_minimal_AIC_model-1]

    
   if plot_best_model == true

        if verbose == true
            println(string("the selected model is: ",model ))
            println(string("the parameters are: ",param_min ))

        end

       # defining inital conditions
       max_t = data[1, end]
       min_t = data[1, 1]
       tspan = (min_t, max_t)
       tsteps = data[1, :] 
            # specific initial condition for system of ODEs


       if smoothing == true
            u0 = [data[2, 1]]

        else
            u0 = [Statistics.mean(data[2, 1:pt_avg])]
        end
        if model == "HPM_3_death" || model == "HPM_3_inhibition"||  model == "HPM_3_death_resistance" ||  model == "dHPM_3_death_resistance" 
            if smoothing == true
    
                u0 = [data[2, 1],0.0,0.0]
        
            else
                u0 = [Statistics.mean(data[2, 1:pt_avg]),0.0,0.0]
            end

            
  
            if model == "HPM_3_death"


                ODE_prob = ODEProblem(HPM_3_death, u0, tspan, param_min)
        
            end
              
            if model == "HPM_3_inhibition"
        
        
                ODE_prob = ODEProblem(HPM_3_inhibition, u0, tspan, param_min)

            end    
     
    

            if model == "HPM_3_death_resistance"


                ODE_prob = ODEProblem(HPM_3_death_resistance, u0, tspan, param_min)
        
            end

            if model == "dHPM_3_death_resistance"


                ODE_prob = ODEProblem(dHPM_3_death_resistance, u0, tspan, param_min)
        
            end
    
    
        
        
        end

        if model == "HPM" || model == "HPM_exp" ||  model == "dHPM"|| model == "HPM_inhibition" ||  model == "dHPM_inhibition" || model== "ODEs_HPM_SR"
            # all the biomass starts as dormient
    
    
            
            if smoothing == true
    
                u0 = [data[2, 1],0.0]
        
            else
                u0 = [Statistics.mean(data[2, 1:pt_avg]),0.0]
            end
        end
       # ode problem def 

        # selection of the model 

        if model == "hyper_gompertz"


            ODE_prob = ODEProblem(hyper_gompertz, u0, tspan, param_min)


        end

        if model == "hyper_logistic"
            ODE_prob = ODEProblem(hyper_logistic, u0, tspan, param_min)

        end

        if model == "bertalanffy_richards"

            ODE_prob = ODEProblem(bertalanffy_richards, u0, tspan, param_min)

        end
            
        if model == "ode_von_bertalanffy"

            ODE_prob = ODEProblem(ODE_von_bertalanffy, u0, tspan, param_of_ode)

        end
        
        if model == "triple_piecewise_bertalanffy_richards"

            ODE_prob = ODEProblem(ODE_triple_piecewise_bertalanffy_richards, u0, tspan, param_min)
    
        end

        if model == "gbsm_piecewise"

            ODE_prob = ODEProblem(ODE_gbsm_piecewise, u0, tspan, param_min)

        end
        
        if model == "logistic"

            ODE_prob = ODEProblem(logistic, u0, tspan, param_min)

        end        
        if model == "dlogistic"

            ODE_prob = ODEProblem(dlogistic, u0, tspan, param_min)

        end
        if model == "exponential"

            ODE_prob = ODEProblem(ODE_exponential, u0, tspan, param_min)

        end
        if model == "gompertz"

            ODE_prob = ODEProblem(gompertz, u0, tspan, param_min)

        end

        if model == "baranyi_richards"

            ODE_prob = ODEProblem(baranyi_richards, u0, tspan, param_min)

        end
        if model == "baranyi_exp"

            ODE_prob = ODEProblem(baranyi_exp, u0, tspan, param_min)
    
        end
        if model == "baranyi_roberts"

            ODE_prob = ODEProblem(baranyi_roberts, u0, tspan, param_min)

        end


        if model == "huang"
            ## attention!!!!!
            u0 = [log(Statistics.mean(data[2, 1:pt_avg]))]
            ODE_prob = ODEProblem(huang, u0, tspan, param_min)


        end

        if model == "piecewise_damped_logistic"


            ODE_prob = ODEProblem(ODE_piecewise_damped_logistic, u0, tspan, param_min)


        end



        if model == "Diauxic_replicator_1"


            ODE_prob = ODEProblem(ODE_Diauxic_replicator_1, u0, tspan, param_min)


        end
    
        if model == "Diauxic_replicator_2"


            ODE_prob = ODEProblem(ODE_Diauxic_replicator_2, u0, tspan, param_min)


        end
    
        if model == "Diauxic_piecewise_damped_logistic"


            ODE_prob = ODEProblem(ODE_Diauxic_piecewise_damped_logistic, u0, tspan, param_min)


        end
    
        if model == "triple_piecewise_damped_logistic"


            ODE_prob = ODEProblem(ODE_triple_piecewise_damped_logistic, u0, tspan, param_min)


        end


        if model == "triple_piecewise"


            ODE_prob = ODEProblem(ODE_triple_piecewise, u0, tspan, param_min)

        end
        if model == "triple_piecewise_sublinear"


            ODE_prob = ODEProblem(ODE_triple_piecewise_sub_linear, u0, tspan, param_min)
    
        end
    

        if model == "four_piecewise"


            ODE_prob = ODEProblem(ODE_four_piecewise, u0, tspan, param_min)

        end


        if model == "HPM_exp"
            
            ODE_prob = ODEProblem(ODEs_HPM_exp, u0, tspan, param_min)
        end


        if model == "HPM"
            
            ODE_prob = ODEProblem(ODEs_McKellar, u0, tspan, param_min)
        end
        
        if model == "ODEs_HPM_SR"
            
            ODE_prob = ODEProblem(ODEs_McKellar, u0, tspan, param_min)
        end

        if model == "dHPM"
            
            ODE_prob = ODEProblem(ODEs_damped_McKellar, u0, tspan, param_min)
        end

        if model == "HPM_inhibition"
            
            ODE_prob = ODEProblem(ODEs_HPM_inhibition, u0, tspan, param_min)
        end

        if model == "dHPM_inhibition"
            
            ODE_prob = ODEProblem(ODEs_dHPM_inhibition, u0, tspan, param_min)
        end

        # solving again the best model
        remade_solution = solve(ODE_prob, KenCarp4(), saveat=tsteps)
        sol_fin = reduce(hcat, remade_solution.u)
        # special case for heterogeneous pop model 
        if (model == "HPM" ||model == "HPM_exp")
        
            sol_t_2 = sol_fin[2,:]
            sol_t_2 =  reshape(sol_t_2, 1,length(sol_t_2))

            sol_t_1 = sol_fin[1,:]
            sol_t_1 =  reshape(sol_t_1, 1,length(sol_t_1))

            sol_fin = sol_t_1 .+ sol_t_2



        end 
        if model == "ODEs_HPM_SR"
        
            sol_t_2 = sol_fin[2,:]
            sol_t_2 =  reshape(sol_t_2, 1,length(sol_t_2))

            sol_t_1 = sol_fin[1,:]
            sol_t_1 =  reshape(sol_t_1, 1,length(sol_t_1))

            sol_fin = sol_t_1 .+ sol_t_2



        end 

        if model == "dHPM"

            sol_t_2 = sol_fin[2,:]
            sol_t_2 =  reshape(sol_t_2, 1,length(sol_t_2))

            sol_t_1 = sol_fin[1,:]
            sol_t_1 =  reshape(sol_t_1, 1,length(sol_t_1))

            sol_fin = sol_t_1 .+ sol_t_2

        end 


        if model == "HPM_3_death" || model == "HPM_3_inhibition" || model == "HPM_3_death_resistance"|| model == "dHPM_3_death_resistance"

            sol_t_2 = sol_fin[2,:]
            sol_t_2 =  reshape(sol_t_2, 1,length(sol_t_2))

            sol_t_1 = sol_fin[1,:]
            sol_t_1 =  reshape(sol_t_1, 1,length(sol_t_1))


            sol_t_3 = sol_fin[3,:]
            sol_t_3 =  reshape(sol_t_3, 1,length(sol_t_3))

            sol_fin = sol_t_3 .+ sol_t_2.+ sol_t_1



        end 



        if model == "dHPM_inhibition"

            sol_t_2 = sol_fin[2,:]
            sol_t_2 =  reshape(sol_t_2, 1,length(sol_t_2))

            sol_t_1 = sol_fin[1,:]
            sol_t_1 =  reshape(sol_t_1, 1,length(sol_t_1))

            sol_fin = sol_t_1 .+ sol_t_2

        end 
        if model == "HPM_inhibition"

            sol_t_2 = sol_fin[2,:]
            sol_t_2 =  reshape(sol_t_2, 1,length(sol_t_2))

            sol_t_1 = sol_fin[1,:]
            sol_t_1 =  reshape(sol_t_1, 1,length(sol_t_1))

            sol_fin = sol_t_1 .+ sol_t_2

        end 
        # plotting it
        display(Plots.scatter(data[1, :], data[2, :], xlabel="Time", ylabel="Arb. Units", label=["Data " nothing], markersize=2, color=:black, title=string(label_exp, " ", name_well)))
        display(Plots.plot!(remade_solution.t, sol_fin[1, 1:end], xlabel="Time", ylabel="Arb. Units", label=[string("Fitting ", model) nothing], c=:red))
        png(string(path_to_plot, label_exp, "_", model, "_", name_well, ".png"))
    end
   return rss_array,df_res_optimization, min_AIC, minimum(rss_array[2,2:end]) ,param_min,model


end











function fitting_custom_ODE(data::Matrix{Float64}, # dataset x times y OD/fluorescence
    name_well::String, # name of the well
    label_exp::String, #label of the experiment
    lb_param::Vector{Float64}, # lower bound param
    ub_param::Vector{Float64}, # upper bound param
    param::Vector{Float64}, # initial guess param
    model::Any, # ode model as a sciml function
    do_plot::Bool, # do plots or no
    path_to_plot::String, # where save plots
    pt_avg::Int, # numebr of the point to generate intial condition
    smoothing::Bool, # the smoothing is done or not?
    type_of_loss::String, # type of used loss 
    blank_array::Vector{Float64}, # data of all blanks
    integrator::Any, # sciml integrator for ODE KenCarp4 suggested 
    optimizator::Any, # any optimizator from OptimizationBBO that support box constrains
    correction::Bool, # if true uses the given calibration curve to fix the data
    calibration_OD_curve ::String  #  the path to calibration curve to fix the data
    
)

    if  correction == true

     data = correction_OD_multiple_scattering(data,calibration_OD_curve)


    end





    max_t = data[1, end]
    min_t = data[1, 1]
    tspan = (min_t, max_t)
    tsteps = data[1, :]


    # "starting condition  using data if smoothing average is used skip this part

    if smoothing == true

        u0 = [data[2, 1]]

    else
        u0 = [Statistics.mean(data[2, 1:pt_avg])]
    end



    function loss_ode_L2(p)
        sol = solve(ODE_prob, integrator, p=p, saveat=tsteps, verbose=false, abstol=1e-10, reltol=1e-10)
        sol_t = reduce(hcat, sol.u)
        if size(sol_t)[2] == size(data)[2]
            lossa = NaNMath.sum(abs2.((data[2, :] .- sol_t[1, 1:end]))) / length(data[2, :])
        else
            lossa = 10.0^9 * length(data[2, :])

        # lossa = NaNMath.sum(abs2.((data[2, :] .- sol_t[1, 1:end]))) / length(data[2, :]) + 1000 * (size(data)[2] - size(sol_t)[2])
        end

        return lossa, sol
    end


    function loss_ode_L2_derivative(p)
        sol = solve(ODE_prob, integrator, p=p, saveat=tsteps, verbose=false, abstol=1e-10, reltol=1e-10)
        sol_t = reduce(hcat, sol.u)
        sol_time = reduce(hcat, sol.t)

        derivative_theo = [(sol_t[1, i] - sol_t[1, i-1]) / (sol_time[1, i] - sol_time[1, i-1]) for i in 2:(size(sol_t)[2])]

        derivative_data = [(data[2, i] - data[2, i-1]) / (data[1, i] - data[1, i-1]) for i in 2:(size(data)[2])]

        if size(derivative_theo)[1] == size(derivative_data)[1]
            lossa = NaNMath.sum(abs2.((derivative_theo[2, :] .- derivative_data[1, 1:end]))) / length(derivative_data[2, :])
        else
            lossa = 10.0^9 * length(data[2, :])
           # lossa = NaNMath.sum(abs2.((derivative_theo[2, :] .- derivative_data[1, 1:end]))) / length(derivative_data[2, :]) + 1000
        end

        return lossa, sol
    end


    function loss_ode_RE(p)
        sol = solve(ODE_prob, integrator, p=p, saveat=tsteps, verbose=false, abstol=1e-10, reltol=1e-10)
        sol_t = reduce(hcat, sol.u)
        if size(sol_t)[2] == size(data)[2]
            lossa = 0.5 * NaNMath.sum(abs2.(1.0 .- (data[2, :] ./ sol_t[1, 1:end]))) / length(data[2, :])
        else
            lossa = 10.0^9 * length(data[2, :])

            #lossa = 0.5 * NaNMath.sum(abs2.(1.0 .- (data[2, 1:size(sol_t)[2]] ./ sol_t[1, 1:end]))) / length(data[2, 1:size(sol_t)[2]]) + 1000 * (size(data)[2] - size(sol_t)[2])
        end

        return lossa, sol
    end





    function loss_ode_blank_weighted_L2(p)
        sol = solve(ODE_prob,integrator, p=p, saveat=tsteps, verbose=false, abstol=1e-10, reltol=1e-10)
        sol_t = reduce(hcat, sol.u)
        empirical_blank_distrib = (blank_array .- mean(blank_array))
        # generation of the empirica distrib respect the mean of the noise 
        test = StatsBase.fit(Histogram, empirical_blank_distrib)
        binning_distrib = test.edges
        probabiltity_distrib = test.weights ./ sum(test.weights)
        if size(sol_t)[2] == size(data)[2]
            lossa = 0.0



            for ll in 1:size(sol_t)[2]
                dist = (data[2, ll] - sol_t[1, ll])
                index = findfirst(x -> (x > dist), binning_distrib[1])
                if (typeof(index) == Nothing || index > length(probabiltity_distrib))

                    lossa = lossa + abs2.(dist) / length(data[2, :])
                else

                    prob = probabiltity_distrib[index]


                    lossa = lossa + abs2.((1 - prob) * dist) / length(data[2, :])

                end

            end
        else

            lossa = 10.0^9 * length(data[2, :])

        end

        return lossa, sol
    end




   

    # selection of the model 


    ODE_prob = ODEProblem(model, u0, tspan, nothing)

    adtype = Optimization.AutoZygote()

    if type_of_loss == "L2"
         optf = Optimization.OptimizationFunction((x, p) -> loss_ode_L2(x))
    end

    if type_of_loss == "RE"
        optf = Optimization.OptimizationFunction((x, p) -> loss_ode_RE(x))
    end

    if type_of_loss == "blank_weighted_L2"
        optf = Optimization.OptimizationFunction((x, p) -> loss_ode_blank_weighted_L2(x))
    end

    if type_of_loss == "L2_derivative"
        optf = Optimization.OptimizationFunction((x, p) -> loss_ode_L2_derivative(x))
    end

    optprob_const = Optimization.OptimizationProblem(optf, param, u0, lb=lb_param, ub=ub_param)

    res = Optimization.solve(optprob_const, optimizator)


        



 
    #revalution of solution for plot an loss evaluation 


    remade_solution = solve(remake(ODE_prob, p=res.u), integrator, saveat=tsteps)
    sol_fin = reduce(hcat, remade_solution.u)









    if type_of_loss == "L2"

        if size(sol_fin)[2] == size(data)[2]
            loss = NaNMath.sum(abs2.((data[2, :] - sol_fin[1, 1:end]))) / length(data[2, :])
        else
            loss = NaNMath.sum(abs2.((data[2, 1:size(sol_fin)[2]] - sol_fin[1, 1:end]))) / length(data[2, 1:size(sol_fin)[2]]) + 1000 * (size(data)[2] - size(sol_fin)[2])
        end
    end


    if type_of_loss == "RE"

        if size(sol_fin)[2] == size(data)[2]
            loss = 0.5 * NaNMath.sum(abs2.(1.0 .- (data[2, :] ./ sol_fin[1, 1:end]))) / length(data[2, :])
        else
            loss = 0.5 * NaNMath.sum(abs2.(1.0 .- (data[2, 1:size(sol_fin)[2]] ./ sol_fin[1, 1:end]))) / length(data[2, 1:size(sol_fin)[2]]) + 1000 * (size(data)[2] - size(sol_fin)[2])
        end
    end

    if type_of_loss == "L2_derivative"

        derivative_theo = [(sol_fin[1, i] - sol_fin[1, i-1]) / (data[1, i] - data[1, i-1]) for i in 2:(size(sol_fin)[2])]

        derivative_data = [(data[2, i] - data[2, i-1]) / (data[1, i] - data[1, i-1]) for i in 2:(size(data)[2])]

        if size(derivative_theo)[1] == size(derivative_data)[1]
            loss = NaNMath.sum(abs2.((derivative_theo[2, :] .- derivative_data[1, 1:end]))) / length(derivative_data[2, :])
        else
            loss = NaNMath.sum(abs2.((derivative_theo[2, :] .- derivative_data[1, 1:end]))) / length(derivative_data[2, :]) + 1000
        end

    end

    if type_of_loss == "blank_weighted_L2"


        empirical_blank_distrib = (blank_array .- mean(blank_array))
        # generation of the empirica distrib respect the mean of the noise 
        test = StatsBase.fit(Histogram, empirical_blank_distrib)
        binning_distrib = test.edges
        probabiltity_distrib = test.weights ./ sum(test.weights)

        if size(sol_fin)[2] == size(data)[2]
            loss = 0.0



            for ll in 1:size(sol_fin)[2]
                dist = (data[2, ll] - sol_fin[1, ll])
                index = findfirst(x -> (x > dist), binning_distrib[1])
                if (typeof(index) == Nothing || index > length(probabiltity_distrib))

                    loss = loss + (abs2.(dist)) / length(data[2, :])
                else
                    prob = probabiltity_distrib[index]

                    loss = loss + (abs2.((1 - prob) * dist)) / length(data[2, :])

                end

            end
        else

            loss = 10.0^9 * length(data[2, :])

        end


    end


    # parameters have different form in function of the used model

    res_temp = res.u


    res_temp =push!(res_temp,loss)

    # plotting if required
    if do_plot == true
        mkpath(path_to_plot)
        display(Plots.scatter(data[1, :], data[2, :], xlabel="Time", ylabel="Arb. Units", label=["Data " nothing], markersize=2, color=:black, title=string(label_exp, " ", name_well)))
        display(Plots.plot!(remade_solution.t, sol_fin[1, 1:end], xlabel="Time", ylabel="Arb. Units", label=[string("Fitting custom model") nothing], c=:red))
        png(string(path_to_plot, label_exp, "_custom_model_", name_well, ".png"))
    end



    return res_temp

end



function correction_OD_multiple_scattering(data::Matrix{Float64},
    calibration_curve::String)




    od_calib = CSV.File(calibration_curve)
    names_of_cols = propertynames(od_calib)
    Od_real = od_calib[names_of_cols[1]]
    od_calib_array = Matrix(transpose(hcat(Od_real, od_calib[names_of_cols[2]])))
    soterd_calib =sort!(od_calib_array, rev = false,dims=2)
    itp = interpolate(soterd_calib[1,:], soterd_calib[2,:], SteffenMonotonicInterpolation())
    extrap_spline = extrapolate( itp,0)
    corrected_data = [ extrap_spline(k) for k in data[2,:] ]
    data_fin = Matrix(transpose(hcat(data[1,:],corrected_data)))
    return data_fin


end   



function plot_data( label_exp::String, #label of the experiment
path_to_data::String, # path to the folder to analyze
path_to_annotation::String,# path to the annotation of the wells
path_to_plot::String, # path where to save Plots
display_plots::Bool, # display plots in julia or not
save_plot::Bool, # save the plot or not
overlay_plots::Bool, # true a single plot for all dataset false one plot per well
blank_subtraction::String, # string on how to use blank (NO,avg_subctraction,time_avg)
average_replicate::Bool, # if true the average between replicates 
correct_negative::String, # if "thr_correction" it put a thr on the minimum value of the data with blank subracted, if "blank_correction" uses blank distrib to impute negative values
thr_negative::Float64,  # used only if correct_negative == "thr_correction"
)


    annotation = CSV.File(string(path_to_annotation))
    names_of_annotated_df = [annotation[l][1] for l in 1:length(annotation)]
    # selcting blank wells
    properties_of_annotation = [annotation[l][2] for l in 1:length(annotation)]
    list_of_blank = names_of_annotated_df[findall(x -> x == "b", properties_of_annotation)]
    list_of_discarded = names_of_annotated_df[findall(x -> x == "X", properties_of_annotation)]
    list_of_blank = Symbol.(list_of_blank)
    list_of_blank = Symbol.(list_of_blank)



    # reading files
    dfs_data = CSV.File(path_to_data)


    # shaping df for the inference

    names_of_cols = propertynames(dfs_data)


    # excluding blank data and discarded wells 
    names_of_cols = filter!(e -> !(e in list_of_blank), names_of_cols)
    if length(list_of_discarded)>0
        names_of_cols = filter!(e -> !(e in list_of_discarded), names_of_cols)
    end

    times_data = dfs_data[names_of_cols[1]]
    if length(list_of_blank) > 0
        blank_array = reduce(vcat, [(dfs_data[k]) for k in list_of_blank])
        blank_array = convert(Vector{Float64}, blank_array)
    else
        blank_array = 0.0
    end    
    ## BLANK ANALYSIS HERE 
    if blank_subtraction == "avg_blank"

        blank_value = mean([mean(dfs_data[k]) for k in list_of_blank])

    elseif blank_subtraction == "time_blank"

        blank_value = [mean([dfs_data[k][j] for k in list_of_blank]) for j in 1:length(times_data)]

    else

        blank_value = zeros(length(times_data))

    end




    ## considering replicates
    list_replicate = unique(properties_of_annotation)
    list_replicate = filter!(e -> e != "b", list_replicate)

    if average_replicate == true
        new_data = times_data

        list_replicate = unique(properties_of_annotation)
        list_replicate = filter!(e -> e != "b", list_replicate)


        for replicate_temp in list_replicate

            names_of_replicate_temp = Symbol.(names_of_annotated_df[findall(x -> x == replicate_temp, properties_of_annotation)])
            replicate_mean = [mean([dfs_data[k][j] for k in names_of_replicate_temp]) for j in 1:length(times_data)]

            new_data = hcat(new_data, replicate_mean)

        end
        new_data = DataFrame(new_data, :auto)
        rename!(new_data, vcat(:Time, reduce(vcat, Symbol.(list_replicate))))
        names_of_cols = propertynames(new_data)
        dfs_data = new_data


    end



    # creating the folder of data if one wants to save
    if save_plot == true

        mkpath(path_to_plot)


    end   

    for well_name in names_of_cols[2:end]


        name_well = string(well_name)

        if average_replicate == true

            data_values = copy(dfs_data[!, well_name])

        else
            data_values = copy(dfs_data[well_name])
        end

        # blank subtraction 
        data_values = data_values .- blank_value
        data = Matrix(transpose(hcat(times_data, data_values)))


        if correct_negative == "thr_correction"

            data = thr_negative_correction(data, thr_negative)

        end

        if correct_negative == "blank_correction"

            data = blank_distrib_negative_correction(data, blank_array)

        end
        # save & not plot overlayed plot

        if save_plot == true  && display_plots  == false && overlay_plots == true

            if well_name == names_of_cols[2]
                Plots.plot(data[1, :], data[2, :], xlabel="Time", ylabel="Arb. Units", label=[ name_well], title=string(label_exp),legend = :outertopright)


            else
                Plots.plot!(data[1, :], data[2, :], xlabel="Time", ylabel="Arb. Units", label=[ name_well],   title=string(label_exp),legend = :outertopright)

            end
            png(string(path_to_plot, label_exp, ".png"))
        end

        # save & not plot single plot
        if save_plot == true  && display_plots  == false && overlay_plots == false
            Plots.plot(data[1, :], data[2, :], xlabel="Time", ylabel="Arb. Units", label=["Data " nothing], color=:black, title=string(label_exp, " ", name_well))
            png(string(path_to_plot, label_exp, "_", name_well, ".png"))
        end
        # not save &  plot overlayed plot

        if save_plot == false  && display_plots  == true && overlay_plots == true

            if well_name == names_of_cols[2]
                display(Plots.plot(data[1, :], data[2, :], xlabel="Time", ylabel="Arb. Units", label=[name_well],  title=string(label_exp),legend = :outertopright))


            else
                display(Plots.plot!(data[1, :], data[2, :], xlabel="Time", ylabel="Arb. Units", label=[ name_well],legend = :outertopright))

            end
       
        end
        # not save &  plot single plot

        if save_plot == false  && display_plots  == true && overlay_plots == false
            display(Plots.plot(data[1, :], data[2, :], xlabel="Time", ylabel="Arb. Units", label=["Data " nothing], markersize=2, color=:black, title=string(label_exp, " ", name_well)))
        end

        #  save &  plot overlayed plot
        if save_plot == true  && display_plots  == true && overlay_plots == true

            if well_name == names_of_cols[2]
                display(Plots.plot(data[1, :], data[2, :], xlabel="Time", ylabel="Arb. Units", label=[ name_well] , title=string(label_exp, " ") ,legend = :outertopright ) )


            else
                display(Plots.plot!(data[1, :], data[2, :], xlabel="Time", ylabel="Arb. Units", label=[ name_well],legend = :outertopright))

            end
            png(string(path_to_plot, label_exp, ".png"))
            
       
        end

        #  save &  plot single plot

        if save_plot == true  && display_plots  == true && overlay_plots == false
            display(Plots.plot(data[1, :], data[2, :], xlabel="Time", ylabel="Arb. Units", label=["Data " nothing], markersize=2, color=:black, title=string(label_exp, " ", name_well)))
            png(string(path_to_plot, label_exp, "_", name_well, ".png"))
        end

  
    end
 

end










function fitting_one_well_constrained_for_changepoint_detection(data::Matrix{Float64}, # dataset x times y OD/fluorescence
    name_well::String, # name of the well
    lb_param::Vector{Float64}, # lower bound param
    ub_param::Vector{Float64}, # upper bound param
    param::Vector{Float64}, # initial guess param
    model::String, # ode model to use 
    pt_avg::Int, # numebr of the point to generate intial condition
    smoothing::Bool, # the smoothing is done or not?
    type_of_loss::String, # type of used loss 
    blank_array::Vector{Float64}, # data of all blanks
    correction::Bool, # if true uses the given calibration curve to fix the data
    calibration_OD_curve ::String  #  the path to calibration curve to fix the data
    )


    if  correction == true

        data = correction_OD_multiple_scattering(data,calibration_OD_curve)
    
    
    end





    max_t = data[1, end]
    min_t = data[1, 1]
    tspan = (min_t, max_t)
    tsteps = data[1, :]


    # "starting condition  using data if smoothing average is used skip this part

    if smoothing == true

        u0 = [data[2, 1]]

    else
        u0 = [Statistics.mean(data[2, 1:pt_avg])]
    end

    #  definition the possible loss (not optimal here but its working)



    if model ∉ ["HPM","HPM_exp", "dHPM","HPM_inhibition","dHPM_inhibition","ODEs_HPM_SR","HPM_3_death","HPM_3_inhibition","HPM_3_death_resistance","dHPM_3_death_resistance"]
              
        ODE_prob =model_selector(model,u0,tspan)

        function loss_ode_L2(p)
            sol = solve(ODE_prob, KenCarp4(autodiff=true), p=p, saveat=tsteps, verbose=false, abstol=1e-10, reltol=1e-10)
            sol_t = reduce(hcat, sol.u)
            if size(sol_t)[2] == size(data)[2]
                lossa = NaNMath.sum(abs2.((data[2, :] .- sol_t[1, 1:end]))) / length(data[2, :])
            else
                lossa = 10.0^9 * length(data[2, :])

            # lossa = NaNMath.sum(abs2.((data[2, :] .- sol_t[1, 1:end]))) / length(data[2, :]) + 1000 * (size(data)[2] - size(sol_t)[2])
            end

            return lossa, sol
        end


        function loss_ode_L2_derivative(p)
            sol = solve(ODE_prob, KenCarp4(autodiff=true), p=p, saveat=tsteps, verbose=false, abstol=1e-10, reltol=1e-10)
            sol_t = reduce(hcat, sol.u)
            sol_time = reduce(hcat, sol.t)

            derivative_theo = [(sol_t[1, i] - sol_t[1, i-1]) / (sol_time[1, i] - sol_time[1, i-1]) for i in 2:(size(sol_t)[2])]

            derivative_data = [(data[2, i] - data[2, i-1]) / (data[1, i] - data[1, i-1]) for i in 2:(size(data)[2])]

            if size(derivative_theo)[1] == size(derivative_data)[1]
                lossa = NaNMath.sum(abs2.((derivative_theo[2, :] .- derivative_data[1, 1:end]))) / length(derivative_data[2, :])
            else
                lossa = 10.0^9 * length(data[2, :])
               # lossa = NaNMath.sum(abs2.((derivative_theo[2, :] .- derivative_data[1, 1:end]))) / length(derivative_data[2, :]) + 1000
            end

            return lossa, sol
        end


        function loss_ode_RE(p)
            sol = solve(ODE_prob, KenCarp4(autodiff=true), p=p, saveat=tsteps, verbose=false, abstol=1e-10, reltol=1e-10)
            sol_t = reduce(hcat, sol.u)
            if size(sol_t)[2] == size(data)[2]
                lossa = 0.5 * NaNMath.sum(abs2.(1.0 .- (data[2, :] ./ sol_t[1, 1:end]))) / length(data[2, :])
            else
                lossa = 10.0^9 * length(data[2, :])

                #lossa = 0.5 * NaNMath.sum(abs2.(1.0 .- (data[2, 1:size(sol_t)[2]] ./ sol_t[1, 1:end]))) / length(data[2, 1:size(sol_t)[2]]) + 1000 * (size(data)[2] - size(sol_t)[2])
            end

            return lossa, sol
        end





        function loss_ode_blank_weighted_L2(p)
            sol = solve(ODE_prob, KenCarp4(autodiff=true), p=p, saveat=tsteps, verbose=false, abstol=1e-10, reltol=1e-10)
            sol_t = reduce(hcat, sol.u)
            empirical_blank_distrib = (blank_array .- mean(blank_array))
            # generation of the empirica distrib respect the mean of the noise 
            test = StatsBase.fit(Histogram, empirical_blank_distrib)
            binning_distrib = test.edges
            probabiltity_distrib = test.weights ./ sum(test.weights)
            if size(sol_t)[2] == size(data)[2]
                lossa = 0.0



                for ll in 1:size(sol_t)[2]
                    dist = (data[2, ll] - sol_t[1, ll])
                    index = findfirst(x -> (x > dist), binning_distrib[1])
                    if (typeof(index) == Nothing || index > length(probabiltity_distrib))

                        lossa = lossa + abs2.(dist) / length(data[2, :])
                    else

                        prob = probabiltity_distrib[index]


                        lossa = lossa + abs2.((1 - prob) * dist) / length(data[2, :])

                    end

                end
            else

                lossa = 10.0^9 * length(data[2, :])

            end

            return lossa, sol
        end



    end
    
    ## function for model predicifion

   
    # solving models that are not the huang model


    

    if model ∉ ["huang"   , "HPM","HPM_exp",  "dHPM","HPM_inhibition","dHPM_inhibition", "ODEs_HPM_SR","HPM_3_death","HPM_3_inhibition","HPM_3_death_resistance","dHPM_3_death_resistance"]

        ODE_prob =model_selector(model,u0,tspan)


# adtype = Optimization.AutoZygote()
# adtype = Optimization.AutoForwardDiff()
        adtype = ADTypes.AutoForwardDiff()

        if type_of_loss == "L2"
            optf = Optimization.OptimizationFunction((x, p) -> loss_ode_L2(x))
        end

        if type_of_loss == "RE"
            optf = Optimization.OptimizationFunction((x, p) -> loss_ode_RE(x))
        end

        if type_of_loss == "blank_weighted_L2"
            optf = Optimization.OptimizationFunction((x, p) -> loss_ode_blank_weighted_L2(x))
        end

        if type_of_loss == "L2_derivative"
            optf = Optimization.OptimizationFunction((x, p) -> loss_ode_L2_derivative(x))
        end

        optprob_const = Optimization.OptimizationProblem(optf, param, u0, lb=lb_param, ub=ub_param)
        res = Optimization.solve(optprob_const, BBO_adaptive_de_rand_1_bin_radiuslimited())


    end
    
    ## evaluation for HPM 
    if model == "HPM_3_death" || model == "HPM_3_inhibition"|| model == "HPM_3_death_resistance"|| model == "dHPM_3_death_resistance"
        if smoothing == true

            u0 = [data[2, 1],0.0,0.0]
    
        else
            u0 = [Statistics.mean(data[2, 1:pt_avg]),0.0,0.0]
        end
        # function selction here
        ODE_prob =model_selector(model,u0,tspan)
        
        function loss_ode_L2_3eq(p)
            sol = solve(ODE_prob, KenCarp4(autodiff=true), p=p, saveat=tsteps, verbose=false, abstol=1e-10, reltol=1e-10)
            sol_t = reduce(hcat, sol.u)
            sol_t_2 = sol_t[2,:]
            sol_t_2 =  reshape(sol_t_2, 1,length(sol_t_2))
            sol_t_3 = sol_t[3,:]
            sol_t_3 =  reshape(sol_t_3, 1,length(sol_t_3))
            sol_t_1 = sol_t[1,:]
            sol_t_1 =  reshape(sol_t_1, 1,length(sol_t_1))

            sol_t = sol_t_1 .+ sol_t_2 .+sol_t_3
 
            if size(sol_t)[2] == size(data)[2]
                lossa = NaNMath.sum(abs2.((data[2, :] .- sol_t[1, 1:end]))) / length(data[2, :])
            else
                lossa = 10.0^9 * length(data[2, :])

                            #                lossa = NaNMath.sum(abs2.((data[2, :] .- sol_t[1, 1:end]))) / length(data[2, :]) + 1000 * (size(data)[2] - size(sol_t)[2])
            end
        
            return lossa, sol
        end



        

        function loss_ode_RE_3eq(p)
            sol = solve(ODE_prob, KenCarp4(autodiff=true), p=p, saveat=tsteps, verbose=false, abstol=1e-10, reltol=1e-10)

            sol_t = reduce(hcat, sol.u)
            sol_t_2 = sol_t[2,:]
            sol_t_2 =  reshape(sol_t_2, 1,length(sol_t_2))

            sol_t_1 = sol_t[1,:]
            sol_t_1 =  reshape(sol_t_1, 1,length(sol_t_1))


            sol_t_3 = sol_t[3,:]
            sol_t_3 =  reshape(sol_t_3, 1,length(sol_t_3))

            sol_t = sol_t_1 .+ sol_t_2 .+sol_t_3
 

            if size(sol_t)[2] == size(data)[2]
                lossa = 0.5 * NaNMath.sum(abs2.(1.0 .- (data[2, :] ./ sol_t[1, 1:end]))) / length(data[2, :])
            else
                lossa = 10.0^9 * length(data[2, :])

#                lossa = 0.5 * NaNMath.sum(abs2.(1.0 .- (data[2, 1:size(sol_t)[2]] ./ sol_t[1, 1:end]))) / length(data[2, 1:size(sol_t)[2]]) + 1000 * (size(data)[2] - size(sol_t)[2]) 
            end

            return lossa, sol
        end        
        


        function loss_ode_blank_weighted_L2_3eq(p)
            sol = solve(ODE_prob, KenCarp4(autodiff=true), p=p, saveat=tsteps, verbose=false, abstol=1e-10, reltol=1e-10)
            sol_t = reduce(hcat, sol.u)
            sol_t_2 = sol_t[2,:]
            sol_t_2 =  reshape(sol_t_2, 1,length(sol_t_2))

            sol_t_1 = sol_t[1,:]
            sol_t_1 =  reshape(sol_t_1, 1,length(sol_t_1))

            sol_t_3 = sol_t[3,:]
            sol_t_3 =  reshape(sol_t_3, 1,length(sol_t_3))

            sol_t = sol_t_1 .+ sol_t_2 .+sol_t_3          
            empirical_blank_distrib = (blank_array .- mean(blank_array))
            # generation of the empirica distrib respect the mean of the noise 
            test = StatsBase.fit(Histogram, empirical_blank_distrib)
            binning_distrib = test.edges
            probabiltity_distrib = test.weights ./ sum(test.weights)
            if size(sol_t)[2] == size(data)[2]
                lossa = 0.0
        
        
        
                for ll in 1:size(sol_t)[2]
                    dist = (data[2, ll] - sol_t[1, ll])
                    index = findfirst(x -> (x > dist), binning_distrib[1])
                    if (typeof(index) == Nothing || index > length(probabiltity_distrib))
       
                        lossa = lossa + abs2.(dist) / length(data[2, :])
                    else
        
                        prob = probabiltity_distrib[index]
        
        
                        lossa = lossa + abs2.((1 - prob) * dist) / length(data[2, :])
        
                    end
        
                end
                
            
            else
        
                lossa = 10.0^9 * length(data[2, :])
        
            end
        
            return lossa, sol
        end      


        function loss_ode_L2_derivative_3eq(p)
          sol = solve(ODE_prob, KenCarp4(autodiff=true), p=p, saveat=tsteps, verbose=false, abstol=1e-10, reltol=1e-10)

          sol_t = reduce(hcat, sol.u)
          sol_t_2 = sol_t[2,:]
          sol_t_2 =  reshape(sol_t_2, 1,length(sol_t_2))

          sol_t_1 = sol_t[1,:]
          sol_t_1 =  reshape(sol_t_1, 1,length(sol_t_1))

          sol_t_3 = sol_t[3,:]
          sol_t_3 =  reshape(sol_t_3, 1,length(sol_t_3))

          sol_t = sol_t_1 .+ sol_t_2 .+sol_t_3        
            if size(sol_t)[2] == size(data)[2]
                lossa = 0.5 * NaNMath.sum(abs2.(1.0 .- (data[2, :] ./ sol_t[1, 1:end]))) / length(data[2, :])
            else
                lossa = 10.0^9 * length(data[2, :])

#                lossa = 0.5 * NaNMath.sum(abs2.(1.0 .- (data[2, 1:size(sol_t)[2]] ./ sol_t[1, 1:end]))) / length(data[2, 1:size(sol_t)[2]]) + 1000 * (size(data)[2] - size(sol_t)[2])
            end
        
            return lossa, sol
        end      
        

        
# adtype = Optimization.AutoZygote()
# adtype = Optimization.AutoForwardDiff()
##        adtype = ADTypes.AutoForwardDiff()

        if type_of_loss == "L2"
            optf = Optimization.OptimizationFunction((x, p) -> loss_ode_L2_3eq(x))
        end

        if type_of_loss == "RE"
            optf = Optimization.OptimizationFunction((x, p) -> loss_ode_RE_3eq(x))
        end

        if type_of_loss == "blank_weighted_L2"
            optf = Optimization.OptimizationFunction((x, p) -> loss_ode_blank_weighted_L2_3eq(x))
        end

        if type_of_loss == "L2_derivative"
            optf = Optimization.OptimizationFunction((x, p) -> loss_ode_L2_derivative_3eq(x))
        end
        optprob_const = Optimization.OptimizationProblem(optf, param, u0, lb=lb_param, ub=ub_param)
        res = Optimization.solve(optprob_const, BBO_adaptive_de_rand_1_bin_radiuslimited())
        


    
    
    end

    if model == "HPM" ||  model == "HPM_exp" || model == "dHPM" ||  model == "HPM_inhibition" ||  model == "dHPM_inhibition" || model == "ODEs_HPM_SR"
        # specific initial condition for system of ODEs
        # all the biomass starts as dormient


        
        if smoothing == true

            u0 = [data[2, 1],0.0]
    
        else
            u0 = [Statistics.mean(data[2, 1:pt_avg]),0.0]
        end
         #function selection
         ODE_prob =model_selector(model,u0,tspan)


        function loss_ode_L2_2eq(p)
            sol = solve(ODE_prob, KenCarp4(autodiff=true), p=p, saveat=tsteps, verbose=false, abstol=1e-10, reltol=1e-10)
            sol_t = reduce(hcat, sol.u)
            sol_t_2 = sol_t[2,:]
            sol_t_2 =  reshape(sol_t_2, 1,length(sol_t_2))

            sol_t_1 = sol_t[1,:]
            sol_t_1 =  reshape(sol_t_1, 1,length(sol_t_1))

            sol_t = sol_t_1 .+ sol_t_2
 
            if size(sol_t)[2] == size(data)[2]
                lossa = NaNMath.sum(abs2.((data[2, :] .- sol_t[1, 1:end]))) / length(data[2, :])
            else
                lossa = 10.0^9 * length(data[2, :])

                            #                lossa = NaNMath.sum(abs2.((data[2, :] .- sol_t[1, 1:end]))) / length(data[2, :]) + 1000 * (size(data)[2] - size(sol_t)[2])
            end
        
            return lossa, sol
        end



        

        function loss_ode_RE_2eq(p)
            sol = solve(ODE_prob, KenCarp4(autodiff=true), p=p, saveat=tsteps, verbose=false, abstol=1e-10, reltol=1e-10)

            sol_t = reduce(hcat, sol.u)
            sol_t_2 = sol_t[2,:]
            sol_t_2 =  reshape(sol_t_2, 1,length(sol_t_2))

            sol_t_1 = sol_t[1,:]
            sol_t_1 =  reshape(sol_t_1, 1,length(sol_t_1))

            sol_t = sol_t_1 .+ sol_t_2
 

            if size(sol_t)[2] == size(data)[2]
                lossa = 0.5 * NaNMath.sum(abs2.(1.0 .- (data[2, :] ./ sol_t[1, 1:end]))) / length(data[2, :])
            else
                lossa = 10.0^9 * length(data[2, :])

#                lossa = 0.5 * NaNMath.sum(abs2.(1.0 .- (data[2, 1:size(sol_t)[2]] ./ sol_t[1, 1:end]))) / length(data[2, 1:size(sol_t)[2]]) + 1000 * (size(data)[2] - size(sol_t)[2]) 
            end

            return lossa, sol
        end        
        


        function loss_ode_blank_weighted_L2_2eq(p)
            sol = solve(ODE_prob, KenCarp4(autodiff=true), p=p, saveat=tsteps, verbose=false, abstol=1e-10, reltol=1e-10)
            sol_t = reduce(hcat, sol.u)
            sol_t_2 = sol_t[2,:]
            sol_t_2 =  reshape(sol_t_2, 1,length(sol_t_2))

            sol_t_1 = sol_t[1,:]
            sol_t_1 =  reshape(sol_t_1, 1,length(sol_t_1))

            sol_t = sol_t_1 .+ sol_t_2
          
            empirical_blank_distrib = (blank_array .- mean(blank_array))
            # generation of the empirica distrib respect the mean of the noise 
            test = StatsBase.fit(Histogram, empirical_blank_distrib)
            binning_distrib = test.edges
            probabiltity_distrib = test.weights ./ sum(test.weights)
            if size(sol_t)[2] == size(data)[2]
                lossa = 0.0
        
        
        
                for ll in 1:size(sol_t)[2]
                    dist = (data[2, ll] - sol_t[1, ll])
                    index = findfirst(x -> (x > dist), binning_distrib[1])
                    if (typeof(index) == Nothing || index > length(probabiltity_distrib))
       
                        lossa = lossa + abs2.(dist) / length(data[2, :])
                    else
        
                        prob = probabiltity_distrib[index]
        
        
                        lossa = lossa + abs2.((1 - prob) * dist) / length(data[2, :])
        
                    end
        
                end
                
            
            else
        
                lossa = 10.0^9 * length(data[2, :])
        
            end
        
            return lossa, sol
        end      


        function loss_ode_L2_derivative_2eq(p)
          sol = solve(ODE_prob, KenCarp4(autodiff=true), p=p, saveat=tsteps, verbose=false, abstol=1e-10, reltol=1e-10)

          sol_t = reduce(hcat, sol.u)
          sol_t_2 = sol_t[2,:]
          sol_t_2 =  reshape(sol_t_2, 1,length(sol_t_2))

          sol_t_1 = sol_t[1,:]
          sol_t_1 =  reshape(sol_t_1, 1,length(sol_t_1))

          sol_t = sol_t_1 .+ sol_t_2         
        
            if size(sol_t)[2] == size(data)[2]
                lossa = 0.5 * NaNMath.sum(abs2.(1.0 .- (data[2, :] ./ sol_t[1, 1:end]))) / length(data[2, :])
            else
                lossa = 10.0^9 * length(data[2, :])

#                lossa = 0.5 * NaNMath.sum(abs2.(1.0 .- (data[2, 1:size(sol_t)[2]] ./ sol_t[1, 1:end]))) / length(data[2, 1:size(sol_t)[2]]) + 1000 * (size(data)[2] - size(sol_t)[2])
            end
        
            return lossa, sol
        end      
        

        
# adtype = Optimization.AutoZygote()
# adtype = Optimization.AutoForwardDiff()
        #adtype = ADTypes.AutoForwardDiff()

        if type_of_loss == "L2"
            optf = Optimization.OptimizationFunction((x, p) -> loss_ode_L2_2eq(x))
        end

        if type_of_loss == "RE"
            optf = Optimization.OptimizationFunction((x, p) -> loss_ode_RE_2eq(x))
        end

        if type_of_loss == "blank_weighted_L2"
            optf = Optimization.OptimizationFunction((x, p) -> loss_ode_blank_weighted_L2_2eq(x))
        end

        if type_of_loss == "L2_derivative"
            optf = Optimization.OptimizationFunction((x, p) -> loss_ode_L2_derivative_2eq(x))
        end
        optprob_const = Optimization.OptimizationProblem(optf, param, u0, lb=lb_param, ub=ub_param)
        res = Optimization.solve(optprob_const, BBO_adaptive_de_rand_1_bin_radiuslimited())
        


    end




 
    #revalution of solution for plot an loss evaluation 


        remade_solution = solve(remake(ODE_prob, p=res.u), KenCarp4(), saveat=tsteps)
        sol_fin = reduce(hcat, remade_solution.u)


    if model == "HPM_3_death" || model == "HPM_3_inhibition"|| model == "HPM_3_death_resistance"|| model == "dHPM_3_death_resistance"

        sol_t_2 = sol_fin[2,:]
        sol_t_2 =  reshape(sol_t_2, 1,length(sol_t_2))

        sol_t_1 = sol_fin[1,:]
        sol_t_1 =  reshape(sol_t_1, 1,length(sol_t_1))


        sol_t_3 = sol_fin[3,:]
        sol_t_3 =  reshape(sol_t_3, 1,length(sol_t_3))

        sol_fin = sol_t_3 .+ sol_t_2.+ sol_t_1



    end 


    if model == "HPM" || model == "HPM_exp"

        sol_t_2 = sol_fin[2,:]
        sol_t_2 =  reshape(sol_t_2, 1,length(sol_t_2))

        sol_t_1 = sol_fin[1,:]
        sol_t_1 =  reshape(sol_t_1, 1,length(sol_t_1))

        sol_fin = sol_t_1 .+ sol_t_2



    end 
    if model == "ODEs_HPM_SR"

        sol_t_2 = sol_fin[2,:]
        sol_t_2 =  reshape(sol_t_2, 1,length(sol_t_2))

        sol_t_1 = sol_fin[1,:]
        sol_t_1 =  reshape(sol_t_1, 1,length(sol_t_1))

        sol_fin = sol_t_1 .+ sol_t_2
    



    end 

    if model == "dHPM"

        sol_t_2 = sol_fin[2,:]
        sol_t_2 =  reshape(sol_t_2, 1,length(sol_t_2))

        sol_t_1 = sol_fin[1,:]
        sol_t_1 =  reshape(sol_t_1, 1,length(sol_t_1))

        sol_fin = sol_t_1 .+ sol_t_2

    end 


    if model == "HPM_inhibition"

        sol_t_2 = sol_fin[2,:]
        sol_t_2 =  reshape(sol_t_2, 1,length(sol_t_2))

        sol_t_1 = sol_fin[1,:]
        sol_t_1 =  reshape(sol_t_1, 1,length(sol_t_1))

        sol_fin = sol_t_1 .+ sol_t_2

    end 

    if model == "dHPM_inhibition"

        sol_t_2 = sol_fin[2,:]
        sol_t_2 =  reshape(sol_t_2, 1,length(sol_t_2))

        sol_t_1 = sol_fin[1,:]
        sol_t_1 =  reshape(sol_t_1, 1,length(sol_t_1))

        sol_fin = sol_t_1 .+ sol_t_2

    end 


  






    if type_of_loss == "L2"

        if size(sol_fin)[2] == size(data)[2]
            loss = NaNMath.sum(abs2.((data[2, :] - sol_fin[1, 1:end]))) / length(data[2, :])
        else
            loss = NaNMath.sum(abs2.((data[2, 1:size(sol_fin)[2]] - sol_fin[1, 1:end]))) / length(data[2, 1:size(sol_fin)[2]]) + 1000 * (size(data)[2] - size(sol_fin)[2])
        end
    end


    if type_of_loss == "RE"

        if size(sol_fin)[2] == size(data)[2]
            loss = 0.5 * NaNMath.sum(abs2.(1.0 .- (data[2, :] ./ sol_fin[1, 1:end]))) / length(data[2, :])
        else
            loss = 0.5 * NaNMath.sum(abs2.(1.0 .- (data[2, 1:size(sol_fin)[2]] ./ sol_fin[1, 1:end]))) / length(data[2, 1:size(sol_fin)[2]]) + 1000 * (size(data)[2] - size(sol_fin)[2])
        end
    end

    if type_of_loss == "L2_derivative"

        derivative_theo = [(sol_fin[1, i] - sol_fin[1, i-1]) / (data[1, i] - data[1, i-1]) for i in 2:(size(sol_fin)[2])]

        derivative_data = [(data[2, i] - data[2, i-1]) / (data[1, i] - data[1, i-1]) for i in 2:(size(data)[2])]

        if size(derivative_theo)[1] == size(derivative_data)[1]
            loss = NaNMath.sum(abs2.((derivative_theo[2, :] .- derivative_data[1, 1:end]))) / length(derivative_data[2, :])
        else
            loss = NaNMath.sum(abs2.((derivative_theo[2, :] .- derivative_data[1, 1:end]))) / length(derivative_data[2, :]) + 1000
        end

    end

    if type_of_loss == "blank_weighted_L2"


        empirical_blank_distrib = (blank_array .- mean(blank_array))
        # generation of the empirica distrib respect the mean of the noise 
        test = StatsBase.fit(Histogram, empirical_blank_distrib)
        binning_distrib = test.edges
        probabiltity_distrib = test.weights ./ sum(test.weights)

        if size(sol_fin)[2] == size(data)[2]
            loss = 0.0



            for ll in 1:size(sol_fin)[2]
                dist = (data[2, ll] - sol_fin[1, ll])
                index = findfirst(x -> (x > dist), binning_distrib[1])
                if (typeof(index) == Nothing || index > length(probabiltity_distrib))

                    loss = loss + (abs2.(dist)) / length(data[2, :])
                else
                    prob = probabiltity_distrib[index]

                    loss = loss + (abs2.((1 - prob) * dist)) / length(data[2, :])

                end

            end
        else

            loss = 10.0^9 * length(data[2, :])

        end


    end


    # parameters have different form in function of the used model
   
    #max_theoretical gr
    max_th_gr =0.0 #maximum( specific_gr_evaluation(  Matrix(transpose(hcat(data[1,:],  sol_fin[1, 1:end]))), pt_smooth_derivative))
    # max empirical gr
    max_em_gr = maximum( specific_gr_evaluation(data, pt_smooth_derivative))

    res_temp = res.u

    res_param = vectorize_df_results(name_well,
        model,
        res_temp,
        max_th_gr,
        max_em_gr,
        loss
    )



    



    return res_param

end



function fitting_one_ode_iterate(data::Matrix{Float64}, # dataset x times y OD/fluorescence
    name_well::String, # name of the well
    lb_param::Vector{Float64}, # lower bound param
    ub_param::Vector{Float64}, # upper bound param
    param::Vector{Float64}, # initial guess param
    model::String, # ode model to use 
    pt_avg::Int, # numebr of the point to generate intial condition
    smoothing::Bool, # the smoothing is done or not?
    type_of_loss::String, # type of used loss 
    blank_array::Vector{Float64}, # data of all blanks
    correction::Bool, # if true uses the given calibration curve to fix the data
    calibration_OD_curve ::String  #  the path to calibration curve to fix the data
    )


    if  correction == true

        data = correction_OD_multiple_scattering(data,calibration_OD_curve)
    
    
    end





    max_t = data[1, end]
    min_t = data[1, 1]
    tspan = (min_t, max_t)
    tsteps = data[1, :]


    # "starting condition  using data if smoothing average is used skip this part

    if smoothing == true

        u0 = [data[2, 1]]

    else
        u0 = [Statistics.mean(data[2, 1:pt_avg])]
    end

    #  definition the possible loss (not optimal here but its working)



    if model ∉ ["HPM", "HPM_exp", "dHPM","HPM_inhibition","dHPM_inhibition","ODEs_HPM_SR","HPM_3_death","HPM_3_inhibition","HPM_3_death_resistance","dHPM_3_death_resistance"]
              
        ODE_prob =model_selector(model,u0,tspan)

        function loss_ode_L2(p)
            sol = solve(ODE_prob, KenCarp4(autodiff=true), p=p, saveat=tsteps, verbose=false, abstol=1e-10, reltol=1e-10)
            sol_t = reduce(hcat, sol.u)
            if size(sol_t)[2] == size(data)[2]
                lossa = NaNMath.sum(abs2.((data[2, :] .- sol_t[1, 1:end]))) / length(data[2, :])
            else
                lossa = 10.0^9 * length(data[2, :])

            # lossa = NaNMath.sum(abs2.((data[2, :] .- sol_t[1, 1:end]))) / length(data[2, :]) + 1000 * (size(data)[2] - size(sol_t)[2])
            end

            return lossa, sol
        end


        function loss_ode_L2_derivative(p)
            sol = solve(ODE_prob, KenCarp4(autodiff=true), p=p, saveat=tsteps, verbose=false, abstol=1e-10, reltol=1e-10)
            sol_t = reduce(hcat, sol.u)
            sol_time = reduce(hcat, sol.t)

            derivative_theo = [(sol_t[1, i] - sol_t[1, i-1]) / (sol_time[1, i] - sol_time[1, i-1]) for i in 2:(size(sol_t)[2])]

            derivative_data = [(data[2, i] - data[2, i-1]) / (data[1, i] - data[1, i-1]) for i in 2:(size(data)[2])]

            if size(derivative_theo)[1] == size(derivative_data)[1]
                lossa = NaNMath.sum(abs2.((derivative_theo[2, :] .- derivative_data[1, 1:end]))) / length(derivative_data[2, :])
            else
                lossa = 10.0^9 * length(data[2, :])
               # lossa = NaNMath.sum(abs2.((derivative_theo[2, :] .- derivative_data[1, 1:end]))) / length(derivative_data[2, :]) + 1000
            end

            return lossa, sol
        end


        function loss_ode_RE(p)
            sol = solve(ODE_prob, KenCarp4(autodiff=true), p=p, saveat=tsteps, verbose=false, abstol=1e-10, reltol=1e-10)
            sol_t = reduce(hcat, sol.u)
            if size(sol_t)[2] == size(data)[2]
                lossa = 0.5 * NaNMath.sum(abs2.(1.0 .- (data[2, :] ./ sol_t[1, 1:end]))) / length(data[2, :])
            else
                lossa = 10.0^9 * length(data[2, :])

                #lossa = 0.5 * NaNMath.sum(abs2.(1.0 .- (data[2, 1:size(sol_t)[2]] ./ sol_t[1, 1:end]))) / length(data[2, 1:size(sol_t)[2]]) + 1000 * (size(data)[2] - size(sol_t)[2])
            end

            return lossa, sol
        end





        function loss_ode_blank_weighted_L2(p)
            sol = solve(ODE_prob, KenCarp4(autodiff=true), p=p, saveat=tsteps, verbose=false, abstol=1e-10, reltol=1e-10)
            sol_t = reduce(hcat, sol.u)
            empirical_blank_distrib = (blank_array .- mean(blank_array))
            # generation of the empirica distrib respect the mean of the noise 
            test = StatsBase.fit(Histogram, empirical_blank_distrib)
            binning_distrib = test.edges
            probabiltity_distrib = test.weights ./ sum(test.weights)
            if size(sol_t)[2] == size(data)[2]
                lossa = 0.0



                for ll in 1:size(sol_t)[2]
                    dist = (data[2, ll] - sol_t[1, ll])
                    index = findfirst(x -> (x > dist), binning_distrib[1])
                    if (typeof(index) == Nothing || index > length(probabiltity_distrib))

                        lossa = lossa + abs2.(dist) / length(data[2, :])
                    else

                        prob = probabiltity_distrib[index]


                        lossa = lossa + abs2.((1 - prob) * dist) / length(data[2, :])

                    end

                end
            else

                lossa = 10.0^9 * length(data[2, :])

            end

            return lossa, sol
        end



    end
    
    ## function for model predicifion

   
    # solving models that are not the huang model


    

    if model ∉ ["huang"   , "HPM","HPM_exp",  "dHPM","HPM_inhibition","dHPM_inhibition", "ODEs_HPM_SR","HPM_3_death","HPM_3_inhibition","HPM_3_death_resistance","dHPM_3_death_resistance"]

        ODE_prob =model_selector(model,u0,tspan)



        if type_of_loss == "L2"
            optf = Optimization.OptimizationFunction((x, p) -> loss_ode_L2(x))
        end

        if type_of_loss == "RE"
            optf = Optimization.OptimizationFunction((x, p) -> loss_ode_RE(x))
        end

        if type_of_loss == "blank_weighted_L2"
            optf = Optimization.OptimizationFunction((x, p) -> loss_ode_blank_weighted_L2(x))
        end

        if type_of_loss == "L2_derivative"
            optf = Optimization.OptimizationFunction((x, p) -> loss_ode_L2_derivative(x))
        end

        optprob_const = Optimization.OptimizationProblem(optf, param, u0, lb=lb_param, ub=ub_param)
        res = Optimization.solve(optprob_const, BBO_adaptive_de_rand_1_bin_radiuslimited())


    end
    
    ## evaluation for HPM 
    if model == "HPM_3_death" || model == "HPM_3_inhibition"|| model == "HPM_3_death_resistance"|| model == "dHPM_3_death_resistance"
        if smoothing == true

            u0 = [data[2, 1],0.0,0.0]
    
        else
            u0 = [Statistics.mean(data[2, 1:pt_avg]),0.0,0.0]
        end
        # function selction here
        ODE_prob =model_selector(model,u0,tspan)
        
        function loss_ode_L2_3eq(p)
            sol = solve(ODE_prob, KenCarp4(autodiff=true), p=p, saveat=tsteps, verbose=false, abstol=1e-10, reltol=1e-10)
            sol_t = reduce(hcat, sol.u)
            sol_t_2 = sol_t[2,:]
            sol_t_2 =  reshape(sol_t_2, 1,length(sol_t_2))
            sol_t_3 = sol_t[3,:]
            sol_t_3 =  reshape(sol_t_3, 1,length(sol_t_3))
            sol_t_1 = sol_t[1,:]
            sol_t_1 =  reshape(sol_t_1, 1,length(sol_t_1))

            sol_t = sol_t_1 .+ sol_t_2 .+sol_t_3
 
            if size(sol_t)[2] == size(data)[2]
                lossa = NaNMath.sum(abs2.((data[2, :] .- sol_t[1, 1:end]))) / length(data[2, :])
            else
                lossa = 10.0^9 * length(data[2, :])

                            #                lossa = NaNMath.sum(abs2.((data[2, :] .- sol_t[1, 1:end]))) / length(data[2, :]) + 1000 * (size(data)[2] - size(sol_t)[2])
            end
        
            return lossa, sol
        end



        

        function loss_ode_RE_3eq(p)
            sol = solve(ODE_prob, KenCarp4(autodiff=true), p=p, saveat=tsteps, verbose=false, abstol=1e-10, reltol=1e-10)

            sol_t = reduce(hcat, sol.u)
            sol_t_2 = sol_t[2,:]
            sol_t_2 =  reshape(sol_t_2, 1,length(sol_t_2))

            sol_t_1 = sol_t[1,:]
            sol_t_1 =  reshape(sol_t_1, 1,length(sol_t_1))


            sol_t_3 = sol_t[3,:]
            sol_t_3 =  reshape(sol_t_3, 1,length(sol_t_3))

            sol_t = sol_t_1 .+ sol_t_2 .+sol_t_3
 

            if size(sol_t)[2] == size(data)[2]
                lossa = 0.5 * NaNMath.sum(abs2.(1.0 .- (data[2, :] ./ sol_t[1, 1:end]))) / length(data[2, :])
            else
                lossa = 10.0^9 * length(data[2, :])

            end

            return lossa, sol
        end        
        


        function loss_ode_blank_weighted_L2_3eq(p)
            sol = solve(ODE_prob, KenCarp4(autodiff=true), p=p, saveat=tsteps, verbose=false, abstol=1e-10, reltol=1e-10)
            sol_t = reduce(hcat, sol.u)
            sol_t_2 = sol_t[2,:]
            sol_t_2 =  reshape(sol_t_2, 1,length(sol_t_2))

            sol_t_1 = sol_t[1,:]
            sol_t_1 =  reshape(sol_t_1, 1,length(sol_t_1))

            sol_t_3 = sol_t[3,:]
            sol_t_3 =  reshape(sol_t_3, 1,length(sol_t_3))

            sol_t = sol_t_1 .+ sol_t_2 .+sol_t_3          
            empirical_blank_distrib = (blank_array .- mean(blank_array))
            # generation of the empirica distrib respect the mean of the noise 
            test = StatsBase.fit(Histogram, empirical_blank_distrib)
            binning_distrib = test.edges
            probabiltity_distrib = test.weights ./ sum(test.weights)
            if size(sol_t)[2] == size(data)[2]
                lossa = 0.0
        
        
        
                for ll in 1:size(sol_t)[2]
                    dist = (data[2, ll] - sol_t[1, ll])
                    index = findfirst(x -> (x > dist), binning_distrib[1])
                    if (typeof(index) == Nothing || index > length(probabiltity_distrib))
       
                        lossa = lossa + abs2.(dist) / length(data[2, :])
                    else
        
                        prob = probabiltity_distrib[index]
        
        
                        lossa = lossa + abs2.((1 - prob) * dist) / length(data[2, :])
        
                    end
        
                end
                
            
            else
        
                lossa = 10.0^9 * length(data[2, :])
        
            end
        
            return lossa, sol
        end      


        function loss_ode_L2_derivative_3eq(p)
          sol = solve(ODE_prob, KenCarp4(autodiff=true), p=p, saveat=tsteps, verbose=false, abstol=1e-10, reltol=1e-10)

          sol_t = reduce(hcat, sol.u)
          sol_t_2 = sol_t[2,:]
          sol_t_2 =  reshape(sol_t_2, 1,length(sol_t_2))

          sol_t_1 = sol_t[1,:]
          sol_t_1 =  reshape(sol_t_1, 1,length(sol_t_1))

          sol_t_3 = sol_t[3,:]
          sol_t_3 =  reshape(sol_t_3, 1,length(sol_t_3))

          sol_t = sol_t_1 .+ sol_t_2 .+sol_t_3        
            if size(sol_t)[2] == size(data)[2]
                lossa = 0.5 * NaNMath.sum(abs2.(1.0 .- (data[2, :] ./ sol_t[1, 1:end]))) / length(data[2, :])
            else
                lossa = 10.0^9 * length(data[2, :])

            end
        
            return lossa, sol
        end      
        


        if type_of_loss == "L2"
            optf = Optimization.OptimizationFunction((x, p) -> loss_ode_L2_3eq(x))
        end

        if type_of_loss == "RE"
            optf = Optimization.OptimizationFunction((x, p) -> loss_ode_RE_3eq(x))
        end

        if type_of_loss == "blank_weighted_L2"
            optf = Optimization.OptimizationFunction((x, p) -> loss_ode_blank_weighted_L2_3eq(x))
        end

        if type_of_loss == "L2_derivative"
            optf = Optimization.OptimizationFunction((x, p) -> loss_ode_L2_derivative_3eq(x))
        end
        optprob_const = Optimization.OptimizationProblem(optf, param, u0, lb=lb_param, ub=ub_param)
        res = Optimization.solve(optprob_const, BBO_adaptive_de_rand_1_bin_radiuslimited())
        


    
    
    end

    if model == "HPM" ||  model == "HPM_exp" ||model == "dHPM" ||  model == "HPM_inhibition" ||  model == "dHPM_inhibition" || model == "ODEs_HPM_SR"
        # specific initial condition for system of ODEs
        # all the biomass starts as dormient


        
        if smoothing == true

            u0 = [data[2, 1],0.0]
    
        else
            u0 = [Statistics.mean(data[2, 1:pt_avg]),0.0]
        end
         #function selection
         ODE_prob =model_selector(model,u0,tspan)


        function loss_ode_L2_2eq(p)
            sol = solve(ODE_prob, KenCarp4(autodiff=true), p=p, saveat=tsteps, verbose=false, abstol=1e-10, reltol=1e-10)
            sol_t = reduce(hcat, sol.u)
            sol_t_2 = sol_t[2,:]
            sol_t_2 =  reshape(sol_t_2, 1,length(sol_t_2))

            sol_t_1 = sol_t[1,:]
            sol_t_1 =  reshape(sol_t_1, 1,length(sol_t_1))

            sol_t = sol_t_1 .+ sol_t_2
 
            if size(sol_t)[2] == size(data)[2]
                lossa = NaNMath.sum(abs2.((data[2, :] .- sol_t[1, 1:end]))) / length(data[2, :])
            else
                lossa = 10.0^9 * length(data[2, :])

                            #                lossa = NaNMath.sum(abs2.((data[2, :] .- sol_t[1, 1:end]))) / length(data[2, :]) + 1000 * (size(data)[2] - size(sol_t)[2])
            end
        
            return lossa, sol
        end



        

        function loss_ode_RE_2eq(p)
            sol = solve(ODE_prob, KenCarp4(autodiff=true), p=p, saveat=tsteps, verbose=false, abstol=1e-10, reltol=1e-10)

            sol_t = reduce(hcat, sol.u)
            sol_t_2 = sol_t[2,:]
            sol_t_2 =  reshape(sol_t_2, 1,length(sol_t_2))

            sol_t_1 = sol_t[1,:]
            sol_t_1 =  reshape(sol_t_1, 1,length(sol_t_1))

            sol_t = sol_t_1 .+ sol_t_2
 

            if size(sol_t)[2] == size(data)[2]
                lossa = 0.5 * NaNMath.sum(abs2.(1.0 .- (data[2, :] ./ sol_t[1, 1:end]))) / length(data[2, :])
            else
                lossa = 10.0^9 * length(data[2, :])

            end

            return lossa, sol
        end        
        


        function loss_ode_blank_weighted_L2_2eq(p)
            sol = solve(ODE_prob, KenCarp4(autodiff=true), p=p, saveat=tsteps, verbose=false, abstol=1e-10, reltol=1e-10)
            sol_t = reduce(hcat, sol.u)
            sol_t_2 = sol_t[2,:]
            sol_t_2 =  reshape(sol_t_2, 1,length(sol_t_2))

            sol_t_1 = sol_t[1,:]
            sol_t_1 =  reshape(sol_t_1, 1,length(sol_t_1))

            sol_t = sol_t_1 .+ sol_t_2
          
            empirical_blank_distrib = (blank_array .- mean(blank_array))
            # generation of the empirica distrib respect the mean of the noise 
            test = StatsBase.fit(Histogram, empirical_blank_distrib)
            binning_distrib = test.edges
            probabiltity_distrib = test.weights ./ sum(test.weights)
            if size(sol_t)[2] == size(data)[2]
                lossa = 0.0
        
        
        
                for ll in 1:size(sol_t)[2]
                    dist = (data[2, ll] - sol_t[1, ll])
                    index = findfirst(x -> (x > dist), binning_distrib[1])
                    if (typeof(index) == Nothing || index > length(probabiltity_distrib))
       
                        lossa = lossa + abs2.(dist) / length(data[2, :])
                    else
        
                        prob = probabiltity_distrib[index]
        
        
                        lossa = lossa + abs2.((1 - prob) * dist) / length(data[2, :])
        
                    end
        
                end
                
            
            else
        
                lossa = 10.0^9 * length(data[2, :])
        
            end
        
            return lossa, sol
        end      


        function loss_ode_L2_derivative_2eq(p)
          sol = solve(ODE_prob, KenCarp4(autodiff=true), p=p, saveat=tsteps, verbose=false, abstol=1e-10, reltol=1e-10)

          sol_t = reduce(hcat, sol.u)
          sol_t_2 = sol_t[2,:]
          sol_t_2 =  reshape(sol_t_2, 1,length(sol_t_2))

          sol_t_1 = sol_t[1,:]
          sol_t_1 =  reshape(sol_t_1, 1,length(sol_t_1))

          sol_t = sol_t_1 .+ sol_t_2         
        
            if size(sol_t)[2] == size(data)[2]
                lossa = 0.5 * NaNMath.sum(abs2.(1.0 .- (data[2, :] ./ sol_t[1, 1:end]))) / length(data[2, :])
            else
                lossa = 10.0^9 * length(data[2, :])

            end
        
            return lossa, sol
        end      
        

        

        if type_of_loss == "L2"
            optf = Optimization.OptimizationFunction((x, p) -> loss_ode_L2_2eq(x))
        end

        if type_of_loss == "RE"
            optf = Optimization.OptimizationFunction((x, p) -> loss_ode_RE_2eq(x))
        end

        if type_of_loss == "blank_weighted_L2"
            optf = Optimization.OptimizationFunction((x, p) -> loss_ode_blank_weighted_L2_2eq(x))
        end

        if type_of_loss == "L2_derivative"
            optf = Optimization.OptimizationFunction((x, p) -> loss_ode_L2_derivative_2eq(x))
        end
        optprob_const = Optimization.OptimizationProblem(optf, param, u0, lb=lb_param, ub=ub_param)
        res = Optimization.solve(optprob_const, BBO_adaptive_de_rand_1_bin_radiuslimited())
        


    end




 
    #revalution of solution for plot an loss evaluation 


        remade_solution = solve(remake(ODE_prob, p=res.u), KenCarp4(), saveat=tsteps)
        sol_fin = reduce(hcat, remade_solution.u)


    if model == "HPM_3_death" || model == "HPM_3_inhibition"|| model == "HPM_3_death_resistance"|| model == "dHPM_3_death_resistance"

        sol_t_2 = sol_fin[2,:]
        sol_t_2 =  reshape(sol_t_2, 1,length(sol_t_2))

        sol_t_1 = sol_fin[1,:]
        sol_t_1 =  reshape(sol_t_1, 1,length(sol_t_1))


        sol_t_3 = sol_fin[3,:]
        sol_t_3 =  reshape(sol_t_3, 1,length(sol_t_3))

        sol_fin = sol_t_3 .+ sol_t_2.+ sol_t_1



    end 


    if (model == "HPM" || model == "HPM_exp")

        sol_t_2 = sol_fin[2,:]
        sol_t_2 =  reshape(sol_t_2, 1,length(sol_t_2))

        sol_t_1 = sol_fin[1,:]
        sol_t_1 =  reshape(sol_t_1, 1,length(sol_t_1))

        sol_fin = sol_t_1 .+ sol_t_2



    end 
    if model == "ODEs_HPM_SR"

        sol_t_2 = sol_fin[2,:]
        sol_t_2 =  reshape(sol_t_2, 1,length(sol_t_2))

        sol_t_1 = sol_fin[1,:]
        sol_t_1 =  reshape(sol_t_1, 1,length(sol_t_1))

        sol_fin = sol_t_1 .+ sol_t_2
    



    end 

    if model == "dHPM"

        sol_t_2 = sol_fin[2,:]
        sol_t_2 =  reshape(sol_t_2, 1,length(sol_t_2))

        sol_t_1 = sol_fin[1,:]
        sol_t_1 =  reshape(sol_t_1, 1,length(sol_t_1))

        sol_fin = sol_t_1 .+ sol_t_2

    end 


    if model == "HPM_inhibition"

        sol_t_2 = sol_fin[2,:]
        sol_t_2 =  reshape(sol_t_2, 1,length(sol_t_2))

        sol_t_1 = sol_fin[1,:]
        sol_t_1 =  reshape(sol_t_1, 1,length(sol_t_1))

        sol_fin = sol_t_1 .+ sol_t_2

    end 

    if model == "dHPM_inhibition"

        sol_t_2 = sol_fin[2,:]
        sol_t_2 =  reshape(sol_t_2, 1,length(sol_t_2))

        sol_t_1 = sol_fin[1,:]
        sol_t_1 =  reshape(sol_t_1, 1,length(sol_t_1))

        sol_fin = sol_t_1 .+ sol_t_2

    end 


  






    if type_of_loss == "L2"

        if size(sol_fin)[2] == size(data)[2]
            loss = NaNMath.sum(abs2.((data[2, :] - sol_fin[1, 1:end]))) / length(data[2, :])
        else
            loss = NaNMath.sum(abs2.((data[2, 1:size(sol_fin)[2]] - sol_fin[1, 1:end]))) / length(data[2, 1:size(sol_fin)[2]]) + 1000 * (size(data)[2] - size(sol_fin)[2])
        end
    end


    if type_of_loss == "RE"

        if size(sol_fin)[2] == size(data)[2]
            loss = 0.5 * NaNMath.sum(abs2.(1.0 .- (data[2, :] ./ sol_fin[1, 1:end]))) / length(data[2, :])
        else
            loss = 0.5 * NaNMath.sum(abs2.(1.0 .- (data[2, 1:size(sol_fin)[2]] ./ sol_fin[1, 1:end]))) / length(data[2, 1:size(sol_fin)[2]]) + 1000 * (size(data)[2] - size(sol_fin)[2])
        end
    end

    if type_of_loss == "L2_derivative"

        derivative_theo = [(sol_fin[1, i] - sol_fin[1, i-1]) / (data[1, i] - data[1, i-1]) for i in 2:(size(sol_fin)[2])]

        derivative_data = [(data[2, i] - data[2, i-1]) / (data[1, i] - data[1, i-1]) for i in 2:(size(data)[2])]

        if size(derivative_theo)[1] == size(derivative_data)[1]
            loss = NaNMath.sum(abs2.((derivative_theo[2, :] .- derivative_data[1, 1:end]))) / length(derivative_data[2, :])
        else
            loss = NaNMath.sum(abs2.((derivative_theo[2, :] .- derivative_data[1, 1:end]))) / length(derivative_data[2, :]) + 1000
        end

    end

    if type_of_loss == "blank_weighted_L2"


        empirical_blank_distrib = (blank_array .- mean(blank_array))
        # generation of the empirica distrib respect the mean of the noise 
        test = StatsBase.fit(Histogram, empirical_blank_distrib)
        binning_distrib = test.edges
        probabiltity_distrib = test.weights ./ sum(test.weights)

        if size(sol_fin)[2] == size(data)[2]
            loss = 0.0



            for ll in 1:size(sol_fin)[2]
                dist = (data[2, ll] - sol_fin[1, ll])
                index = findfirst(x -> (x > dist), binning_distrib[1])
                if (typeof(index) == Nothing || index > length(probabiltity_distrib))

                    loss = loss + (abs2.(dist)) / length(data[2, :])
                else
                    prob = probabiltity_distrib[index]

                    loss = loss + (abs2.((1 - prob) * dist)) / length(data[2, :])

                end

            end
        else

            loss = 10.0^9 * length(data[2, :])

        end


    end


    # parameters have different form in function of the used model
   
    #max_theoretical gr
    max_th_gr = 0.0
    # max empirical gr
    max_em_gr = 0.0 

    res_temp = res.u

    res_param = vectorize_df_results(name_well,
        model,
        res_temp,
        max_th_gr,
        max_em_gr,
        loss
    )



    



    return res_param

end







function  selection_ODE_fixed_change_points(data_testing::Matrix{Float64}, # dataset x times y OD/fluorescence
    name_well::String, # name of the well
    label_exp::String, #label of the experiment
    list_lb_param::Any, # lower bound param
    list_ub_param::Any, # upper bound param
    param_guess::Any,# list of param
    list_of_models::Vector{String}, # ode model to use 
    do_plot::Bool, # do plots or no
    path_to_plot::String, # where save plots
    n_change_points::Int,
    win_size::Int, # numebr of the point to generate intial condition
    pt_smooth_derivative::Int,
    correction::Bool, # if true uses the given calibration curve to fix the data
    calibration_OD_curve ::String  #  the path to calibration curve to fix the data
    )
    # inizialization penality  function

  







    if  correction == true

        data_testing = correction_OD_multiple_scattering(data_testing,calibration_OD_curve)
    

    end

    list_change_points_dev = detect_list_change_points_derivative( data_testing,n_change_points,win_size,pt_smooth_derivative)
    interval_changepoints = push!( list_change_points_dev[2],data_testing[1,1])
    interval_changepoints = push!( list_change_points_dev[2],data_testing[1,end])
    interval_changepoints = sort(interval_changepoints)
    bc = [data_testing[1,1],data_testing[2,2]]

    temp_res_model_selection= Type{Any}
    param_out =  Vector{Vector{Any}}()
    composed_sol =Type{Any}
    composed_time = Type{Any}

  for i in 2:(length(interval_changepoints))
        	aic_old= 10^9

        if i==2
    
            tspan_array = findall(( data_testing[1,:] .<= interval_changepoints[i] ))
    
            data_temp =  Matrix(transpose(hcat(data_testing[1,tspan_array], data_testing[2,tspan_array])));
    
        else
            
            tspan_array_1 = findall((data_testing[1,:] .> interval_changepoints[i-1]  ))

        
            tspan_array_2 = findall(( data_testing[1,:] .<= interval_changepoints[i] ))
            tspan_array = intersect(tspan_array_1,tspan_array_2)
    
                data_temp =  Matrix(transpose(hcat(data_testing[1,tspan_array], data_testing[2,tspan_array])));
               # imposing_bounduary condition
              data_temp= hcat( bc,data_temp)
        end
    

        for mm in 1:length(list_of_models)
    


            model =  list_of_models[mm]
            temp_temp_res_win = fitting_one_ode_iterate(data_temp, # dataset x times y OD/fluorescence
            name_well, # name of the well
            list_lb_param[mm], # lower bound param
            list_ub_param[mm], # upper bound param
            param_guess[mm], # initial guess param
            model, # ode model to use 
            1, # numebr of the point to generate intial condition
            false, # the smoothing is done or not?
            "L2", # type of used loss 
            zeros(9), # data of all blanks
            false, # if true uses the given calibration curve to fix the data
            "no" #  the path to calibration curve to fix the data
            )
    
            param_fitting = temp_temp_res_win[3:(end-3)]

            aic =  length(data_temp[2,:]) * log(  temp_temp_res_win[end] / length(data_temp[2,:]) ) + 2* length(param_fitting)

            if aic < aic_old


                aic_old =copy(aic)
                temp_temp_res_win = push!(temp_temp_res_win,aic)
                temp_res_model_selection = copy(temp_temp_res_win)

            end


            
        end
           # selection of te model
           model =  temp_res_model_selection[1]

          # param of the best model
           param_fitting = temp_res_model_selection[3:(end-4)]

           temp_res_win = temp_res_model_selection
            u0 = [data_temp[2,1]]
    
        if model == "HPM_3_death" || model == "HPM_3_inhibition"|| model == "HPM_3_death_resistance"|| model == "dHPM_3_death_resistance"
            u0 = [data_temp[2,1],0.0,0.0]
    
        end 
        if model == "HPM" || model == "HPM_exp" || model == "dHPM" || model == "HPM_inhibition" || model == "dHPM_inhibition"  || model=="ODEs_HPM_SR"
            u0 = [data_temp[2,1],0.0]
    
    
        end    

            remade_solution = ODE_sim_for_iterate(model, #string of the model
            u0, # starting condition
            data_temp[1,:], # start time of the sim
            KenCarp4(), # which sciml solver of ode
           param_fitting # parameters of the ODE model
        )


          time_sol = reduce(hcat, remade_solution.t)
            sol_fin = reduce(hcat, remade_solution.u)
    
    
        if model == "HPM_3_death" || model == "HPM_3_inhibition"|| model == "HPM_3_death_resistance"|| model == "dHPM_3_death_resistance"
    
            sol_t_2 = sol_fin[2,:]
            sol_t_2 =  reshape(sol_t_2, 1,length(sol_t_2))
    
            sol_t_1 = sol_fin[1,:]
            sol_t_1 =  reshape(sol_t_1, 1,length(sol_t_1))
    
    
            sol_t_3 = sol_fin[3,:]
            sol_t_3 =  reshape(sol_t_3, 1,length(sol_t_3))
    
            sol_fin = sol_t_3 .+ sol_t_2.+ sol_t_1
    
    
    
        end 
    
    
        if (model == "HPM" || model == "HPM_exp")
    
            sol_t_2 = sol_fin[2,:]
            sol_t_2 =  reshape(sol_t_2, 1,length(sol_t_2))
    
            sol_t_1 = sol_fin[1,:]
            sol_t_1 =  reshape(sol_t_1, 1,length(sol_t_1))
    
            sol_fin = sol_t_1 .+ sol_t_2
    
    
    
        end 
        if model == "ODEs_HPM_SR"
    
            sol_t_2 = sol_fin[2,:]
            sol_t_2 =  reshape(sol_t_2, 1,length(sol_t_2))
    
            sol_t_1 = sol_fin[1,:]
            sol_t_1 =  reshape(sol_t_1, 1,length(sol_t_1))
    
            sol_fin = sol_t_1 .+ sol_t_2
        
    
    
    
        end 
    
        if model == "dHPM"
    
            sol_t_2 = sol_fin[2,:]
            sol_t_2 =  reshape(sol_t_2, 1,length(sol_t_2))
    
            sol_t_1 = sol_fin[1,:]
            sol_t_1 =  reshape(sol_t_1, 1,length(sol_t_1))
    
            sol_fin = sol_t_1 .+ sol_t_2
    
        end 
    
    
        if model == "HPM_inhibition"
    
            sol_t_2 = sol_fin[2,:]
            sol_t_2 =  reshape(sol_t_2, 1,length(sol_t_2))
    
            sol_t_1 = sol_fin[1,:]
            sol_t_1 =  reshape(sol_t_1, 1,length(sol_t_1))
    
            sol_fin = sol_t_1 .+ sol_t_2
    
        end 
    
        if model == "dHPM_inhibition"
    
            sol_t_2 = sol_fin[2,:]
            sol_t_2 =  reshape(sol_t_2, 1,length(sol_t_2))
    
            sol_t_1 = sol_fin[1,:]
            sol_t_1 =  reshape(sol_t_1, 1,length(sol_t_1))
    
            sol_fin = sol_t_1 .+ sol_t_2
    
        end 
    
        value_bonduary = remade_solution.t[end]
        time_bonduary= sol_fin[end]
        bc = [value_bonduary,time_bonduary]
        if (pt_smooth_derivative > (length(data_temp[1,:]) - 2) && length(data_temp[1,:]) >3 )

          emp_max_gr_of_segment =  maximum(specific_gr_evaluation(data_temp, 0))
          temp_res_win[(end-2)]= emp_max_gr_of_segment
  
           data_th =Matrix(vcat(time_sol, sol_fin)) 
          th_max_gr_of_segment = maximum(specific_gr_evaluation(data_th,0))
          temp_res_win[(end-3)]= th_max_gr_of_segment

        elseif length(data_temp[1,:])<=3

            emp_max_gr_of_segment =  "NA"
            temp_res_win[(end-2)]= emp_max_gr_of_segment
            th_max_gr_of_segment ="NA"
            temp_res_win[(end-3)]= th_max_gr_of_segment

      else
          emp_max_gr_of_segment =  maximum(specific_gr_evaluation(data_temp, pt_smooth_derivative))
          temp_res_win[(end-2)]= emp_max_gr_of_segment
  
           data_th =Matrix(vcat(time_sol, sol_fin)) 
          th_max_gr_of_segment = maximum(specific_gr_evaluation(data_th,pt_smooth_derivative))
          temp_res_win[(end-3)]= th_max_gr_of_segment
      end
     

     

        if i == 2
            composed_time = copy(time_sol)
            composed_sol = copy(sol_fin)

            temp_res_win=  push!(temp_res_win,i-1)
            param_out = push!(param_out,temp_res_win)
    
        else
            composed_time = hcat(composed_time,reduce(hcat,time_sol))

            composed_sol =hcat(composed_sol,reduce(hcat,sol_fin))
            temp_res_win=  push!(temp_res_win,i-1)
            param_out =  push!(param_out,temp_res_win)
        
        end
    
    end
   

    if do_plot == true

        mkpath(path_to_plot)
        
        display(Plots.scatter(data_testing[1, :], data_testing[2, :], xlabel="Time", ylabel="Arb. Units", label=["Data " nothing], markersize=1, color=:black, title=string(label_exp, " ", name_well)))
        display( Plots.vline!(interval_changepoints[2:end], c=:black, label=["change points derivative" nothing]))
        display(Plots.plot!(reduce(vcat,composed_time), reduce(vcat,composed_sol), xlabel="Time", ylabel="Arb. Units", label=[" fitting " nothing], color=:red, title=string(label_exp, " fitting ", name_well)))
        png(string(path_to_plot, label_exp, "_model_selecetion_seg_",n_change_points,"_", name_well, ".png"))

    end

    return param_out,interval_changepoints,composed_time,composed_sol
end


function   ODE_selection_NMAX_change_points(data_testing::Matrix{Float64}, # dataset x times y OD/fluorescence
    name_well::String, # name of the well
    label_exp::String, #label of the experiment
    list_lb_param::Any, # lower bound param
    list_ub_param::Any, # upper bound param
    param_guess::Any,# list of param
    list_of_models::Vector{String}, # ode model to use 
    do_plot::Bool, # do plots or no
    path_to_plot::String, # where save plots
    path_to_results::String,
    n_change_points::Int,
    win_size::Int, # numebr of the point to generate intial condition
    pt_smooth_derivative::Int,
    penality_parameter::Float64,
    correction::Bool, # if true uses the given calibration curve to fix the data
    calibration_OD_curve ::String,  #  the path to calibration curve to fix the data
   save_all_model::Bool )

    # fitting single models

    change_point_list =  Vector{Vector{Any}}()

    res = ODE_Model_selection(data_testing, # dataset x times y OD/fluorescence
         name_well, # name of the well
         label_exp, #label of the experiment
         list_lb_param, # lower bound param
         list_ub_param, # upper bound param
         list_of_models, # ode model to use 
         1, # numebr of the point to generate intial condition
        false, # the smoothing is done or not?
         "L2", # type of used loss 
         zeros(9), # data of all blanks
         false, # one wants the results of the best fit to be plotted
         path_to_plot,
         false,
         correction, # if true uses the given calibration curve to fix the data
         calibration_OD_curve   #  the path to calibration curve to fix the data
     
         )

         if save_all_model == true
            mkpath(path_to_results)
            CSV.write(string(path_to_results,label_exp ,"_segmented_ODE_scoring_", name_well,"_seg_0.csv"), Tables.table(Matrix(res[1])))

            CSV.write(string(path_to_results,label_exp ,"_segmented_ODE_param_", name_well,"_seg_0.csv"), Tables.table(Vector(res[2])))

        end

            top_model = vcat(res[6],res[5])


            loss_penality = res[4]
            println(loss_penality)
            n_param = length(res[5] )
            length_dataset =length(data_testing[2,:])

            score_of_the_models  =  length_dataset * log(  loss_penality[1] /length_dataset) + penality_parameter * n_param 
            println(score_of_the_models)

         elbow_plot = beta_scanning(length_dataset,
         loss_penality[1],
         0.1,
         5.0,
         0.1,
         n_param)  

         change_point_list = [0.0]
         change_point_to_plot = [0.0,0.0,0.0]
         time_points_to_plot = [0.0,0.0,0.0]
         sol_to_plot = [0.0,0.0,0.0]
    #fitting all model with change points
    if n_change_points > 0
        for n in 1:n_change_points
                 println("testing n intevarls ", n+1)
            direct_search_results =   selection_ODE_fixed_change_points(data_testing, # dataset x times y OD/fluorescence
                name_well, # name of the well
                label_exp, #label of the experiment
                list_lb_param, # lower bound param
                list_ub_param, # upper bound param
                param_guess,# list of param
                list_of_models, # ode model to use 
                false, # do plots or no
                path_to_plot, # where save plots
                n,
                win_size, # numebr of the point to generate intial condition
                pt_smooth_derivative,
                correction, # if true uses the given calibration curve to fix the data
                calibration_OD_curve  #  the path to calibration curve to fix the data
            )
            # composing piecewise penality
            loss_penality = sum([direct_search_results[1][kk][(end-2)] for kk in 1:length(direct_search_results[1])])
            
            n_param = sum([ length(direct_search_results[1][kk][3:(end-5)] ) for kk in 1:length(direct_search_results[1])]) + n_change_points
            
            
            new_penality  = length_dataset * log(  loss_penality /length_dataset) + penality_parameter * n_param 
       #     println(direct_search_results[1])
       new_elbow =   beta_scanning(length_dataset,
         loss_penality[1],
         0.1,
         5.0,
         0.1,
         n_param) 
            println(new_penality)

            println("######")
            elbow_plot = hcat(elbow_plot,new_elbow)       


            if new_penality <= score_of_the_models
                score_of_the_models = copy(new_penality)
               top_model = copy(direct_search_results[1])
               time_points_to_plot = copy(direct_search_results[3])
               sol_to_plot  = copy(direct_search_results[4])
               change_point_to_plot = copy(direct_search_results[2])
            end
           
           
            if save_all_model == true
                CSV.write(string(path_to_results,label_exp ,"_segmented_ODE_", name_well,"_seg_",n,".csv"), Tables.table(Vector(direct_search_results[1])))

            end

        end
    end

    # plotting best model if required 
    if do_plot == true

        mkpath(path_to_plot)
        
        display(Plots.scatter(data_testing[1, :], data_testing[2, :], xlabel="Time", ylabel="Arb. Units", label=["Data " nothing], markersize=1, color=:black, title=string(label_exp, " ", name_well)))
        display( Plots.vline!(change_point_to_plot[2:end], c=:black, label=["change points derivative" nothing]))
        display(Plots.plot!(reduce(vcat,time_points_to_plot), reduce(vcat,sol_to_plot), xlabel="Time", ylabel="Arb. Units", label=[" fitting " nothing], color=:red, title=string(label_exp, " fitting ", name_well)))
        png(string(path_to_plot, label_exp, "_model_selecetion_seg_",length(change_point_to_plot[2:end]),"_", name_well, ".png"))

    end

    return top_model,elbow_plot,time_points_to_plot,sol_to_plot


end





function beta_scanning(length_dataset,
    loss_penality,
    penality_parameter_min,
    penality_parameter_max,
    delta_pen,
    n_param)  
    
    
    
    score_of_the_models  =  [length_dataset * log(  loss_penality/length_dataset) + beta * n_param for beta in penality_parameter_min:delta_pen:penality_parameter_max ]

    return score_of_the_models
end

