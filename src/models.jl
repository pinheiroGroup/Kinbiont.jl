########################## No CONSTRAIN
# heterogeneous population model (HPM) with dormant active and starving population

struct JModel
  name::String
  func::Function
  params::Vector{String}
end

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

function aHPM_3_death_resistance(du, u, param, t)
  du[1] = - u[1] * param[2]
  du[2] = u[1] * param[2] + param[1] * u[2]  -  param[3]*  u[2]#*(1-u[2])
  du[3] = param[3]*  u[2]  + param[4]*  u[3]* (1 - NaNMath.pow((u[3]+u[1]+u[2])/param[5],param[6]))
end

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

function ODEs_aHPM_inhibition(du, u, param, t)
  du[1] =  - u[1] * (param[2])  +  param[1] *  u[1] 
  du[2] = + u[1] * (param[2]) + param[4] * u[2]* ( 1 -  NaNMath.pow( (u[1] + u[2])/param[3],param[5])) 
end

function ODEs_HPM_SR(du, u, param, t)
  du[1] =  -  param[4]  / (1+ param[3] * exp(- param[2] * t) )*  u[1] +  param[1] *  u[1] -   param[5]  * u[1] 
  du[2] =  +  param[5]  * u[1] 
end

# adjusted heterogeneous population model (aHPM).
function ODEs_adjusted_McKellar(du, u, param, t)
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
function ODE_Diauxic_piecewise_adjusted_logistic(du, u, param, t)
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

function ODE_Diauxic_piecewise_adjusted_logistic_bk(du, u, param, t)
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

function ODE_piecewise_adjusted_logistic(du, u, param, t)
  if t <= param[3]
    du[1] = param[5]
  else
    du[1] = param[1] * u[1] * (1 - NaNMath.pow((u[1] / param[2]), param[4]))
  end
end

function ODE_triple_piecewise_adjusted_logistic(du, u, param, t)
  if t <= param[3]
    du[1] = u[1] * param[5]
  elseif ( t <=  param[6] && t > param[3] )
    du[1] = param[1] * u[1] * (1 - NaNMath.pow((u[1] / param[2]), param[4]))
  elseif ( t >  param[6] && t > param[3] ) 
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

# aLogistic curve 
function alogistic(du, u, param, t)
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

models_list = [
  JModel(
    "HPM_3_inhibition",
    HPM_3_inhibition,
    ["model", "well", "gr", "exit_lag_rate", "inactivation_rate", "th_max_gr","emp_max_gr","loss"]
  ),
  JModel(
    "HPM_3_death",
    HPM_3_death,
    ["model", "well", "gr", "exit_lag_rate", "inactivation_rate","death_rate","th_max_gr","emp_max_gr","loss"]
  ),
  JModel(
    "aHPM_3_death_resistance",
    aHPM_3_death_resistance,
    ["model", "well", "gr", "exit_lag_rate", "inactivation_rate","death_rate","n_res","shape","th_max_gr","emp_max_gr","loss"]
   ),
  JModel(
    "HPM",
    ODEs_McKellar,
    ["model", "well", "gr", "exit_lag_rate", "N_max","th_max_gr","emp_max_gr","loss"]
   ),
  JModel(
    "HPM_exp",
    ODEs_HPM_exp,
    ["model", "well", "gr", "exit_lag_rate","th_max_gr","emp_max_gr","loss"]
   ),
  JModel(
    "HPM_inhibition",
    ODEs_HPM_inhibition,
    ["model", "well", "gr", "inhibition_rate", "gr_inhibition","N_max","th_max_gr","emp_max_gr", "loss"]
   ),
  JModel(
    "aHPM_inhibition",
    ODEs_aHPM_inhibition,
    ["model", "well", "gr", "inhibition_rate", "gr_inhibition","N_max", "shape","th_max_gr","emp_max_gr" , "loss"]
   ),
  JModel(
    "ODEs_HPM_SR",
    ODEs_HPM_SR,
    ["model", "well", "gr", "gr_phage","scale", "death_rate","resistance_rate","th_max_gr","emp_max_gr", "loss"]
   ),
  JModel(
    "aHPM",
    ODEs_adjusted_McKellar,
    ["model", "well", "gr", "exit_lag_rate", "N_max", "shape" ,"th_max_gr","emp_max_gr", "loss"]
   ),
  JModel(
    "Diauxic_replicator_1",
    ODE_Diauxic_replicator_1,
    ["model", "well", "gr", "N_max", "lag", "arbitrary_const", "linear_const","th_max_gr","emp_max_gr", "loss"]
   ),
  JModel(
    "Diauxic_replicator_2",
    ODE_Diauxic_replicator_2,
    ["model", "well", "gr", "N_max", "lag", "arbitrary_const", "linear_const", "growth_stationary","th_max_gr","emp_max_gr", "loss"]
   ),
  JModel(
    "Diauxic_piecewise_adjusted_logistic",
    ODE_Diauxic_piecewise_adjusted_logistic,
    ["model", "well", "gr_1", "N_max", "shape_1", "lag", "linear_const", "t_shift", "gr_2","N_max_2","shape_2","end_second_lag","lag_2_gr","th_max_gr","emp_max_gr", "loss"]
   ),
  JModel(
    "exponential",
    ODE_exponential,
    ["model", "well", "gr", "th_max_gr","emp_max_gr", "loss"]
   ),
  JModel(
    "piecewise_adjusted_logistic",
    ODE_piecewise_adjusted_logistic,
    ["model", "well", "gr", "N_max", "lag", "shape", "linear_const","th_max_gr","emp_max_gr", "loss"]
   ),
  JModel(
    "triple_piecewise_adjusted_logistic",
    ODE_triple_piecewise_adjusted_logistic,
    ["model", "well", "gr", "N_max", "lag", "shape", "linear_const", "t_stationary", "linear_lag","th_max_gr","emp_max_gr", "loss"]
   ),
  JModel(
    "triple_piecewise",
    ODE_triple_piecewise,
    ["model", "well", "gr", "gr_2", "gr_3", "lag", "t_stationary","th_max_gr","emp_max_gr", "loss"]
   ),
  JModel(
    "triple_piecewise_sub_linear",
    ODE_triple_piecewise_sub_linear,
    ["model", "well", "gr", "gr_2", "gr_3", "lag", "t_stationary","N_max","th_max_gr","emp_max_gr", "loss"]
   ),
  JModel(
    "gbsm_piecewise",
    ODE_gbsm_piecewise,
    ["model", "well", "gr", "a_1", "b_1", "c", "a_2", "b_2","th_max_gr","emp_max_gr", "loss"]
   ),
  JModel(
    "ODE_four_piecewise",
    ODE_four_piecewise,
    ["model", "well", "gr", "gr_2", "gr_3", "gr_4", "lag", "t_decay_gr", "t_stationary","th_max_gr","emp_max_gr", "loss"]
   ),
  JModel(
    "hyper_gompertz",
    hyper_gompertz,
    ["model", "well", "gr", "N_max", "shape","th_max_gr","emp_max_gr", "loss"]
   ),
  JModel(
    "hyper_logistic",
    hyper_logistic,
    ["model", "well", "doubling_time", "gr", "N_max", "shape","th_max_gr","emp_max_gr", "loss"]
   ),
  JModel(
    "bertalanffy_richards",
    bertalanffy_richards,
    ["model", "well", "gr", "N_max", "shape","th_max_gr","emp_max_gr", "loss"]
   ),
  JModel(
    "ode_von_bertalanffy",
    ODE_von_bertalanffy,
    ["model", "well", "alpha", "beta", "a","b","th_max_gr","emp_max_gr", "loss"]
   ),
  JModel(
    "triple_piecewise_bertalanffy_richards",
    ODE_triple_piecewise_bertalanffy_richards,
    ["model", "well", "gr", "gr_lag","t_lag","t_stationary","gr_stat", "shape", "N_max","th_max_gr","emp_max_gr","loss"]
   ),
  JModel(
    "logistic",
    logistic,
    ["model", "well", "gr", "N_max","th_max_gr","emp_max_gr", "loss"]
   ),
  JModel(
    "alogistic",
    alogistic,
    ["model", "well", "gr", "N_max","shape","th_max_gr","emp_max_gr", "loss"]
   ),
  JModel(
    "gompertz",
    gompertz,
    ["model", "well", "gr", "N_max","th_max_gr","emp_max_gr", "loss"]
   ),
  JModel(
    "baranyi_richards",
    baranyi_richards,
    ["model", "well", "gr", "N_max", "lag_time", "shape","th_max_gr","emp_max_gr", "loss"]
   ),
  JModel(
    "baranyi_exp",
    baranyi_exp,
    ["model", "well", "gr", "lag_time", "shape","th_max_gr","emp_max_gr", "loss"]
   ),
  JModel(
    "baranyi_roberts",
    baranyi_roberts,
    ["model", "well", "gr", "N_max", "lag_time", "shape_1", "shape_2","th_max_gr","emp_max_gr", "loss"]
   )
]

models = Dict(model.name => model for model in models_list)
