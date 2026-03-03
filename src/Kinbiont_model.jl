using NaNMath

struct new_Kinbiont_generic_model
    name::String
    func::Function
    guess::Union{Function,Nothing}
    params::Vector{String}
    Type_of_model::String

end
Kinbiont_models_list = [
  new_Kinbiont_generic_model(
        "NL_piecewise_lin_logistic",
        NL_piecewise_linear_logistic,
        guess_NL_model_exp,
        ["model", "well", "t_lag", "N_lag", "N_max","growth_rate","th_max_gr", "emp_max_gr", "loss"],
        "NL"
    ),
    new_Kinbiont_generic_model(
        "NL_piecewise_exp_logistic",
        NL_piecewise_exp_logistic,
        guess_NL_model_exp,
        ["model", "well", "t_lag", "N_lag", "N_max","growth_rate","linear_rate","th_max_gr", "emp_max_gr", "loss"],
        "NL"
    ),
    new_Kinbiont_generic_model(
        "NL_exponential",
        NL_model_exp,
        guess_NL_model_exp,
        ["model", "well", "N0", "growth_rate", "th_max_gr", "emp_max_gr", "loss"],
        "NL"
    ),
    new_Kinbiont_generic_model(
        "NL_logistic",
        NL_model_logistic,
        guess_NL_model_logistic,
        ["model", "well", "N_max", "growth_rate", "lag", "th_max_gr", "emp_max_gr", "loss"],
        "NL"
    ),
    new_Kinbiont_generic_model(
        "NL_Gompertz",
        NL_model_Gompertz,
        guess_NL_model_Gompertz,
        ["model", "well", "N_max", "growth_rate", "lag", "th_max_gr", "emp_max_gr", "loss"],
        "NL"
    ),
    new_Kinbiont_generic_model(
        "NL_Bertalanffy",
        NL_model_Bertalanffy,
        guess_NL_model_Bertalanffy,
        ["model", "well", "N_0", "N_max", "growth_rate", "shape", "th_max_gr", "emp_max_gr", "loss"],
        "NL"
    ),
    new_Kinbiont_generic_model(
        "NL_Richards",
        NL_model_Richards,
        guess_NL_model_Richards,
        ["model", "well", "N_max", "shape", "growth_rate", "lag", "th_max_gr", "emp_max_gr", "loss"],
        "NL"
    ),
    new_Kinbiont_generic_model(
        "NL_Morgan",
        NL_model_Morgan,
        guess_NL_model_Morgan,
        ["model", "well", "N_0", "K", "shape", "N_max", "th_max_gr", "emp_max_gr", "loss"],
        "NL"
    ),
    new_Kinbiont_generic_model(
        "NL_Weibull",
        NL_model_Weibull,
        guess_NL_model_Weibull,
        ["model", "well", "N_max", "N_0", "growth_rate", "shape", "th_max_gr", "emp_max_gr", "loss"],
        "NL"
    ),
    new_Kinbiont_generic_model(
    "HPM_3_inhibition",
    HPM_3_inhibition,
    nothing,
    ["label_exp", "well", "model", "gr", "exit_lag_rate", "inactivation_rate", "th_max_gr", "emp_max_gr", "loss"],
    "ODE"
  ),
  new_Kinbiont_generic_model(
    "HPM_3_death",
    HPM_3_death,
    nothing,
    ["label_exp", "well", "model", "gr", "exit_lag_rate", "inactivation_rate", "death_rate", "th_max_gr", "emp_max_gr", "loss"],
    "ODE"
  ),
  new_Kinbiont_generic_model(
    "aHPM_3_death_resistance",
    aHPM_3_death_resistance,
    nothing,
    ["label_exp", "well", "model", "gr", "exit_lag_rate", "inactivation_rate", "death_rate", "n_res", "shape", "th_max_gr", "emp_max_gr", "loss"],
    "ODE"
  ),
  new_Kinbiont_generic_model(
    "HPM",
    ODEs_McKellar,
    nothing,
    ["label_exp", "well", "model", "gr", "exit_lag_rate", "N_max", "th_max_gr", "emp_max_gr", "loss"],
    "ODE"
  ),
  new_Kinbiont_generic_model(
    "HPM_exp",
    ODEs_HPM_exp,
    nothing,
    ["label_exp", "well", "model", "gr", "exit_lag_rate", "th_max_gr", "emp_max_gr", "loss"],
    "ODE"
  ),
  new_Kinbiont_generic_model(
    "HPM_inhibition",
    ODEs_HPM_inhibition,
    nothing,
    ["label_exp", "well", "model", "gr", "inhibition_rate", "gr_inhibition", "N_max", "th_max_gr", "emp_max_gr", "loss"],
    "ODE"
  ),
  new_Kinbiont_generic_model(
    "aHPM_inhibition",
    ODEs_aHPM_inhibition,
    nothing,
    ["label_exp", "well", "model", "gr", "inhibition_rate", "gr_inhibition", "N_max", "shape", "th_max_gr", "emp_max_gr", "loss"],
    "ODE"
  ),
  new_Kinbiont_generic_model(
    "ODEs_HPM_SR",
    ODEs_HPM_SR,    
    nothing,
    ["label_exp", "well", "model", "gr", "gr_phage", "scale", "death_rate", "resistance_rate", "th_max_gr", "emp_max_gr", "loss"],
    "ODE"
  ),
  new_Kinbiont_generic_model(
    "aHPM",
    ODEs_adjusted_McKellar,
    nothing,
    ["label_exp", "well", "model", "gr", "exit_lag_rate", "N_max", "shape", "th_max_gr", "emp_max_gr", "loss"],
    "ODE"
  ),
  new_Kinbiont_generic_model(
    "Diauxic_replicator_1",
    ODE_Diauxic_replicator_1,
    nothing,
    ["label_exp", "well", "model", "gr", "N_max", "lag", "arbitrary_const", "linear_const", "th_max_gr", "emp_max_gr", "loss"],
    "ODE"
  ),
  new_Kinbiont_generic_model(
    "Diauxic_replicator_2",
    ODE_Diauxic_replicator_2,
    nothing,
    ["label_exp", "well", "model", "gr", "N_max", "lag", "arbitrary_const", "linear_const", "growth_stationary", "th_max_gr", "emp_max_gr", "loss"],
    "ODE"
  ),
  new_Kinbiont_generic_model(
    "Diauxic_piecewise_adjusted_logistic",
    ODE_Diauxic_piecewise_adjusted_logistic,
    nothing,
    ["label_exp", "well", "model", "gr_1", "N_max", "shape_1", "lag", "linear_const", "t_shift", "gr_2", "N_max_2", "shape_2", "end_second_lag", "lag_2_gr", "th_max_gr", "emp_max_gr", "loss"],
    "ODE"
  ),
  new_Kinbiont_generic_model(
    "exponential",
    ODE_exponential,
    nothing,
    ["label_exp", "well", "model", "gr", "th_max_gr", "emp_max_gr", "loss"],
    "ODE"
  ),
  new_Kinbiont_generic_model(
    "piecewise_adjusted_logistic",
    ODE_piecewise_adjusted_logistic,
    nothing,
    ["label_exp", "well", "model", "gr", "N_max", "lag", "shape", "linear_const", "th_max_gr", "emp_max_gr", "loss"],
    "ODE"
  ),
  new_Kinbiont_generic_model(
    "triple_piecewise_adjusted_logistic",
    ODE_triple_piecewise_adjusted_logistic, 
       nothing,
    ["label_exp", "well", "model", "gr", "N_max", "lag", "shape", "linear_const", "t_stationary", "linear_lag", "th_max_gr", "emp_max_gr", "loss"],
    "ODE"
  ),
  new_Kinbiont_generic_model(
    "triple_piecewise",
    ODE_triple_piecewise,
    nothing,
    ["label_exp", "well", "model", "gr", "gr_2", "gr_3", "lag", "t_stationary", "th_max_gr", "emp_max_gr", "loss"],
    "ODE"
  ),
  new_Kinbiont_generic_model(
    "triple_piecewise_sub_linear",
    ODE_triple_piecewise_sub_linear,
    nothing,
    ["label_exp", "well", "model", "gr", "gr_2", "gr_3", "lag", "t_stationary", "N_max", "th_max_gr", "emp_max_gr", "loss"],
    "ODE"
  ),
  new_Kinbiont_generic_model(
    "gbsm_piecewise",
    ODE_gbsm_piecewise,
    nothing,
    ["label_exp", "well", "model", "gr", "a_1", "b_1", "c", "a_2", "b_2", "th_max_gr", "emp_max_gr", "loss"],
    "ODE"
  ),
  new_Kinbiont_generic_model(
    "ODE_four_piecewise",
    ODE_four_piecewise,
    nothing,
    ["label_exp", "well", "model", "gr", "gr_2", "gr_3", "gr_4", "lag", "t_decay_gr", "t_stationary", "th_max_gr", "emp_max_gr", "loss"],
    "ODE"
  ),
  new_Kinbiont_generic_model(
    "hyper_gompertz",
    hyper_gompertz,
    nothing,
    ["label_exp", "well", "model", "gr", "N_max", "shape", "th_max_gr", "emp_max_gr", "loss"],
    "ODE"
  ),
  new_Kinbiont_generic_model(
    "hyper_logistic",
    hyper_logistic,
    nothing,
    ["label_exp", "well", "model", "doubling_time", "gr", "N_max", "shape", "th_max_gr", "emp_max_gr", "loss"],
    "ODE"
  ),
  new_Kinbiont_generic_model(
    "bertalanffy_richards",
    bertalanffy_richards,
    nothing,
    ["label_exp", "well", "model", "gr", "N_max", "shape", "th_max_gr", "emp_max_gr", "loss"],
    "ODE"
  ),
  new_Kinbiont_generic_model(
    "ode_von_bertalanffy",
    ODE_von_bertalanffy,
    nothing,
    ["label_exp", "well", "model", "alpha", "beta", "a", "b", "th_max_gr", "emp_max_gr", "loss"],
    "ODE"
  ),
  new_Kinbiont_generic_model(
    "triple_piecewise_bertalanffy_richards",
    ODE_triple_piecewise_bertalanffy_richards,
    nothing,
    ["label_exp", "well", "model", "gr", "gr_lag", "t_lag", "t_stationary", "gr_stat", "shape", "N_max", "th_max_gr", "emp_max_gr", "loss"],
    "ODE"
  ),
  new_Kinbiont_generic_model(
    "logistic",
    logistic,
    nothing,
    ["label_exp", "well", "model", "gr", "N_max", "th_max_gr", "emp_max_gr", "loss"],
    "ODE"
  ),
  new_Kinbiont_generic_model(
    "alogistic",
    alogistic,
    nothing,
    ["label_exp", "well", "model", "gr", "N_max", "shape", "th_max_gr", "emp_max_gr", "loss"],
    "ODE"
  ),
  new_Kinbiont_generic_model(
    "gompertz",
    gompertz,
    nothing,
    ["label_exp", "well", "model", "gr", "N_max", "th_max_gr", "emp_max_gr", "loss"],
    "ODE"
  ),
  new_Kinbiont_generic_model(
    "baranyi_richards",
    baranyi_richards,
    nothing,
    ["label_exp", "well", "model", "gr", "N_max", "lag_time", "shape", "th_max_gr", "emp_max_gr", "loss"],
    "ODE"
  ),
  new_Kinbiont_generic_model(
    "baranyi_exp",
    baranyi_exp,
    nothing,
    ["label_exp", "well", "model", "gr", "lag_time", "shape", "th_max_gr", "emp_max_gr", "loss"],
    "ODE"
  ),
  new_Kinbiont_generic_model(
    "baranyi_roberts",
    baranyi_roberts,
        nothing,
    ["label_exp", "well", "model", "gr", "N_max", "lag_time", "shape_1", "shape_2", "th_max_gr", "emp_max_gr", "loss"],
    "ODE"
  ),
]


Kinbiont_models = Dict(new_Kinbiont_generic_model.name => new_Kinbiont_generic_model for new_Kinbiont_generic_model in Kinbiont_models_list)
