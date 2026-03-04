struct Kinbiont_res_one_well_log_lin
    method::String
    params::Vector{Any}
    fit::Any
    times::Any
    confidence_band::Any
end

struct Kinbiont_res_Cybernetic_model
    method::String
    params::Vector{Any}
    fit::Any
end  
struct  Kinbiont_res_odes_system
    method::String
    params::Vector{Any}
    fit::Any
end
struct Kinbiont_res_one_well
    method::String
    params::Vector{Any}
    fit::Any
    times::Any
end
struct Kinbiont_res_bootstrap_NL
    method::String
    params::Matrix{Any}
    fit::Any
    times::Any
    fin_param::Any
    new_param_fin::Any
    mean_param::Any
    sd_param::Any


end

struct Kinbiont_res_Log_Lin_files
    method::String
    params::Matrix{Any}
    fits::Tuple{Any}
    data::Tuple{Any}
    confidence_bands::Tuple{Any}
end

struct Kinbiont_res_one_file
    method::String
    params::Matrix{Any}
    fits::Tuple{Any}
    data::Tuple{Any}
end

struct Kinbiont_res_one_file_segmentation
    method::String
    params::Matrix{Any}
    fits::Tuple{Any}
    data::Tuple{Any}
    cp::Tuple{Any}
    vector_AIC::Any

end


struct Kinbiont_res_model_selection
    method::String
    params::Vector{Any}
    fit::Vector{Float64}
    times::Vector{Float64}
    rss_array::Any
    min_rss_array::Any
    param_min::Any
    min_AIC::Vector{Any}
    selected_model::String
    full_param::Vector{Any}
end


struct Kinbiont_res_NL_model_selection
    method::String
    params::Vector{Any}
    fit::Vector{Float64}
    times::Vector{Float64}
    score_res::Any
    top_loss::Any
end

struct Kinbiont_res_sensitivity_NL
    method::String
    params::Matrix{Any}
    fit::Any
    times::Any
    combinations::Matrix{Any}

end

struct Kinbiont_res_sensitivity
    method::String
    params::Matrix{Any}
    combinations::Matrix{Any}
end

struct Kinbiont_res_segmentation_ODE
    method::String
    params::Matrix{Any}
    fit::Array{Float64} 
    times::Array{Float64}
    interval_cdp::Array{Any}
    score_of_the_models::Any
end


struct Kinbiont_res_segmentation_NL
    method::String
    params::Matrix{Any}
    fit::Array{Float64} 
    times::Array{Float64}
    interval_cdp::Array{Any}
end






struct Kinbiont_model_to_fit
    model_list::Vector{Union{String,Function}} 
    # vector of models could be string if default kinbiont a string 
    # otherwise a function or an ode model ["logistic",NL_function,ODE_function]

    models_CI::Vector{Union{Nothing,Vector{Float64}}}
    # vector of starting guess of the models parameters  
    # [[p_1,p_2],[p_1,p_2],[p_1,p_2,p_3]] # mandatory

    Type_of_model::Vector{Union{Nothing,String}}
    # vector of type of model two options "ODE", "NL" same lenght model list, 
    # e.g., ["NL","NL","ODE"]
end

#=
CONSTRUCTOR:
function Kinbiont_data(data::Matrix{Float64}, times::Vector{Float64}, labels::Vector{String})
    Kinbiont_data(data, times, labels, nothing, nothing, nothing, nothing)
end
=#
mutable struct Kinbiont_data
    Data_matrix::Matrix{Float64} # vector of models could be string if default 
    # kinbiont a string otherwise a function or an ode model 
    # ["logistic",NL_function,ODE_function]
    Times::Vector{Float64} # vector of starting guess of the models parameters
    # [[p_1,p_2],[p_1,p_2],[p_1,p_2,p_3]] # mandatory
    Labels::Vector{String} # vector of starting upper bouds of the models 
    # parameters [[p_1,p_2],[p_1,p_2],[p_1,p_2,p_3]] # optional could be
    # missing e.g.,  [[missing],[p_1,p_2],[p_1,p_2,p_3]] or missing
    Clustering_labels::Union{Nothing, Vector{Int}} # vector of cluster ids,
    # one id per curve (with curve being a Data_matrix row)
    Clustering_wcss::Union{Nothing, Float64} # Within Cluster Sum of Squares
    # from the K-Means procedure
    Clustering_centroids_norm::Union{Nothing, Vector{Vector{Float64}}} # z-norm
    # centroids as computed by k-means (a K-length vector of vectors (K prototypes time points))
    Clustering_centroids_orig::Union{Nothing, Vector{Vector{Float64}}} # v
    # original space re-computed centroids (a K-length vector of vectors (K prototypes time points))
end

"""
Kinbiont_preprocessing_options:
Contains all configuration parameters for preprocessing and clustering of OD time-series.
"""
mutable struct Kinbiont_preprocessing_options
    # === General preprocessing flags ===
    smooth_data::Bool                     # whether to smooth curves
    type_of_smoothing::String             # "boxcar", "rolling_avg", "lowess", "gaussian", or "none"
    pt_avg::Int                           # window size for rolling average (if applicable)
    pt_smooth_derivative::Int             # used for derivative smoothing in downstream fits
    do_blank_subtraction::String          # "yes"/"no" or "auto"
    avg_replicate::Bool                   # whether to average replicates
    negative_corrections::Union{Float64, Nothing} # value to replace negatives, or nothing
    multiple_scattering_correction::Bool  # enable OD multiple scattering correction
    calibration_OD_curve::String          # path to calibration file (CSV)
    thr_lowess::Float64                   # bandwidth/frac parameter for LOWESS smoothing

    # === Gaussian smoothing parameters ===
    gaussian_h_mult::Union{Float64, Nothing}   # smoothing width multiplier
    gaussian_time_grid::Union{Vector{Float64}, Nothing} # optional target time grid

    # === Boxcar smoothing parameters ===
    pt_boxcar::Union{Int, Nothing}          # window length for boxcar filter

    # === Clustering parameters ===
    cluster_data::Bool                     # perform clustering (true/false)
    K::Union{Int,Nothing}                  # number of clusters
    q_low::Union{Float64, Nothing}          # lower quantile (robust min)
    q_high::Union{Float64, Nothing}         # upper quantile (robust max)
    tol_const::Union{Float64, Nothing}      # tolerance for constant cluster
    const_qbias::Union{Float64, Nothing}    # bias correction for constant cluster
    trend_factor::Union{Float64, Nothing}   # bias correction for growing curves
    use_exp_clust::Union{Bool, Nothing}     # enable exponential cluster
    rng::Union{Random.AbstractRNG, Nothing} # RNG for reproducible initialization
end

struct Kinbiont_generic_results
    Data_matrix::Matrix{Float64} # vector of models could be string if default kinbiont a string otherwise a function or an ode model ["logistic",NL_function,ODE_function]
    method::String
    params::Matrix{Any}
    fit::Vector{Vector{Float64}}
    times::Vector{Vector{Float64}}
    min_AIC::Vector{Any}
    selected_model::Vector{Float64}
    full_param::Vector{Any}
    Labels::Vector{Float64}              
end

#struct Kinbiont_options_fit
 #   loss_type="L2", 
  #  multistart=false,
   # correction_AIC=true,
   # n_restart=50,
   # optimizer=BBO_adaptive_de_rand_1_bin_radiuslimited(),
   # auto_diff_method=nothing,
   # cons=nothing,
   # beta_smoothing_ms=2.0,
   # opt_params...
#end
