struct Kimchi_res_one_well
    method::String
    params::Vector{Any}
    fit::Any
    times::Any
end
struct Kimchi_res_bootstrap_NL
    method::String
    params::Matrix{Any}
    fit::Any
    times::Any
    fin_param::Any
    new_param_fin::Any
    mean_param::Any
    sd_param::Any


end

struct Kimchi_res_Log_Lin
    method::String
    params::Matrix{Any}
end

struct Kimchi_res_one_file
    method::String
    params::Matrix{Any}
    fits::Tuple{Any}
    data::Tuple{Any}
end

struct Kimchi_res_one_file_segmentation
    method::String
    params::Matrix{Any}
    fits::Tuple{Any}
    data::Tuple{Any}
    cp::Tuple{Any}
    vector_AIC::Any

end


struct Kimchi_res_model_selection
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


struct Kimchi_res_NL_model_selection
    method::String
    params::Vector{Any}
    fit::Vector{Float64}
    times::Vector{Float64}
    score_res::Any
    top_loss::Any
end

struct Kimchi_res_sensitivity_NL
    method::String
    params::Matrix{Any}
    fit::Any
    times::Any
    combinations::Matrix{Any}

end

struct Kimchi_res_sensitivity
    method::String
    params::Matrix{Any}
    combinations::Matrix{Any}
end

struct Kimchi_res_segmentation_ODE
    method::String
    params::Matrix{Any}
    fit::Array{Float64} 
    times::Array{Float64}
    interval_cdp::Array{Any}
    score_of_the_models::Any
end


struct Kimchi_res_segmentation_NL
    method::String
    params::Matrix{Any}
    fit::Array{Float64} 
    times::Array{Float64}
    interval_cdp::Array{Any}
end



