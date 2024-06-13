struct Kimchi_res_one_well
    method::String
    fit::Vector{Float64}
    times::Vector{Float64}
    params::Vector{Any}
end

struct Kimchi_res_one_file
    method::String
    params::Matrix{Any}
end

struct Kimchi_res_model_selection
    method::String
    fit::Vector{Float64}
    times::Vector{Float64}
    params::Vector{Any}
    AIC::Vector{Float64}
    selected_model::String
    model_comparison::Matrix{Any}
    params_2::Vector{Any}
end

struct Kimchi_res_sensitivity
    method::String
    combinations::Matrix{Any}
    params::Matrix{Any}
end

struct Kimchi_res_segmentation_ODE
    method::String
    params::Matrix{Any}
    intervals_cdp::Array{Any}
    times_of_the_fit::Array{Float64}
    fit::Array{Float64} 
    total_loss::Float64

end