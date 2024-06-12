struct Kimchi_res_one_well
    method::String
    fit::Vector{Float64}
    params::Vector{Any}
end

struct Kimchi_res_one_file
    method::String
    params::Matrix{Any}
end

Kimchi_one_file_output = Dict(model.name => model for model in models_list)
