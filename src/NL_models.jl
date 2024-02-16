

struct NL_Model
    name::String
    func::Function
    params::Vector{String}
end



# from "Statistical evaluation of mathematical models for microbial growth"
# non linear model 
# exponential



function NL_model_exp(p, times)

    u = p[1] .* exp.(p[2] .* times)

    return u

end

function NL_model_logistic(p, times)

    u = p[1] ./ (1 .+ exp.(.-p[2] .* (times .- p[3])))

    return u

end



function NL_model_Gompertz(p, times)

    u = p[1] .* exp.(-exp.(-p[2] .* (times .- p[3])))

    return u

end


function NL_model_Bertalanffy(p, times)

    u = p[1] .+ (p[2] .- p[1]) .* (1 .- exp.(.-p[3] .* times)) .^ (1 ./ p[4])

    return u

end
#Richards

function NL_model_Richards(p, times)

    u = p[1] ./ (1 .+ p[2] .* exp.(.-p[3] .* (times .- p[4]))) .^ (1 ./ p[2])

    return u

end
#Morgan


function NL_model_Morgan(p, times)

    u = (p[1] .* p[2] .^ p[3] .+ p[4] .* times .^ p[3]) ./ (p[2] .^ p[3] .+ times .^ p[3])


    return u

end


#Weibull


function NL_model_Weibull(p, times)

    u = p[1] .- (p[1] .- p[2]) .* exp.(-(p[3] .* times) .^ p[4])


    return u

end



NL_models_list = [
    NL_Model(
        "NL_exponential",
        NL_model_exp,
        ["model", "well", "N0", "growth_rate", "th_max_gr", "emp_max_gr", "loss"]
    ),
    NL_Model(
        "NL_logistic",
        NL_model_logistic,
        ["model", "well", "N_max", "growth_rate", "lag", "th_max_gr", "emp_max_gr", "loss"]
    ),
    NL_Model(
        "NL_Gompertz",
        NL_model_Gompertz,
        ["model", "well", "N_max", "growth_rate", "lag", "th_max_gr", "emp_max_gr", "loss"]
    ),
    NL_Model(
        "NL_Bertalanffy",
        NL_model_Bertalanffy,
        ["model", "well", "N_0", "N_max", "growth_rate", "shape", "th_max_gr", "emp_max_gr", "loss"]
    ),
    NL_Model(
        "NL_Richards",
        NL_model_Richards,
        ["model", "well", "N_max", "shape", "growth_rate", "lag", "th_max_gr", "emp_max_gr", "loss"]
    ),
    NL_Model(
        "NL_Morgan",
        NL_model_Morgan,
        ["model", "well", "N_0", "K", "shape", "N_max", "th_max_gr", "emp_max_gr", "loss"]
    ),
    NL_Model(
        "NL_Weibull",
        NL_model_Weibull,
        ["model", "well", "N_max", "N_0", "growth_rate", "shape", "th_max_gr", "emp_max_gr", "loss"]
    ),
]

NL_models = Dict(model.name => model for model in NL_models_list)


