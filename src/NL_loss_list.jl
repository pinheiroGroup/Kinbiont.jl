"
List of loss functions to fit ODE
"


function NL_Log_L2(data,model_function,u, p)

    model = model_function(u, data[1, :])
    n_data = length(data[1, :])
    residuals = (model .- data[2, :]) ./ n_data
    return log(sum(residuals .^ 2))
end


function NL_L2(data,model_function,u, p)

    model = model_function(u, data[1, :])
    n_data = length(data[1, :])
    residuals = (model .- data[2, :]) ./ n_data
    return (sum(residuals .^ 2))
end

function NL_RE(data,model_function,u, p)

    model = model_function(u, data[1, :])
    n_data = length(data[1, :])
    
    return (0.5/n_data).* sum( (1 .-data[2,:]./model).^2)
end

function NL_Log_RE(data,model_function,u, p)

    model = model_function(u, data[1, :])
    n_data = length(data[1, :])

    return log((0.5/n_data).*sum( (1 .-data[2,:]./model).^2))
end






function NL_RE_fixed_CI(data,model_function,pen, u, p)

    model = model_function(u, data[1, :])
    n_data = length(data[1, :])
    penality_ci = n_data/pen

    re_ci = (0.5/n_data).* sum( (1 .-data[2,1]./model[1]).^2)
    re = (0.5/n_data).* sum( (1 .-data[2,2:end]./model[2:end]).^2)
    re_t = penality_ci * re_ci + re
    return re_t
end


function NL_L2_fixed_CI(data,model_function,pen, u, p)

    model = model_function(u, data[1, :])
    n_data = length(data[1, :])
    penality_ci = n_data/pen

    residuals_ci = penality_ci .* (model[1] .- data[2, 1]) 
    residuals = (model[2:end] .- data[2, 2:end]) ./ (n_data -1)
    residuals_tot = (sum(residuals_ci .^ 2)) + (sum(residuals .^ 2))
    return residuals_tot
end


function select_loss_function_NL(loss_name, data,pen, model_function)
    loss_functions = Dict(
        "L2" => NL_L2,
        "RE" => NL_RE,
        "L2_log" => NL_Log_L2,
        "RE_log" => NL_Log_RE,
        "L2_fixed_CI" => NL_L2_fixed_CI,
        "RE_fixed_CI" => NL_RE_fixed_CI)



        if loss_name == "L2_fixed_CI" || loss_name == "RE_fixed_CI"
            return (u,p) -> loss_functions[loss_name](data, model_function, pen, u, p)
        
        else
        return (u,p) -> loss_functions[loss_name](data, model_function, u, p)
end
end
