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










function select_loss_function_NL(loss_name, data, model_function)
    loss_functions = Dict(
        "L2" => NL_L2,
        "RE" => NL_RE,
        "L2_log" => NL_Log_L2,
        "RE_log" => NL_Log_RE,)


        return (u,p) -> loss_functions[loss_name](data, model_function, u, p)
end


