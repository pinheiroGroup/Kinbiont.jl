#using DataDrivenDiffEq
#using ModelingToolkit: variable
#using DataDrivenSparse

function Kinbiont_data_generic_fitting(data_vector,
        times, 
        label, # data array could be a matrix A of any, a path in string 
               # "path_1", an array of matrix and string array 
               # ["path_1",Matrix_1,...]
        Kinbiont_model, #  Kinbiont fit model data struct 
        Kinbiont_options_fit
    )
 
    # unpack Kinbiont_model
    model_list = Kinbiont_model.model_list
    Type_of_model = Kinbiont_model.Type_of_model

    # initialize constants for AIC
    score_res = ["AIC"]
    top_score = 10^20
    top_model = Vector{Any}
    top_fitted_sol = Vector{Any}
    top_loss = 10^20
    temp_param = nothing
    count_custom = 0

    if !(Type_of_model in ["NL", "ODE","log_lin","AUC","clustering","DDDE"])
        @error("Model type $Type_of_model not recognized. Supported types are 'NL', 'ODE', 'log_lin', and 'AUC'.")
    end   

    if Type_of_model in ["NL", "ODE"]
        # grab model function and its string from model_list (unused)
        model_function, model_string = generic_model_selector(model_list)

        if model_string == "custom"
            count_custom += 1
            model_string = "custom_$(count_custom)"
        end

        temp_param = Kinbiont_model.models_CI
    end

    if Type_of_model == "ODE"
        res = Kinbiont_ODE_model_generic_fit(
            data_vector,
            times,
            model_string,
            model_function,
            Kinbiont_model,
            Kinbiont_options_fit,
            label,
            temp_param
        )
    end
 
    if Type_of_model == "NL"
        res = Kinbiont_NL_model_generic_fit(
            data_vector,
            times,
            model_string,
            model_function,
            Kinbiont_model,
            Kinbiont_options_fit,
            label,
            temp_param,
            
        )
    end

    if Type_of_model == "log_lin"
        res = Kinbiont_generic_log_lin(
            data_vector,
            times,
            Kinbiont_options_fit,
            label
        )
    end

    #=  
        Data Driven Differential Equation 
        (finds systems of equations from dataset)
    =#
    if Type_of_model == "DDDE"
        ddsol = Kinbiont_generic_DDDE(
            data_vector,
            times,
            Kinbiont_options_fit,
        )
        println("=== DataDrivenSolution summary ===")
        println(ddsol)

        # Extract discovered model and parameters
        system = get_basis(ddsol)
        params = get_parameter_map(system)

        println("\n=== Discovered system (basis) ===")
        println(system)

        println("\n=== Parameter map ===")
        println(params)
    end
    return res
end

# to be tested (DDDE: DataDrivenDiffEq.jl)
function Kinbiont_generic_DDDE(
        data_vector,
        times,
        Kinbiont_options_fit,
    )

    # 1) check data
    t = Float64.(times)
    y = Float64.(data_vector)
    N = length(t) # length of t (should be the same as y)
    if N != length(y) || N < 3
        error("DataDrivenDiffEq: check times and data length,
            should be the same and >= 3")
    end

    # 2) Numerical estimate of dy/dt (finite differences)
    dydt = similar(y) # construct buffer for derivative estimate of y

    for i in 1:N
        if i == 1
            # forward because t[0] doesn't exist
            dt = t[2] - t[1]
            dydt[1] = (y[2] - y[1]) / dt
        elseif i == N 
            # backward because t[N+1] doesn't exist
            dt = t[N] - t[N-1]
            dydt[N] = (y[N] - y[N-1]) / dt
        else
            # central t[i-1] and t[i+1] both exist
            dt = t[i+1] - t[i-1]
            dydt[i] = (y[i+1] - y[i-1]) / dt
        end
    end

    # 3) DataDriven Problem Construction (X -> Y, i.e., u -> du/dt)
    # DataDrivenDiffEq expects a matrix (dim x samples)
    # Raw Data seen as Matrix (num samples x samples):
    X = reshape(y, 1, :) 
    # Derivative of Data as Matrix (num samples x sample's derivatives):
    Y = reshape(dydt, 1, :)

    # degree of the polynomial and lambda (equation complexity tuner) for STLSQ
    # TO DO: add these fields to Kinbiont_options_fit
    max_degree = hasproperty(Kinbiont_options_fit, :ddde_max_degree) ?
        getproperty(Kinbiont_options_fit, :ddde_max_degree) : 4

    lambda_min = hasproperty(Kinbiont_options_fit, :ddde_lambda_min) ?
        getproperty(Kinbiont_options_fit, :ddde_lambda_min) : -5.0

    lambda_max = hasproperty(Kinbiont_options_fit, :ddde_lambda_max) ?
        getproperty(Kinbiont_options_fit, :ddde_lambda_max) : -1.0

    lambda_step = hasproperty(Kinbiont_options_fit, :ddde_lambda_step) ?
        getproperty(Kinbiont_options_fit, :ddde_lambda_step) : 0.5

    # generate lambdas from min to max using the step from options_fit
    lambdas = exp10.(lambda_min:lambda_step:lambda_max)

    # constructs the actual problem to solve
    problem = DirectDataDrivenProblem(X, Y)

    # 4) Polynomial base in u (OD)
    u = variable(:u)
    basis = Basis(polynomial_basis([u], max_degree), [u])

    # 5) STLSQ Solver (sparse regression for the ODE)
    opt = STLSQ(lambdas)
    common = DataDrivenCommonOptions(digits = 3)

    # ddsol is a ::DataDrivenSolution{Float64}
    ddsol = solve(problem, basis, opt; options = common)
    # from which we extract the symbolic model of the ODE, or its parameters...
    return ddsol
end

function generic_model_selector(model_function)

    if typeof(model_function) == String
        model_string = Kinbiont_models[model_function].name
        model_function = Kinbiont_models[model_string].func
    else
        model_string = "custom"    
    end  

    return model_function, model_string

end
    
function Kinbiont_ODE_model_generic_fit(
        data_vector,
        times,
        model_string,
        model_function,
        Kinbiont_model,
        Kinbiont_options_fit,
        label,
        temp_param,
    )

    # unpacking options fit
    type_of_loss = Kinbiont_options_fit.type_of_loss
    optimizer = Kinbiont_options_fit.optimizer
    integrator = Kinbiont_options_fit.integrator
    pt_smooth_derivative = Kinbiont_options_fit.pt_smooth_derivative
    auto_diff_method = Kinbiont_options_fit.auto_diff_method
    multistart=Kinbiont_options_fit.multistart
    n_restart=Kinbiont_options_fit.n_restart
    cons=Kinbiont_options_fit.cons
    opt_params=Kinbiont_options_fit.opt_params
    data = permutedim(hcat(times,data_vector))
    n_equation = detect_ode_dimension(model_function, temp_param; maxdim=100)

    res = fitting_one_well_custom_ODE(
        data, # dataset first row times second row OD
        model_string, # name of the well
        label, #label of the experiment
        model_function, # ode model to use
        temp_param,# initial guess param
        n_equation; # number ode in the system
        integrator=integrator, # selection of sciml integrator
        pt_avg=1, # numebr of the point to generate intial condition
        pt_smooth_derivative=pt_smooth_derivative,
        smoothing=false, # the smoothing is done or not?
        type_of_loss=type_of_loss, # type of used loss
        blank_array=zeros(100), # data of all blanks
        multiple_scattering_correction=false, # if true uses the given
                                              # calibration curve to fix data
        multistart=multistart,
        n_restart=n_restart,
        optimizer=optimizer,
        auto_diff_method=auto_diff_method,
        cons=cons,
        opt_params...
    )

    return res
end

function Kinbiont_NL_model_generic_fit(
        data_vector,
        times,
        model_string,
        model_function,
        Kinbiont_model,
        Kinbiont_options_fit,
        label,
        temp_param,
    )

    # unpacking options fit
    type_of_loss = Kinbiont_options_fit.type_of_loss
    optimizer = Kinbiont_options_fit.optimizer
    integrator = Kinbiont_options_fit.integrator
    pt_smooth_derivative = Kinbiont_options_fit.pt_smooth_derivative
    auto_diff_method = Kinbiont_options_fit.auto_diff_method
    multistart = Kinbiont_options_fit.multistart
    n_restart = Kinbiont_options_fit.n_restart
    cons = Kinbiont_options_fit.cons
    opt_params = Kinbiont_options_fit.opt_params
    data = permutedim(hcat(times,data_vector))  

    res = fit_NL_model(data, # dataset first row times second row OD
        model_string, # name of the well
        label, #label of the experiment
        model_function, #  model to use
        temp_param;# initial guess param    
        pt_smooth_derivative=pt_smooth_derivative,
        smoothing=false, # the smoothing is done or not?
        type_of_loss=type_of_loss, # type of used loss
        multiple_scattering_correction=false, # if true uses the given 
                                              # calibration curve to fix data
        penality_CI=3.0,
        optimizer=optimizer,
        multistart=multistart,
        n_restart=n_restart,
        auto_diff_method=auto_diff_method,
        cons=cons,
        opt_params...
    )

    return res

end

function Kinbiont_generic_log_lin(
        data_vector,
        times,
        Kinbiont_options_fit,
        label
    )

    # unpacking options fit
    pt_smoothing_derivative = Kinbiont_options_fit.pt_smooth_derivative
    pt_min_size_of_win = Kinbiont_options_fit.pt_min_size_of_win
    start_exp_win_thr = Kinbiont_options_fit.start_exp_win_thr
    type_of_win	 = Kinbiont_options_fit.type_of_win
    threshold_of_exp = Kinbiont_options_fit.threshold_of_exp

    data = permutedim(hcat(times,data_vector))

    res = fitting_one_well_Log_Lin(data,
        "log_lin",
        label;
        type_of_smoothing="NO",
        pt_smoothing_derivative=pt_smoothing_derivative,
        pt_min_size_of_win=pt_min_size_of_win,
        type_of_win=type_of_win,
        threshold_of_exp=threshold_of_exp,
        start_exp_win_thr=start_exp_win_thr,
    )

    return res

end

function detect_ode_dimension(odefun, param; maxdim=100)
    for n in 1:maxdim
        u  = zeros(n)
        du = zeros(n)
        try
            odefun(du, u, param, 0.0)
            return n
        catch err
            continue
        end
    end
    error("ERROR: Could not determine state dimension up to maxdim = $maxdim")
end