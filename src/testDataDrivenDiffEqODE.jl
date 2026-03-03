using DataDrivenDiffEq
using ModelingToolkit: variable
using DataDrivenSparse
using Random
using Plots

#============================================================
 Code used to test the following function on synthetic data
============================================================#

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

# ============================================================
# Synthetic logistic OD data
# ============================================================
"""
Generate synthetic logistic OD data:

    y'(t) = r * y * (1 - y/K)

Return:
  y_noisy::Vector{Float64}   measured OD (with noise)
  t_vec::Vector{Float64}     time points
  y_clean::Vector{Float64}   noiseless logistic curve
"""
function make_logistic_OD_data(;
    r::Float64 = 0.8,
    K::Float64 = 1.0,
    y0::Float64 = 0.05,
    t_end::Float64 = 10.0,
    N::Int = 201,
    noise_std::Float64 = 0.01,
    rng = Random.default_rng()
)
    t_vec = collect(range(0.0, t_end; length = N))
    y_clean = Vector{Float64}(undef, N)

    for (i, ti) in enumerate(t_vec)
        # analytic logistic solution
        y_clean[i] = K / (1 + (K / y0 - 1) * exp(-r * ti))
    end

    y_noisy = y_clean .+ noise_std .* randn(rng, N)
    return y_noisy, t_vec, y_clean
end

# ============================================================
# Minimal options object for DDDE
# (hasproperty(...) works on NamedTuple)
# ============================================================
Kinbiont_options_fit_DDDE = (
    ddde_max_degree = 4,   # max polynomial degree in y
    ddde_lambda_min = -5.0,
    ddde_lambda_max = -1.0,
    ddde_lambda_step = 0.5,
    ddde_digits = 3        # rounding in output model
)

# ============================================================
# Test: apply Kinbiont_generic_DDDE to synthetic logistic data
# ============================================================
function demo_DDDE_on_logistic()
    rng = Random.MersenneTwister(1234)

    y_noisy, t_vec, y_clean = make_logistic_OD_data(
        r = 0.8,
        K = 1.0,
        y0 = 0.05,
        t_end = 10.0,
        N = 201,
        noise_std = 0.01,
        rng = rng
    )

    # Call your DDDE routine (label is not needed here, so we omit it)
    ddsol = Kinbiont_generic_DDDE(
        y_noisy,
        t_vec,
        Kinbiont_options_fit_DDDE,
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

    # Optional: quick plot of clean vs noisy data
    plt = plot(
        t_vec, y_clean;
        label = "clean logistic",
        xlabel = "time",
        ylabel = "OD",
        title = "Synthetic logistic OD (clean vs noisy)",
    )
    plot!(plt, t_vec, y_noisy; seriestype = :scatter, label = "noisy data")
    display(plt)

    return ddsol
end

# If you run this as a script: execute the demo
demo_DDDE_on_logistic()
