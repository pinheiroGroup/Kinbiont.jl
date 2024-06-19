



using Kimchi
using Optimization
using OptimizationOptimJL
using OptimizationMOI
using OptimizationBBO
using OptimizationNLopt
using Ipopt
 using Random
using Distributions
using DifferentialEquations
# model 1 logistic
model = "logistic"

ub_log = [0.1,2.0]
lb_log = [0.01,0.4]

u0= [0.1]
n_start= [0.1]
noise_value = 0.05


tstart = 0.0
tmax = 600.0
delta_t = 10.0

param_of_ode = [0.0,0.9]
param_of_ode[1] = rand(Uniform(lb_log[1],ub_log[1]),1)[1]
sim = ODE_sim(model, n_start, tstart, tmax, delta_t, param_of_ode)
noise_unifom = rand(Uniform(-noise_value,noise_value),length(sim.t))


data_t = reduce(hcat,sim.t)
data_o = reduce(hcat,sim.u)
data_OD = vcat(data_t,data_o)
data_OD[2,:] = data_OD[2,:] .+ noise_unifom

# symbolic problem
param_0 = lb_log .+ (ub_log .- lb_log)./2
ODE_prob = ODEProblem(logistic, u0, (data_OD[1,1],data_OD[1,end]), param_0)

loss_function = select_loss_function("RE", data_OD, ODE_prob, Tsit5(), data_OD[1,:], [0.0])

function KimchiSolve( loss_function ,
                      u0,
                      p;
                      opt,
                      auto_diff_method = nothing,
                      cons = nothing,
                      opt_params...)

    # generation of optimization fuction

    if isnothing(auto_diff_method) == true && isnothing(cons) == true

        optf = Optimization.OptimizationFunction((x, p) -> loss_function(x))

    elseif isnothing(auto_diff_method) == false && isnothing(cons) == true

        optf = Optimization.OptimizationFunction((x, p) -> loss_function(x),auto_diff_method)


    elseif isnothing(auto_diff_method) == true && isnothing(cons) == false

        optf = Optimization.OptimizationFunction((x, p) -> loss_function(x), cons = cons)


    else isnothing(auto_diff_method) == false && isnothing(cons) == false
        
        
        optf = Optimization.OptimizationFunction((x, p) -> loss_function(x),auto_diff_method, cons = cons)


    end   

    
    prob = OptimizationProblem(optf, p, u0; opt_params...)
    
    sol = solve(prob, opt)
    
    return sol
end

# resolution 1
KimchiSolve(loss_function, u0, p; opt=BBO_adaptive_de_rand_1_bin_radiuslimited(), lb = lb_log, ub = ub_log)
# resolution 2
KimchiSolve(loss_function, u0, p; opt=NelderMead(),maxiters=10,sense=4)
# resolution 3
KimchiSolve(loss_function, u0, p; opt=NLopt.LN_PRAXIS(),   abstol =0.01)
# resolution 4
KimchiSolve(loss_function, u0, p; opt=BFGS(),   auto_diff_method = Optimization.AutoFiniteDiff() ,abstol =0.01)
# resolution 5
cons_(res, x, p) = (res .= [p[1], p[1]])
KimchiSolve(loss_function, u0, p; opt=IPNewton(), auto_diff_method = nothing, cons = cons_, lcons = [-1.0, -1.0], ucons = [0.8, 2.0] )





function KimchiSolve( loss_function ,
    u0,
    p;
    opt,
    auto_diff_method = nothing,
    cons = nothing,
    opt_params...)

# generation of optimization fuction

if isnothing(auto_diff_method) == true && isnothing(cons) == true

optf = Optimization.OptimizationFunction((x, p) -> loss_function(x))

elseif isnothing(auto_diff_method) == false && isnothing(cons) == true

optf = Optimization.OptimizationFunction((x, p) -> loss_function(x),auto_diff_method)


elseif isnothing(auto_diff_method) == true && isnothing(cons) == false

optf = Optimization.OptimizationFunction((x, p) -> loss_function(x), cons = cons)


else isnothing(auto_diff_method) == false && isnothing(cons) == false


optf = Optimization.OptimizationFunction((x, p) -> loss_function(x),auto_diff_method, cons = cons)


end   


prob = OptimizationProblem(optf, p, u0; opt_params...)

sol = solve(prob, opt)

return sol
end

function KimchiSolve2( loss_function ,
    u0,
    p;
    opt,
    opt_func_params...,
    opt_params...)

# generation of optimization fuction


optf = Optimization.OptimizationFunction((x, p) -> loss_function(x),opt_func_params)



prob = OptimizationProblem(optf, p, u0; opt_params...)

sol = solve(prob, opt)

return sol
end


# resolution 1
KimchiSolve2(loss_function, u0, p; opt=BBO_adaptive_de_rand_1_bin_radiuslimited(), lb = lb_log, ub = ub_log)
# resolution 2
KimchiSolve2(loss_function, u0, p; opt=NelderMead(),maxiters=10,sense=4)
# resolution 3
KimchiSolve2(loss_function, u0, p; opt=NLopt.LN_PRAXIS(),   abstol =0.01)\
# resolution 4
KimchiSolve2(loss_function, u0, p; opt=BFGS(),   auto_diff_method = Optimization.AutoFiniteDiff() ,abstol =0.01)
# resolution 5
cons_(res, x, p) = (res .= [p[1], p[1]])
KimchiSolve2(loss_function, u0, p; opt=IPNewton(), auto_diff_method = nothing, cons = cons_, lcons = [-1.0, -1.0], ucons = [0.8, 2.0] )
