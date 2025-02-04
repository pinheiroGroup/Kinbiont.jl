using Kinbiont
using DifferentialEquations
using CSV
using SymbolicRegression
using Plots
using StatsBase
using SymbolicRegression
using Distributions
# Generate a dataset with an unknown dependence on a feature 



function model_1(du, u, param, t)
    # Define the ODEs
    du[1] = param[1] * u[1] * u[4]
    du[2] = param[4] * du[1] - param[3] * u[2] - param[2] * u[2]
    du[3] = param[3] * u[2] - param[2] * u[3]
    du[4] = -du[1]
end

u0 = [0.1, 0.0, 0.0,1.0]  # Initial conditions for the variables

# Parameters
param = [0.1, 0.01, 0.5, 0.42]
lb1 =  [0.01, 0.0001, 0.0, 0.01]
ub1=  [0.2, 0.3, 1.1, 1.0]
param_guess = lb1 .+ (ub1 .- lb1) ./ 2

param0 = param[4]
noise_value = 0.01

 # Calling the simulation function
Simulation =  ODEs_system_sim(
    model_1, #string of the model
    u0, # starting condition
    0.0, # start time of the sim
    50.0, # final time of the sim
    1.0, # delta t for poisson approx
    param; # parameters of the ODE model
)
# Plotting scatterplot of data without noise

#adding uniform random noise
sol_time = reduce(hcat, Simulation.t)
sol_t = reduce(hcat,Simulation.u)
noise_unifom = rand(Uniform(-0.05,0.05),size(sol_t)[2])
sol_t_noise = [sol_t[i,:] .+ rand(Uniform(-0.05,0.05),size(sol_t)[2]) for i in 1:size(sol_t)[1]]
sol_t_noise =permutedims(reduce(hcat,sol_t_noise))

data = vcat(sol_time,sol_t_noise)
# Plot data with noise
display(scatter(data[1,:],data[2,:]))
display( scatter!(data[1,:],data[3,:]))
display(scatter!(data[1,:],data[4,:]))
display( scatter!(data[1,:],data[5,:]))



fit = fit_ODEs_System(data,
                "test",
                model_1, 
                param_guess,
                u0;
                lb= lb1,
               ub = ub1
)

plot!(fit[3])




