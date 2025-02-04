using Kinbiont
using DifferentialEquations
using OptimizationBBO
using NaNMath
using Plots
using Distributions
# examples from Genome-Scale Reconstruction of Microbial Dynamic Phenotype: Successes and Challenges
# Extremely High Proteins Concentration: Macromolecular Crowding
function enzyme_aggregation(du, u, param, t)
    e, x, y, m = u
    k1, k2, k3, k4, k_cat, n = param
    
    du[1] = k4 * x - k3 * m * e + k2 * y^n - k1 * e + k_cat * x  # Free active enzymes
    du[2] = k3 * m * e - k4 * x - k_cat * x                      # Enzyme-substrate complex
    du[3] = k1 * e - k2 * y^n                                    # Inactive aggregates
    du[4] = -du[1]                                              # Substrate degradation rate
end

u0 = [1.0, 0.1, 0.1, 1.0]  # Initial conditions [e, x, y, m]
param = [0.1, 0.1, 0.05, 0.05, 0.02, 2]  # [k1, k2, k3, k4, k_cat, n, e0]


Simulation =  ODEs_system_sim(
    enzyme_aggregation, #string of the model
    u0, # starting condition
    0.0, # start time of the sim
    30.0, # final time of the sim
    1.0, # delta t for poisson approx
    param; # parameters of the ODE model
)

# Plot the data


scatter(Simulation)


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
                enzyme_aggregation, 
                param,
                u0;

)

plot!(fit[3])

