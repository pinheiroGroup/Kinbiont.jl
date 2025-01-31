using Kinbiont
using DifferentialEquations
using OptimizationBBO
using NaNMath
using Plots
using Distributions

# generate Simulated data with SIR model
Simulation =  ODEs_system_sim(
    "SIR", #string of the model
   [0.9,0.1,0.0], # starting condition
    0.0, # start time of the sim
    30.0, # final time of the sim
    1.0, # delta t for poisson approx
    [0.5,0.3]; # parameters of the ODE model
)


# Plot the data


scatter(Simulation)

# Add the noise and formating the data as required by Kibiont


#adding uniform random noise

sol_time = reduce(hcat, Simulation.t)
sol_t = reduce(hcat,Simulation.u)
noise_unifom = rand(Uniform(-0.05,0.05),size(sol_t)[2])
sol_t_noise = [sol_t[i,:] .+ rand(Uniform(-0.05,0.05),size(sol_t)[2]) for i in 1:size(sol_t)[1]]
sol_t_noise =permutedims(reduce(hcat,sol_t_noise))

data = vcat(sol_time,sol_t_noise)
# Plot data with noise
scatter(data[1,:],data[2,:])
scatter!(data[1,:],data[3,:])
scatter!(data[1,:],data[4,:])

Start_IC = [0.9,0.1,0.0]
# fit all dataset  

fit = fit_ODEs_System(data,
                    "test",
                    "SIR", 
                    [0.1,0.5],
                    Start_IC;
)

plot!(fit[3])

# Now we remove a data form the system the measura=ment of R
data_reduced = hcat(data[1,:],data[2,:])
data_reduced = permutedims(hcat(data_reduced,data[3,:]))


fit = fit_ODEs_System(data_reduced,
                    "test",
                    "SIR", 
                    [0.1,0.5],
                    Start_IC;
                    set_of_equation_to_fit = [1,2]
)



scatter(data[1,:],data[2,:])
scatter!(data[1,:],data[3,:])
scatter!(data[1,:],data[4,:])
plot!(fit[3])


