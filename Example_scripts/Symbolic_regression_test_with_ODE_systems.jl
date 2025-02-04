using Kinbiont
using DifferentialEquations
using CSV
using SymbolicRegression
using Plots
using StatsBase
using SymbolicRegression
using Distributions
# Generate a dataset with an unknown dependence on a feature 
function unknown_response(feature)

    response = (1 - feature)^2
    return response

end



function model_1(du, u, param, t)
    # Define the ODEs
    du[1] = param[1] * u[1] * u[4]
    du[2] = param[4] * du[1] - param[3] * u[2] - param[2] * u[2]
    du[3] = param[3] * u[2] - param[2] * u[3]
    du[4] = -du[1]
end
u0 = [0.1, 0.0, 0.0,1.0]  # Initial conditions for the variables

# Parameters
param = [0.1, 0.001, 0.5, 0.42]
lb1 =  [0.01, 0.0001, 0.1, 0.0]
ub1=  [0.2, 0.3, 1.1,1.0]
param_guess = lb1 .+ (ub1 .- lb1) ./ 2

param0 = param[4]
noise_value = 0.01

# defining the range of the perturbation on feature

feature_range = 0.0:0.1:2.0


plot(0, 0)
for f in feature_range

    # changing the parameters with unknown perturbation 
    param[4] = param0 * unknown_response(f)


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
                    string(f),
                    model_1, 
                    param_guess,
                    u0;
                    lb= lb1,
                   ub = ub1,
                  # set_of_equations_to_fit = [1,3,4]
    )

    display(plot!(fit[3]))


    if f == feature_range[1]
        results_fit = fit[2]
    else
        results_fit = vcat(results_fit, reduce(hcat,fit[2][2,:]))
    end



end

scatter(results_fit[2:end,1],results_fit[2:end,6],xlabel="Feature value", ylabel="Growth rate",)


# setting option for symbolic regression
options = SymbolicRegression.Options(
    binary_operators=[+, /, *, -],
    unary_operators=[square],
    constraints=nothing,
    elementwise_loss=nothing,
    loss_function=nothing,
    tournament_selection_n=12,
    tournament_selection_p=0.86,
    topn=12,
    complexity_of_operators=nothing,
    complexity_of_constants=nothing,
    complexity_of_variables=nothing,
    parsimony=0.05,
    dimensional_constraint_penalty=nothing,
    alpha=0.100000,
    maxsize=10,
    maxdepth=nothing
)
# generating feature matrix
# the first column is the label as a string of the feature value we used for the fitting labeling

feature_matrix = [[string(f),f] for f in feature_range]
feature_matrix = permutedims(reduce(hcat,feature_matrix))
results_fit[:,2] =results_fit[:,1] 
results_fit = permutedims(results_fit)
gr_sy_reg = Kinbiont.downstream_symbolic_regression(results_fit, feature_matrix, 6; options=options)

scatter(results_fit[2,2:end],results_fit[6,2:end,],xlabel="Feature value", ylabel="Growth rate",)
hline!(unique(gr_sy_reg[3][:, 1]), label=["Eq. 1" nothing], line=(3, :green, :dash))
plot!(unique(results_fit[2,2:end]), unique(gr_sy_reg[3][:, 2]), label=["Eq. 2" nothing], line=(3, :red))
plot!(unique(results_fit[2,2:end]), unique(gr_sy_reg[3][:, 3]), label=["Eq. 3" nothing], line=(3, :blue, :dashdot))
