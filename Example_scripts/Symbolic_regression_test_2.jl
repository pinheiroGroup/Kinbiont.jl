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


# defining the used ODE model 
results_fit =  Any

ODE_models = "baranyi_richards"

ub_1 = [0.1, 5.1, 500.0, 5.0]
lb_1 = [0.0001, 0.1, 0.00, 0.2]
p1_guess = lb_1 .+ (ub_1 .- lb_1) ./ 2


# defining the range of the perturbation on feature

feature_range = 0.0:0.2:2.0

# defining the parameters values for the simulation 
p_sim = [0.05, 1.0, 50.0, 1.0]
psim_1_0 =  p_sim[1]
t_min = 0.0
t_max = 800.0
n_start = [0.1]
delta_t = 5.0
noise_value = 0.03

plot(0, 0)
for f in feature_range

    # changing the parameters with unknown perturbation 
    p_sim[1] = psim_1_0 * unknown_response(f) .+ 0.01


    # Calling the simulation function
    sim = Kinbiont.ODE_sim("baranyi_richards", n_start, t_min, t_max, delta_t, p_sim)

    # Plotting scatterplot of data without noise

    #adding uniform random noise
    noise_unifom = rand(Uniform(-noise_value, noise_value), length(sim.t))


    data_t = reduce(hcat, sim.t)
    data_o = reduce(hcat, sim.u)
    data_OD = vcat(data_t, data_o)
    data_OD[2, :] = data_OD[2, :] .+ noise_unifom
    # ploting scatterplot of data with noise

    display(Plots.scatter!(data_OD[1, :], data_OD[2, :], xlabel="Time", ylabel="Arb. Units", label= nothing, color=:red, markersize=2, size=(300, 300)))

    results_ODE_fit = fitting_one_well_ODE_constrained(
        data_OD,
        string(f),
        "test_ODE",
        "baranyi_richards",
        p1_guess;
        lb=lb_1,
        ub=ub_1
    )
    display(Plots.plot!(results_ODE_fit[4], results_ODE_fit[3], xlabel="Time", ylabel="Arb. Units", label= nothing, color=:red, markersize=2, size=(300, 300)))


    if f == feature_range[1]
        results_fit = results_ODE_fit[2]
    else
        results_fit = hcat(results_fit, results_ODE_fit[2])
    end



end

scatter(results_fit[2,:],results_fit[4,:,],xlabel="Feature value", ylabel="Growth rate",)


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


gr_sy_reg = Kinbiont.downstream_symbolic_regression(results_fit, feature_matrix, 4; options=options)

scatter(results_fit[2,:],results_fit[4,:,],xlabel="Feature value", ylabel="Growth rate",)
hline!(unique(gr_sy_reg[3][:, 1]), label=["Eq. 1" nothing], line=(3, :green, :dash))
plot!(unique(results_fit[2,:]), unique(gr_sy_reg[3][:, 2]), label=["Eq. 2" nothing], line=(3, :red))
plot!(unique(results_fit[2,:]), unique(gr_sy_reg[3][:, 3]), label=["Eq. 3" nothing], line=(3, :blue, :dashdot))

