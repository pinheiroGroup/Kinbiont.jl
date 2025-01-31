using Kinbiont
using DifferentialEquations
using OptimizationBBO
using Plots
using Catalyst
using DiffEqParamEstim



u0 = [:S => 301, :E => 100, :SE => 0, :P => 0]
ps = [:kB => 0.00166, :kD => 0.0001, :kP => 0.1]

model_Michaelis_Menten_enzyme_kinetics= @reaction_network begin
    kB, S + E --> SE
    kD, SE --> S + E
    kP, SE --> P + E
end 






Simulation =  Kinbiont.Kinbiont_Reaction_network_sim("Michaelis_Menten",
    u0,
   0.0, # start time of the sim
   10.0, # final time of the sim
    0.1, # delta t for poisson approx
    ps; # parameters of the ODE model
    )

plot(Simulation)


sol_time = reduce(hcat, Simulation.t)
sol_t = reduce(hcat,Simulation.u)
noise_unifom = rand(Uniform(-0.01,0.05),size(sol_t)[2])
sol_t_noise = [sol_t[i,:] .+ rand(Uniform(-0.05,10.05),size(sol_t)[2]) for i in 1:size(sol_t)[1]]
sol_t_noise =permutedims(reduce(hcat,sol_t_noise))

data = vcat(sol_time,sol_t_noise)
scatter(data[1,:],data[2,:])
scatter!(data[1,:],data[3,:])
scatter!(data[1,:],data[4,:])
scatter!(data[1,:],data[5,:])


fit = RN_fit(data, 
    model_Michaelis_Menten_enzyme_kinetics,
    u0,
    [:kB => 0.00166, :kD => 0.0001, :kP => 0.1];
)

plot!(fit[4])

# example on how create a custum model


model_Glycolysis = @reaction_network begin
    (kf1, kr1), Glucose + ATP <--> Glucose6P + ADP
    (kf2, kr2), Glucose6P <--> Fructose6P
    (kf3, kr3), Fructose6P + ATP <--> Fructose16BP + ADP
    (kf4, kr4), Fructose16BP <--> DHAP + GAP
    (kf5, kr5), DHAP <--> GAP
    (kf6, kr6), GAP + NADplus + Pi <--> BPG13 + NADH + Hplus
    (kf7, kr7), BPG13 + ADP <--> PG3 + ATP
    (kf8, kr8), PG3 <--> PG2
    (kf9, kr9), PG2 <--> PEP + H2O
    (kf10, kr10), PEP + ADP <--> Pyruvate + ATP
end

u0 = [
    :Glucose => 100, 
    :ATP => 200, 
    :Glucose6P => 0, 
    :Fructose6P => 0, 
    :Fructose16BP => 0, 
    :DHAP => 0, 
    :GAP => 0, 
    :NADplus => 100, 
    :NADH => 0, 
    :Pi => 100, 
    :BPG13 => 0, 
    :PG3 => 0, 
    :PG2 => 0, 
    :PEP => 0, 
    :Pyruvate => 0, 
    :ADP => 0, 
    :H2O => 0, 
    :Hplus => 0
]


ps = [
    :kf1 => 0.0011, :kr1 => 0.005, 
    :kf2 => 0.005, :kr2 => 0.002, 
    :kf3 => 0.02, :kr3 => 0.01, 
    :kf4 => 0.015, :kr4 => 0.007, 
    :kf5 => 0.01, :kr5 => 0.005, 
    :kf6 => 0.02, :kr6 => 0.01, 
    :kf7 => 0.03, :kr7 => 0.015, 
    :kf8 => 0.01, :kr8 => 0.005, 
    :kf9 => 0.008, :kr9 => 0.004, 
    :kf10 => 0.04, :kr10 => 0.02
]


Simulation =  Kinbiont.Kinbiont_Reaction_network_sim(model_Glycolysis,
    u0,
   0.0, # start time of the sim
   30.0, # final time of the sim
    0.1, # delta t for poisson approx
    ps; # parameters of the ODE model
    )

plot(Simulation)
