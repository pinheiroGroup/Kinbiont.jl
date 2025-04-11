struct Kinbiont_ODEs_system
    name::String
    func::Function
    params::Vector{String}
end

function SIR(du, u, param, t)
    du[1] = -u[1] * u[2] * param[1]
    du[2] = u[1] * u[2] * param[1] - u[2] * param[2]
    du[3] = u[2] * param[2]
end

function SIR_BD(du, u, param, t)
    du[1] = -u[1] * u[2] * param[1] + param[3] * u[1] - u[1] * param[4]
    du[2] = u[1] * u[2] * param[1] - u[2] * param[2] - u[2] * param[4]
    du[3] = u[2] * param[2] - u[3] * param[4]
end

function SIR_B(du, u, param, t)
    du[1] = -u[1] * u[2] * param[1] + param[3] * u[1]
    du[2] = u[1] * u[2] * param[1] - u[2] * param[2]
    du[3] = u[2] * param[2] - u[3] * param[4]
end

function SIS(du, u, param, t)
    du[1] = -u[1] * u[2] * param[1] + param[2] * u[2]
    du[2] = u[1] * u[2] * param[1] - u[2] * param[2]
end

function Lotka_Volterra(du, u, param, t)
    du[1] = u[1] * param[1] - param[2] * u[1] * u[2]
    du[2] = -param[3] * u[2] + param[4] * u[2] * u[1]
end

function Lotka_Volterra_with_substrate(du, u, param, t)
    du[1] = u[1] * param[1] * (u[3] / (u[3] + param[5])) - param[2] * u[1] * u[2]
    du[2] = -param[3] * u[2] + param[4] * u[2] * u[1]
    du[3] = -u[1] * param[1] * (u[3] / (u[3] + param[5]))

end




function monod_chemostat(du, u, param, t)
    # State variables:
    # u[1] = x (Biomass concentration)
    # u[2] = s (Substrate concentration)

    # Parameters:
    # param[1] = Ks      (Substrate affinity constant)
    # param[2] = m       (Maintenance coefficient)
    # param[3] = Ymax    (Maximum yield coefficient)
    # param[4] = mum     (Maximum specific growth rate)
    # param[5] = D       (Dilution rate)
    # param[6] = sr      (Inflow substrate concentration)

    x = u[1]  # Biomass concentration
    s = u[2]  # Substrate concentration

    # Critical substrate concentration s*
    s_star = param[1] * param[2] * param[3] / param[4]

    # Growth rate μ
    μ = param[4] * (s - s_star) / (param[1] + s)

    # Yield Y
    Y = (param[3] * param[5]) / (param[5] + param[2] * param[3])

    # Biomass dynamics
    du[1] = μ * x - param[5] * x

    # Substrate dynamics
    du[2] = param[5] * (param[6] - s) - (μ * x / Y) - param[2] * x
end



function monod_ierusalimsky(du, u, param, t)
    # State variables:
    # u[1] = x (Biomass concentration)
    # u[2] = s (Substrate concentration)
    # u[3] = p (Product concentration)

    # Parameters:
    # param[1] = Ks      (Substrate affinity constant)
    # param[2] = Kp      (Product inhibition constant)
    # param[3] = mum     (Maximum specific growth rate)
    # param[4] = Ymax    (Maximum yield coefficient)
    # param[5] = Yp      (Product yield coefficient)
    # param[6] = D       (Dilution rate)
    # param[7] = sr      (Inflow substrate concentration)
    # param[8] = m       (Maintenance coefficient)

    x = u[1]  # Biomass concentration
    s = u[2]  # Substrate concentration
    p = u[3]  # Product concentration

    # Growth rate μ with Monod-Ierusalimsky kinetics
    μ = param[3] * (s / (param[1] + s)) * (param[2] / (param[2] + p))

    # Yield Y
    Y = (param[4] * param[6]) / (param[6] + param[8] * param[4])

    # Biomass dynamics
    du[1] = μ * x - param[6] * x

    # Substrate dynamics
    du[2] = param[6] * (param[7] - s) - (μ * x / Y) - param[8] * x

    # Product dynamics
    du[3] = param[5] * μ * x - param[6] * p
end


function droop_model(du, u, param, t)
    # State variables:
    # u[1] = X (Biomass concentration)
    # u[2] = S (Substrate concentration)
    # u[3] = Q (Intracellular quota)

    # Parameters:
    # param[1] = mum     (Maximum growth rate, μm)
    # param[2] = rho_m   (Maximum nitrogen absorption rate, ρm)
    # param[3] = Ks      (Half-saturation constant, Ks)
    # param[4] = D       (Dilution rate, D)
    # param[5] = Sin     (Inflow substrate concentration, S_in)
    # param[6] = Q0      (Minimum cell quota, Q0)

    X = u[1]  # Biomass concentration
    S = u[2]  # Substrate concentration
    Q = u[3]  # Intracellular quota

    # Uptake rate ρ(S) based on Michaelis-Menten kinetics
    rho = param[2] * (S / (S + param[3]))

    # Growth rate μ(Q) as a function of intracellular quota Q
    mu = param[1] * (1 - param[6] / Q)

    # Biomass dynamics
    du[1] = mu * X - param[4] * X

    # Substrate dynamics
    du[2] = rho * X - param[4] * S + param[4] * param[5]

    # Intracellular quota dynamics
    du[3] = rho - mu * Q
end
function Synthetic_chemostat_model(du, u, p, t)
    S, N, R = u  # state variables: Substrate, Biomass, Reporter
    
    D, sr,  k_r, k_s, k_s_prime, Q, Q_prime, Y, a_0 = p  # Parameters
    
    # Calculate the substrate consumption rate q_s
    q_s_value = R * (Q * S) / (k_s + S) + (1 - R) * (Q_prime * S) / (k_s_prime + S)
    
    # Growth rate μ
    μ_value = Y * q_s_value - a_0 * R
    
    # Differential equations
    du[1] = D * (sr - S) - q_s_value * N  # dS/dt
    du[2] = μ_value * N - D * N  # dN/dt
    du[3] = μ_value * (S / (k_r + S) - R)  # dR/dt
end

# Define the Synthetic Batch Model ODEs
function  Synthetic_batch_model(du, u, p, t)
    S, N, R = u  # state variables: Substrate, Biomass, Reporter
    
    k_r, k_s, k_s_prime, Q, Q_prime, Y, a_0 = p  # Parameters
    
    # Calculate the substrate consumption rate q_s
    q_s_value = R * (Q * S) / (k_s + S) + (1 - R) * (Q_prime * S) / (k_s_prime + S)
    
    # Growth rate μ
    μ_value = Y * q_s_value - a_0 * R
    
    # Differential equations
    du[1] = -q_s_value * N  # dS/dt
    du[2] = μ_value * N  # dN/dt
    du[3] = μ_value * (S / (k_r + S) - R)  # dR/dt
end





# Monod Chemostat Model


ODEs_system_list = [
    Kinbiont_ODEs_system(
        "Monod_Chemostat",  # Model name
        monod_chemostat,    # The ODE function
        ["Ks", "m", "Ymax", "mum", "D", "sr"]  # Parameters list
    ),
    Kinbiont_ODEs_system(
        "Synthetic_Chemostat",  # Model name
        Synthetic_chemostat_model,    # The ODE function
        ["D", "sr", "Kr", "Ks", "ks'", "Q", "Q'", "Y", "a_0"]  # Parameters list
    ),
    Kinbiont_ODEs_system(
        "Synthetic_batch",  # Model name
        Synthetic_batch_model,    # The ODE function
        [ "Kr", "Ks", "ks'", "Q", "Q'", "Y", "a_0"]  # Parameters list
    ),
    # Monod Ierusalimsky Model
    Kinbiont_ODEs_system(
        "Monod_Ierusalimsky",  # Model name
        monod_ierusalimsky,    # The ODE function
        ["Ks", "Kp", "mum", "Ymax", "Yp", "D", "sr", "m"]  # Parameters list
    ),

    # Droop Model
    Kinbiont_ODEs_system(
        "Droop",  # Model name
        droop_model,  # The ODE function
        ["mum", "rho_m", "Ks", "D", "Sin", "Q0"]  # Parameters list
    ),

    Kinbiont_ODEs_system(
        "SIR",
        SIR,
        ["Infection_rate", "recovery_rate"]
    ),
    Kinbiont_ODEs_system(
        "SIR_BD",
        SIR_BD,
        ["Infection_rate", "recovery_rate", "Birth_rate", "Death_rate"]
    ),
    Kinbiont_ODEs_system(
        "SIR_B",
        SIR_B,
        ["Infection_rate", "recovery_rate", "Birth_rate"]
    ),
    Kinbiont_ODEs_system(
        "SIS",
        SIS,
        ["Infection_rate", "recovery_rate"]
    ), Kinbiont_ODEs_system(
        "Lotka_Volterra",
        Lotka_Volterra,
        ["prey_GR", "Prey_DR", "Predator_DR", "Predator_GR"]
    ), Kinbiont_ODEs_system(
        "Lotka_Volterra_with_substrate",
        Lotka_Volterra_with_substrate,
        ["prey_GR", "Prey_DR", "Predator_DR", "Predator_GR", "Mu_max", "Monod_constant"]
    ),]


ODEs_system_models = Dict(ODEs_system_list.name => ODEs_system_list for ODEs_system_list in ODEs_system_list)

export SIR  
export SIR_BD  
export SIR_B  
export SIS  
export Lotka_Volterra  
export Lotka_Volterra_with_substrate  
export droop_model
export monod_ierusalimsky
export monod_chemostat
export Synthetic_batch_model
export Synthetic_chemostat_model