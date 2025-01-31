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






# Define the Synthetic Chemostat Model ODE system with P and U components
function synthetic_chemostat_model(du, u, param, t)
    # Extract state variables
    x = u[1]      # Biomass concentration (x)
    s = u[2]      # Substrate concentration (s)
    r = u[3]      # biological inertia

    # Extract parameters
    D = param[1]     # Dilution rate (D)
    s_r = param[2]   # Inflow substrate concentration (s_r)
    Q_s = param[3]   # Q_s (nutrient uptake rate coefficient)
    Q_s_prime = param[4]  # Q_s' (second nutrient uptake coefficient)
    K_s = param[5]   # K_s (half-saturation constant)
    K_s_prime = param[6]  # K_s' (half-saturation constant)
    Y = param[7]     # Yield coefficient (Y)
    a_0 = param[8]   # Biological inertia constant (a_0)
    K_r = param[9]   # Half-saturation constant for r (K_r)

    # Calculate q_s (nutrient uptake rate)
    q_s = r * Q_s * K_s / (K_s + s) + (1 - r) * Q_s_prime * K_s_prime / (K_s_prime + s)

    # Define the ODEs
    du[1] = Y * q_s - a_0 * r * x - D * x        # d(x)/dt (biomass dynamics)
    du[2] = D * (s_r - s) - q_s * x               # d(s)/dt (substrate dynamics)
    du[3] = (Y * q_s - a_0 * r) * (s / (K_r + s) - r)  # d(r)/dt (r dynamics)

    return du
end

# Monod Chemostat Model


ODEs_system_list = [
    Kinbiont_ODEs_system(
        "Monod_Chemostat",  # Model name
        monod_chemostat,    # The ODE function
        ["Ks", "m", "Ymax", "mum", "D", "sr"]  # Parameters list
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


    # Synthetic Chemostat Model
    Kinbiont_ODEs_system(
        "Synthetic_Chemostat",  # Model name
        synthetic_chemostat_model,  # The ODE function
        ["D", "s_r", "Q_s", "Q_s_prime", "K_s", "K_s_prime", "Y", "a_0", "K_r"]  # Parameters list
    ), Kinbiont_ODEs_system(
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
export synthetic_chemostat_model
export droop_model
export monod_ierusalimsky
export monod_chemostat