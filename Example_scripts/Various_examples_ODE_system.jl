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
    k1, k2, k3, k4, k_cat, n, e0 = param
    
    du[1] = k4 * x - k3 * m * e + k2 * y^n - k1 * e + k_cat * x  # Free active enzymes
    du[2] = k3 * m * e - k4 * x - k_cat * x                      # Enzyme-substrate complex
    du[3] = k1 * e - k2 * y^n                                    # Inactive aggregates
    du[4] = -du[1]                                              # Substrate degradation rate
end

u0 = [1.0, 0.1, 0.1, 1.0]  # Initial conditions [e, x, y, m]
param = [0.1, 0.1, 0.05, 0.05, 0.02, 2, 1.0]  # [k1, k2, k3, k4, k_cat, n, e0]


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

# Examples of structured model
function structured_model(du, u, param, t)
    s, x, m1, m2, m3, m4, p1, p2 = u
    qs, Y, v1, v2, v3, v4, v5 = param
    
    du[1] = -qs * u[2] * u[1]                                 # Extracellular substrate
    du[2] = Y * qs * u[2] * u[1]                                         # Cell mass
    du[3] = qs * x - v1 * x                                    # Intracellular metabolite M1
    du[4] = v1 * x - v2 * x - v4 * x                           # Intracellular metabolite M2
    du[5] = v2 * x - v3 * x                                    # Intracellular metabolite M3
    du[6] = v4 * x - v5 * x                                    # Intracellular metabolite M4
    du[7] = v3 * x                                             # Secreted product P1
    du[8] = v5 * x                                             # Secreted product P2
end

u0_structured = [1.0, 0.5, 0.1, 0.1, 0.1, 0.1, 0.0, 0.0]  # Initial conditions for metabolism
param_structured = [0.1, 0.5, 0.05, 0.05, 0.02, 0.02, 0.01]  # Parameters for metabolism



Simulation =  ODEs_system_sim(
    structured_model, #string of the model
    u0_structured, # starting condition
    0.0, # start time of the sim
    50.0, # final time of the sim
    1.0, # delta t for poisson approx
    param_structured; # parameters of the ODE model
)

# Plot the data


scatter(Simulation)

############################# From Baldazzi et al. eLife 2023;12:e79815. DOI: https://doi.org/10.7554/eLife.79815 Resource allocation accounts for the  large variability of rate- yield phenotypes    across bacterial strains

function ode_system_1(du, u, p, t)
    # Unpack variables
    C, U, Mu, R, Mc, Mer, Mef = u
    Vmc, Vmer, ρmef, Vmef, ρru, Vr, Vmu, γ, χu, χr, χc, χer, χef = p

    # Define the ODEs
    du[1] = Vmc - Vmer - ρmef * Vmef - ρru * (Vr + Vmu) - γ * C
    du[2] = Vmu - γ * U
    du[3] = χu * Vr - γ * Mu
    du[4] = χr * Vr - γ * R
    du[5] = χc * Vr - γ * Mc
    du[6] = χer * Vr - γ * Mer
    du[7] = χef * Vr - γ * Mef
end


u0 = [1.0, 0.5, 0.1, 0.1, 0.1, 0.1, 0.0]  # Initial conditions for metabolism
param = [0.1, 0.5, 0.05, 0.05, 0.02, 0.02, 0.01,0.1, 0.5, 0.05, 0.05, 0.02, 0.02, 0.0,0.1]  # Parameters for metabolism



Simulation =  ODEs_system_sim(
    ode_system_1, #string of the model
    u0, # starting condition
    0.0, # start time of the sim
    50.0, # final time of the sim
    1.0, # delta t for poisson approx
    param; # parameters of the ODE model
)
scatter(Simulation)




function ode_system2(du, u, p, t)
    # Unpack variables
    c, u, mu, r, mc, mer, mef, a_star, a = u
    vmc, vmer, ρmef, vmef, ρru, vr, vmu, μ, γ, χu, χr, χc, χer, χef, nmer, nmef, nr, nmu, vd = p

    # Define the ODEs
    du[1] = vmc - vmer - ρmef * vmef - ρru * (vr + vmu) - (μ + γ) * c
    du[2] = vmu - (μ + γ) * u
    du[3] = χu * vr - (μ + γ) * mu
    du[4] = χr * vr - (μ + γ) * r
    du[5] = χc * vr - (μ + γ) * mc
    du[6] = χer * vr - (μ + γ) * mer
    du[7] = χef * vr - (μ + γ) * mef
    du[8] = nmer * vmer + nmef * vmef - nr * vr - nmu * vmu - vd
    du[9] = -nmer * vmer - nmef * vmef + nr * vr + nmu * vmu + vd
end

# Initial conditions
u0 = [1.0, 0.5, 0.1, 0.1, 0.1, 0.1, 0.0, 0.0, 0.0]  # Initial conditions for metabolism and energy

# Parameters
param = [0.1, 0.5, 0.05, 0.05, 0.02, 0.02, 0.01, 0.1, 0.1, 0.5, 0.05, 0.05, 0.02, 0.02, 0.0, 0.1, 0.1, 0.1, 0.1, 0.1]


Simulation =  ODEs_system_sim(
    ode_system2, #string of the model
    u0, # starting condition
    0.0, # start time of the sim
    50.0, # final time of the sim
    1.0, # delta t for poisson approx
    param; # parameters of the ODE model
)
scatter(Simulation)



function vd(a_star, kd)
    return kd * a_star
end

function vmc(mc)
    return 49.6 # Placeholder, replace with actual function if needed
end

function vmer(mer, c, a)
    return 4.6 # Placeholder, replace with actual function if needed
end

function vmef(mef, c, a)
    return 9.8 # Placeholder, replace with actual function if needed
end

function vr(r, c, a_star)
    return 19.2 # Placeholder, replace with actual function if needed
end

function vmu(mu, c, a_star)
    return 6.5 # Placeholder, replace with actual function if needed
end

function ode_system4(du, u, p, t)
    # Unpack variables
    c, u, mu, r, mc, mer, mef, a_star, a = u
    kd, ρmef, ρru, γ, χu, χr, χc, χer, χef, nmer, nmef, nr, nmu, β = p

    # Calculate μ
    μ = β * (vmc(mc) - vmer(mer, c, a) - ρmef * vmef(mef, c, a) - (ρru - 1) * (vr(r, c, a_star) + vmu(mu, c, a_star))) - γ

    # Define the ODEs
    du[1] = vmc(mc) - vmer(mer, c, a) - ρmef * vmef(mef, c, a) - ρru * (vr(r, c, a_star) + vmu(mu, c, a_star)) - (μ + γ) * c
    du[2] = vmu(mu, c, a_star) - (μ + γ) * u
    du[3] = χu * vr(r, c, a_star) - (μ + γ) * mu
    du[4] = χr * vr(r, c, a_star) - (μ + γ) * r
    du[5] = χc * vr(r, c, a_star) - (μ + γ) * mc
    du[6] = χer * vr(r, c, a_star) - (μ + γ) * mer
    du[7] = χef * vr(r, c, a_star) - (μ + γ) * mef
    du[8] = nmer * vmer(mer, c, a) + nmef * vmef(mef, c, a) - nr * vr(r, c, a_star) - nmu * vmu(mu, c, a_star) - vd(a_star, kd)
    du[9] = -nmer * vmer(mer, c, a) - nmef * vmef(mef, c, a) + nr * vr(r, c, a_star) + nmu * vmu(mu, c, a_star) + vd(a_star, kd)
end

# Initial conditions
u0 = [0.35, 10.2, 11.1, 13.2, 2.7, 1.9, 1.1, 0.009, 0.011]  # Initial conditions for metabolism and energy

# Parameters
param = [0.1, 0.05, 0.02, 0.027, 0.5, 0.05, 0.05, 0.02, 0.02, 0.1, 0.1, 0.1, 0.1, 40.65]

Simulation =  ODEs_system_sim(
    ode_system4, #string of the model
    u0, # starting condition
    0.0, # start time of the sim
    50.0, # final time of the sim
    1.0, # delta t for poisson approx
    param; # parameters of the ODE model
)
scatter(Simulation)


### from  Abacterial size law revealed by a coarse grained model of cell physiolog



function ode_system5(du, u, p, t)
    # Unpack variables
    A, E, Ra, Q, X, U, Ri = u
    kE, fE, fR, kcm_on, kcm_off, fQ, fX, fU, s, asat = p

    # Define the ODEs
    du[1] = kE - s*Ra * A / (A + asat)
    du[2] = fE * s*Ra * A / (A + asat)
    du[3] = fR * s*Ra * A / (A + asat) - kcm_on * Ra + kcm_off * Ri
    du[4] = fQ * s*Ra * A / (A + asat)
    du[5] = fX * s*Ra * A / (A + asat)
    du[6] = fU * s*Ra * A / (A + asat)
    du[7] = kcm_on * Ra - kcm_off * Ri
end

# Initial conditions
u0 = [1.0, 0.5, 0.1, 0.1, 0.1, 0.1, 0.0]  # Initial conditions for the variables

# Parameters
param = [0.1, 0.5, 0.05, 0.05, 0.02, 0.02, 0.01, 0.1, 0.5, 0.05]



Simulation =  ODEs_system_sim(
    ode_system5, #string of the model
    u0, # starting condition
    0.0, # start time of the sim
    50.0, # final time of the sim
    1.0, # delta t for poisson approx
    param; # parameters of the ODE model
)
scatter(Simulation)

##########################################################



function Central_dogma(du, u, param, t)
    # Unpack variables
    DNA, RNA, P, nutrients = u
    yield, ktr, ktl, kdeg_RNA, kdeg_P  = param


    # Define the ODEs
    # du[1] ODE DNA replication using nutrients
    du[1] = + yield * DNA * nutrients
    # du[2] ODE rna dynamics
    du[2] = ktr * DNA - kdeg_RNA * RNA
    # du[3] protein translation
    du[3] = ktl * RNA - kdeg_P * P
    # nutrient consumption
    du[4] = -yield * DNA * nutrients

end

u0 = [1.0, 0.0, 0.0,  1.1]  # Initial conditions for the variables

# Parameters
param = [0.1, 0.05, 0.05, 0.05, 0.002, 0.002]



Simulation =  ODEs_system_sim(
    Central_dogma, #string of the model
    u0, # starting condition
    0.0, # start time of the sim
    50.0, # final time of the sim
    1.0, # delta t for poisson approx
    param; # parameters of the ODE model
)
scatter(Simulation)

# concentration version 



##########################################################



function Central_dogma_with_ribo(du, u, param, t)
    # Unpack variables
    DNA, RNA, P, R, nutrients = u
    yield, ktr, ktl, kdeg_RNA, kdeg_P , R_size = param


    # Define the ODEs
    # du[1] ODE DNA replication using nutrients
    du[1] = + yield * DNA * nutrients
    # du[2] ODE rna dynamics
    du[2] = ktr * DNA - kdeg_RNA * RNA
    # du[3] protein translation
    du[3] = (1-R_size)  * ktl * R* RNA - (1-R_size)  * kdeg_P * P
    # du[4] ribosome dynamics   
    du[4] = R_size * ktl * R * RNA - R * kdeg_P
    # nutrient consumption
    du[5] = -yield * DNA * nutrients

end



u0 = [0.2, 0.0, 0.0,  0.1,1.1]  # Initial conditions for the variables

# Parameters
param = [0.1, 0.05, 0.05, 0.002, 0.002, 0.6]



Simulation =  ODEs_system_sim(
    Central_dogma_with_ribo, #string of the model
    u0, # starting condition
    0.0, # start time of the sim
    60.0, # final time of the sim
    1.0, # delta t for poisson approx
    param; # parameters of the ODE model
)
scatter(Simulation)