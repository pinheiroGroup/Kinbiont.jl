struct Kinbiont_ODEs_system
    name::String
    func::Function
    params::Vector{String}
end

function SIR(du, u, param, t)
    du[1] = - u[1]*u[2] * param[1]
    du[2] = u[1]*u[2] * param[1] - u[2] * param[2]
    du[3] = u[2] * param[2]
end

function SIR_BD(du, u, param, t)
    du[1] = -u[1]*u[2] * param[1] + param[3] * u[1] - u[1] *param[4]
    du[2] = u[1]*u[2] * param[1] - u[2] * param[2]- u[2] *param[4]
    du[3] = u[2] * param[2]- u[3] *param[4]
end

function SIR_B(du, u, param, t)
    du[1] = -u[1]*u[2] * param[1] + param[3] * u[1] 
    du[2] = u[1]*u[2] * param[1] - u[2] * param[2]
    du[3] = u[2] * param[2]- u[3] *param[4]
end

function SIS(du, u, param, t)
    du[1] = -u[1]*u[2] * param[1] + param[2] * u[2] 
    du[2] = u[1]*u[2] * param[1] - u[2] * param[2]
end

function Lotka_Volterra(du, u, param, t)
    du[1] = u[1] * param[1] - param[2] * u[1] * u[2] 
    du[2] = - param[3] * u[2] + param[4]* u[2] * u[1]
end

function Lotka_Volterra_with_substrate(du, u, param, t)
    du[1] = u[1] *  param[1] * (u[3] / (u[3]+param[5])) -  param[2] * u[1] * u[2] 
    du[2] = - param[3] * u[2] + param[4]* u[2] * u[1]
    du[3] = - u[1] *  param[1] * ( u[3] / (u[3] + param[5]))

end

ODEs_system_list = [

    Kinbiont_ODEs_system(
        "SIR",
        SIR,
        ["Infection_rate", "recovery_rate"]
    ),
    Kinbiont_ODEs_system(
        "SIR_BD",
        SIR_BD,
        ["Infection_rate", "recovery_rate","Birth_rate","Death_rate"]
    ),
    Kinbiont_ODEs_system(
        "SIR_B",
        SIR_B,
        ["Infection_rate", "recovery_rate","Birth_rate"]
    ),
    Kinbiont_ODEs_system(
        "SIS",
        SIS,
        ["Infection_rate", "recovery_rate"]
    ),

    Kinbiont_ODEs_system(
        "Lotka_Volterra",
        Lotka_Volterra,
        ["prey_GR", "Prey_DR","Predator_DR","Predator_GR"]
    ),

    Kinbiont_ODEs_system(
        "Lotka_Volterra_with_substrate",
        Lotka_Volterra_with_substrate,
        ["prey_GR", "Prey_DR","Predator_DR","Predator_GR","Mu_max","Monod_constant"]
    ),

   ]


ODEs_system_models = Dict(ODEs_system_list.name => ODEs_system_list for ODEs_system_list in ODEs_system_list)
