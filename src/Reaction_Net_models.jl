struct Kinbiont_Reaction_network
    name::String
    network::Any
    param::Any
end




model_Michaelis_Menten_enzyme_kinetics= @reaction_network begin
   kB, S + E --> SE
   kD, SE --> S + E
   kP, SE --> P + E
end




Reaction_networks_list = [

    Kinbiont_Reaction_network(
        "Michaelis_Menten",
        model_Michaelis_Menten_enzyme_kinetics,
        parameters(model_Michaelis_Menten_enzyme_kinetics)
    )
    ]


Kinbiont_Reaction_network_models = Dict(Reaction_networks_list.name => Reaction_networks_list for Reaction_networks_list in Reaction_networks_list)



export model_Michaelis_Menten_enzyme_kinetics

