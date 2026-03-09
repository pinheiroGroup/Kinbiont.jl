@testset "MODEL_REGISTRY" begin
    @test MODEL_REGISTRY isa Dict{String, AbstractGrowthModel}

    @test haskey(MODEL_REGISTRY, "NL_logistic")
    @test MODEL_REGISTRY["NL_logistic"] isa NLModel

    @test haskey(MODEL_REGISTRY, "logistic")
    @test MODEL_REGISTRY["logistic"] isa ODEModel

    @testset "Every ODEModel has n_eq >= 1" begin
        for (name, m) in MODEL_REGISTRY
            m isa ODEModel && @test m.n_eq >= 1
        end
    end

    @test MODEL_REGISTRY["NL_logistic"].name == "NL_logistic"

    nl_count  = count(m -> m isa NLModel,  values(MODEL_REGISTRY))
    ode_count = count(m -> m isa ODEModel, values(MODEL_REGISTRY))
    @test nl_count  >= 5
    @test ode_count >= 10
end
