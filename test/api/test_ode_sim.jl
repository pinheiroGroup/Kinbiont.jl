# Group 3 — ODE simulation utilities
# Tests for ODE_sim, model_selector, and basic simulation behaviour.

@testset "ODE simulation utilities" begin

    @testset "3.1 ODE_sim logistic — monotone and bounded" begin
        # True params: gr=0.5, N_max=1.5, N0=0.05
        sim = ODE_sim("logistic", [0.05], 0.0, 20.0, 0.5, [0.5, 1.5])
        od  = sim[1, :]

        @test all(diff(od) .>= -1e-6)      # non-decreasing (small ODE solver tolerance)
        @test od[1]   ≈ 0.05 atol=1e-3     # initial condition recovered
        @test od[end] ≈ 1.5  atol=0.05     # converges to N_max
    end

    @testset "3.2 ODE_sim output length matches time grid" begin
        dt  = 0.5
        sim = ODE_sim("logistic", [0.05], 0.0, 10.0, dt, [0.5, 1.5])
        expected = length(0.0:dt:10.0)

        @test length(sim.t) == expected
        @test length(sim.u) == expected
    end

    @testset "3.3 ODE_sim gompertz runs and produces finite values" begin
        sim = ODE_sim("gompertz", [0.05], 0.0, 10.0, 0.5, [0.5, 1.5])

        @test length(sim.t) > 0
        @test all(isfinite.(sim[1, :]))
    end

    @testset "3.4 model_selector returns an ODEProblem" begin
        using SciMLBase: AbstractODEProblem
        prob = model_selector("logistic", [0.05], (0.0, 10.0), [0.5, 1.5])

        @test prob isa AbstractODEProblem
        @test prob.u0 ≈ [0.05]
        @test prob.tspan == (0.0, 10.0)
    end

    @testset "3.5 model_selector covers all registered ODE models" begin
        # Smoke test: every model in the registry can produce an ODEProblem
        for (name, m) in MODEL_REGISTRY
            m isa ODEModel || continue
            u0 = fill(0.05, m.n_eq)
            @test_nowarn model_selector(name, u0, (0.0, 1.0), ones(length(m.param_names)))
        end
    end
end
