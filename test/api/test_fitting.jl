@testset "kinbiont_fit" begin
    # Shared synthetic data: 2 logistic curves
    Random.seed!(42)
    sim   = ODE_sim("logistic", [0.05], 0.0, 24.0, 1.0, [0.3, 1.2])
    times = collect(0.0:1.0:24.0)
    c1    = sim[1, :]
    c2    = c1 .* 0.9   # slight variant for a second curve
    curves = vcat(reshape(c1, 1, :), reshape(c2, 1, :))
    data   = GrowthData(curves, times, ["curve1", "curve2"])

    @testset "Single NL model" begin
        spec = ModelSpec([MODEL_REGISTRY["NL_logistic"]], [[1.5, 0.3, 0.0]])
        res  = kinbiont_fit(data, spec, FitOptions(loss="RE"))

        @test res isa GrowthFitResults
        @test length(res) == 2
        @test res[1].best_model isa NLModel
        @test length(res[1].best_params) == 3
        @test res[1].best_aic isa Float64
        @test isfinite(res[1].best_aic)
        @test length(res[1].fitted_curve) <= length(times)
        @test res[1].loss isa Float64
    end

    @testset "Multi-model NL selection" begin
        spec = ModelSpec(
            [MODEL_REGISTRY["NL_logistic"], MODEL_REGISTRY["NL_Gompertz"]],
            [[1.5, 0.3, 0.0], [1.5, 0.3, 20.0]],
        )
        res = kinbiont_fit(data, spec, FitOptions())

        @test length(res[1].all_results) == 2

        aics = [r.aic for r in res[1].all_results]
        best_idx = argmin(aics)
        @test res[1].best_aic == aics[best_idx]
    end

    @testset "ODE model" begin
        spec = ModelSpec([MODEL_REGISTRY["logistic"]], [[0.3, 1.2]])
        res  = kinbiont_fit(data, spec)

        @test res[1].best_model isa ODEModel
        @test length(res[1].best_params) >= 1
        @test res[1].loss >= 0
    end

    @testset "LogLin model" begin
        spec = ModelSpec([LogLinModel()], [Float64[]])
        res  = kinbiont_fit(data, spec)

        @test res[1].best_model isa LogLinModel
        @test res[1].fitted_curve isa Vector{Float64}
        @test length(res[1].fitted_curve) <= length(times)
    end

    @testset "DDDE model" begin
        spec = ModelSpec([DDDEModel()], [Float64[]])
        res  = kinbiont_fit(data, spec)

        @test res[1].best_model isa DDDEModel
        @test res[1].fitted_curve isa Vector{Float64}
        @test length(res[1].fitted_curve) == length(times)
        @test res[1].loss isa Float64
        @test res[1].loss >= 0
        @test isfinite(res[1].best_aic)
        # params are the discovered polynomial coefficients
        @test length(res[1].best_params) == DDDEModel().max_degree + 1
    end
end
