using CSV
using DataFrames

@testset "save_results I/O" begin
    Random.seed!(42)
    sim   = ODE_sim("logistic", [0.05], 0.0, 24.0, 1.0, [0.3, 1.2])
    times = collect(0.0:1.0:24.0)
    c1    = sim[1, :]
    c2    = c1 .* 0.9
    curves = vcat(reshape(c1, 1, :), reshape(c2, 1, :))
    data   = GrowthData(curves, times, ["well_A", "well_B"])

    spec = ModelSpec([MODEL_REGISTRY["NL_logistic"]], [[1.5, 0.3, 0.0]])
    res  = kinbiont_fit(data, spec, FitOptions(loss="RE"))

    tmpdir = mktempdir()
    paths  = save_results(res, tmpdir)

    @testset "File creation" begin
        @test hasproperty(paths, :summary)
        @test hasproperty(paths, :fitted_curves)
        @test hasproperty(paths, :all_models)
        @test isfile(paths.summary)
        @test isfile(paths.fitted_curves)
        @test isfile(paths.all_models)
    end

    @testset "Summary CSV structure" begin
        df = CSV.read(paths.summary, DataFrame)
        @test nrow(df) == 2
        @test "label"      in names(df)
        @test "best_model" in names(df)
        @test "n_params"   in names(df)
        @test "aic"        in names(df)
        @test "loss"       in names(df)
        @test "param_1"    in names(df)
    end

    @testset "Fitted curves CSV structure" begin
        df = CSV.read(paths.fitted_curves, DataFrame)
        @test "label"    in names(df)
        @test "time"     in names(df)
        @test "observed" in names(df)
        @test "fitted"   in names(df)
        @test nrow(df) >= 2 * length(times)
        @test Set(df.label) == Set(["well_A", "well_B"])
    end

    @testset "All-models CSV structure" begin
        df = CSV.read(paths.all_models, DataFrame)
        @test "label"      in names(df)
        @test "model_name" in names(df)
        @test "aic"        in names(df)
        @test "loss"       in names(df)
        @test "is_best"    in names(df)
        # one model, two curves → 2 rows
        @test nrow(df) == 2 * length(spec.models)
        # exactly one is_best per label
        for lbl in ["well_A", "well_B"]
            subset_df = filter(row -> row.label == lbl, df)
            @test sum(subset_df.is_best) == 1
        end
    end
end
