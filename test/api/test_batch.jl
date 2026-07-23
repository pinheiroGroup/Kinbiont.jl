using CSV
using DataFrames

@testset "GUI-compatible batch fitting" begin
    Random.seed!(42)
    sim   = ODE_sim("logistic", [0.05], 0.0, 24.0, 1.0, [0.3, 1.2])
    times = collect(0.0:1.0:24.0)
    c1    = sim[1, :]
    c2    = c1 .* 0.9
    # A near-flat curve that should be skipped by skip_flat_threshold
    flat  = fill(0.05, length(times))
    curves = vcat(reshape(c1, 1, :), reshape(c2, 1, :), reshape(flat, 1, :))
    data   = GrowthData(curves, times, ["well_A", "well_B", "flat_well"])

    @testset "kinbiont_batch_fit single model" begin
        batch = kinbiont_batch_fit(
            data;
            experiment="smoke",
            labels=["well_A", "well_B", "flat_well", "missing_well"],
            model_name="NL_logistic",
            optimizer="LN_BOBYQA",
            maxiters=2000,
            skip_flat_threshold=0.05,
        )

        @test batch.experiment == "smoke"
        @test batch.model == "NL_logistic"
        @test batch.model_names == ["NL_logistic"]

        # well_A and well_B should fit
        @test length(batch.results) == 2
        for r in batch.results
            @test r["experiment"] == "smoke"
            @test r["model"] == "NL_logistic"
            @test length(r["parameters"]) == length(r["param_names"])
            @test isfinite(r["aic"])
            @test isfinite(r["loss_rmse"])
            @test r["optimizer_used"] == "LN_BOBYQA"
            @test r["preprocessing"]["smooth"] == false
            @test r["preprocessing"]["cut_stationary_phase"] == true
            @test r["stationary_phase_start"] ≈ last(r["fit_time"])
        end

        # flat_well is below skip threshold
        @test length(batch.skipped) == 1
        @test batch.skipped[1]["well"] == "flat_well"

        # missing_well raises an error string
        @test any(occursin("missing_well", e) for e in batch.errors)
    end

    @testset "kinbiont_batch_fit preprocessing is optimizer-independent" begin
        batch = kinbiont_batch_fit(
            data;
            experiment="optimizer_preprocessing",
            labels=["well_A"],
            model_name="NL_logistic",
            optimizer="LN_COBYLA",
            maxiters=2000,
        )

        @test length(batch.results) == 1
        r = only(batch.results)
        @test r["optimizer_used"] == "LN_COBYLA"
        @test r["preprocessing"]["smooth"] == false
        @test r["preprocessing"]["cut_stationary_phase"] == true
        @test r["stationary_phase_start"] ≈ last(r["fit_time"])
    end

    @testset "kinbiont_batch_fit centered smoothing" begin
        batch = kinbiont_batch_fit(
            data;
            experiment="smoothed",
            labels=["well_A"],
            model_name="NL_logistic",
            optimizer="LN_BOBYQA",
            maxiters=2000,
            smooth=true,
            smooth_window=3,
        )

        @test length(batch.results) == 1
        r = only(batch.results)
        @test r["preprocessing"]["smooth"] == true
        @test r["preprocessing"]["smooth_method"] == "boxcar"
        @test r["preprocessing"]["smooth_window"] == 3
        @test r["smoothed_time"] == times
        @test r["smoothed_od"][2] ≈ mean(c1[1:3])
        @test_throws ArgumentError kinbiont_batch_fit(data; smooth=true, smooth_window=4)
    end

    @testset "kinbiont_batch_fit multi-model" begin
        batch = kinbiont_batch_fit(
            data;
            experiment="multi",
            labels=["well_A"],
            model_names=["NL_logistic", "NL_Gompertz"],
            optimizer="LN_BOBYQA",
            maxiters=2000,
        )

        @test batch.model == "multi"
        @test batch.model_names == ["NL_logistic", "NL_Gompertz"]
        @test length(batch.results) == 1
        @test batch.results[1]["model"] in ("NL_logistic", "NL_Gompertz")
    end

    @testset "kinbiont_batch_fit with compute_loglin companion" begin
        batch = kinbiont_batch_fit(
            data;
            experiment="loglin_companion",
            labels=["well_A"],
            model_name="NL_logistic",
            optimizer="LN_BOBYQA",
            maxiters=2000,
            compute_loglin=true,
        )

        @test length(batch.results) == 1
        r = batch.results[1]
        @test haskey(r, "gr_loglin")
        @test haskey(r, "R_squared_loglin")
        @test haskey(r, "loglin_converged")
    end

    @testset "save_gui_batch_results roundtrip" begin
        batch = kinbiont_batch_fit(
            data;
            experiment="rt",
            labels=["well_A", "well_B"],
            model_name="NL_logistic",
            optimizer="LN_BOBYQA",
            maxiters=2000,
        )
        tmpdir = mktempdir()
        paths  = save_gui_batch_results(batch, tmpdir; prefix="rt_batch")

        @test hasproperty(paths, :summary)
        @test hasproperty(paths, :fitted_curves)
        @test isfile(paths.summary)
        @test isfile(paths.fitted_curves)

        summary = CSV.read(paths.summary, DataFrame)
        @test nrow(summary) == 2
        @test "experiment" in names(summary)
        @test "well" in names(summary)
        @test "model" in names(summary)
        @test "aic" in names(summary)
        @test "loss_rmse" in names(summary)
        @test "optimizer_used" in names(summary)

        curves_df = CSV.read(paths.fitted_curves, DataFrame)
        @test nrow(curves_df) == 2
        @test "experiment" in names(curves_df)
        @test "well" in names(curves_df)
    end

    @testset "kinbiont_batch_loglin" begin
        batch = kinbiont_batch_loglin(
            data;
            experiment="ll",
            labels=["well_A", "well_B", "flat_well"],
            skip_flat_threshold=0.05,
        )

        @test batch.experiment == "ll"
        @test batch.model == "log_lin"
        @test batch.model_names == ["log_lin"]
        @test length(batch.results) == 2
        @test length(batch.skipped) == 1

        for r in batch.results
            @test r["method"] == "Log-lin"
            @test haskey(r, "gr_loglin")
            @test haskey(r, "loglin_converged")
            @test r["loglin_converged"] === true
            @test isfinite(r["gr_loglin"])
            @test isfinite(r["R_squared_loglin"])
        end
    end

    @testset "kinbiont_fit_loglin single-curve API" begin
        one_curve = data[["well_A"]]
        r = kinbiont_fit_loglin(one_curve; experiment="single")

        @test r["well"] == "well_A"
        @test r["loglin_converged"] === true
        @test isfinite(r["gr_loglin"])
        @test isfinite(r["N_max_emp"])

        opts = FitOptions(negative_threshold=0.01)
        r_from_opts = kinbiont_fit_loglin(one_curve, opts; experiment="single")
        @test r_from_opts["N_max_emp"] == r["N_max_emp"]
    end

    @testset "save_gui_batch_loglin_results roundtrip" begin
        batch = kinbiont_batch_loglin(data; experiment="llrt", labels=["well_A"])
        tmpdir = mktempdir()
        paths  = save_gui_batch_loglin_results(batch, tmpdir; prefix="llrt_batch")

        @test hasproperty(paths, :summary)
        @test isfile(paths.summary)

        df = CSV.read(paths.summary, DataFrame)
        @test nrow(df) == 1
        @test "method" in names(df)
        @test "gr_loglin" in names(df)
        @test "R_squared_loglin" in names(df)
    end

    @testset "Empty labels defaults to all wells" begin
        batch = kinbiont_batch_fit(
            data;
            experiment="all",
            model_name="NL_logistic",
            optimizer="LN_BOBYQA",
            maxiters=2000,
            skip_flat_threshold=0.05,
        )
        # well_A and well_B fit, flat_well skipped
        @test length(batch.results) + length(batch.skipped) == 3
    end
end
