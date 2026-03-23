@testset "Preprocessing" begin
    # Small synthetic data: 5 curves × 10 timepoints
    Random.seed!(42)
    times  = collect(0.0:1.0:9.0)
    curves = abs.(rand(5, 10)) .+ 0.1   # all positive
    data   = GrowthData(curves, times, ["c$i" for i in 1:5])

    @testset "No-op (all steps disabled)" begin
        processed = preprocess(data, FitOptions())
        @test processed.curves ≈ data.curves
        @test processed.times  == data.times
        @test processed.clusters === nothing
    end

    @testset "Blank subtraction" begin
        opts = FitOptions(blank_subtraction=true, blank_value=0.05)
        processed = preprocess(data, opts)
        @test processed.curves ≈ data.curves .- 0.05
    end

    @testset "Auto blank detection from labels" begin
        blank_od  = 0.08
        # Build a dataset with two blank wells (label "b") and three signal wells
        n_tp      = length(times)
        blank_row = fill(blank_od, n_tp)
        all_curves = vcat(curves, blank_row', blank_row')
        all_labels = ["c1","c2","c3","c4","c5","b","b"]
        data_with_blanks = GrowthData(all_curves, times, all_labels)

        opts = FitOptions(blank_subtraction=true, blank_from_labels=true)
        processed = preprocess(data_with_blanks, opts)
        # every curve should be shifted by the mean blank (= blank_od)
        @test processed.curves ≈ data_with_blanks.curves .- blank_od
    end

    @testset "Auto blank: warns and uses 0.0 when no 'b' labels" begin
        opts = FitOptions(blank_subtraction=true, blank_from_labels=true)
        @test_warn r"no wells labelled" preprocess(data, opts)
        processed = preprocess(data, opts)
        @test processed.curves ≈ data.curves   # subtracting 0.0
    end
  
    @testset "Replicate averaging collapses duplicate labels" begin
        # 4 curves: A appears twice, B appears twice
        rep_curves = vcat(curves[1:2, :], curves[3:4, :])
        rep_labels = ["A", "A", "B", "B"]
        rep_data   = GrowthData(rep_curves, times, rep_labels)

        opts = FitOptions(average_replicates=true)
        processed = preprocess(rep_data, opts)

        @test size(processed.curves, 1) == 2          # collapsed to 2 unique labels
        @test processed.labels == ["A", "B"]
        @test processed.curves[1, :] ≈ mean(rep_curves[1:2, :], dims=1)[:]
        @test processed.curves[2, :] ≈ mean(rep_curves[3:4, :], dims=1)[:]
    end

    @testset "Replicate averaging drops 'b' and 'X' wells" begin
        rep_curves = vcat(curves[1:2, :], curves[3:4, :], curves[5:5, :])
        rep_labels = ["A", "A", "b", "X", "B"]
        rep_data   = GrowthData(rep_curves, times, rep_labels)

        opts = FitOptions(average_replicates=true)
        processed = preprocess(rep_data, opts)

        @test processed.labels == ["A", "B"]
        @test size(processed.curves, 1) == 2
    end

    @testset "Replicate averaging is a no-op when disabled" begin
        opts = FitOptions(average_replicates=false)
        processed = preprocess(data, opts)
        @test processed.curves ≈ data.curves
        @test processed.labels == data.labels
    end

    @testset "Clustering assigns cluster ids" begin
        opts = FitOptions(cluster=true, n_clusters=2, cluster_trend_test=false)
        processed = preprocess(data, opts)
        @test processed.clusters isa Vector{Int}
        @test length(processed.clusters) == size(data.curves, 1)
    end

    @testset "Pure function — original data not mutated" begin
        original_copy = copy(data.curves)
        _ = preprocess(data, FitOptions(blank_subtraction=true, blank_value=0.1))
        @test data.curves == original_copy
    end

    @testset "Smoothing (rolling_avg) does not error" begin
        opts = FitOptions(smooth=true, smooth_method=:rolling_avg, smooth_pt_avg=3)
        processed = preprocess(data, opts)
        # rolling_avg may shorten the time dimension; number of curves is preserved
        @test size(processed.curves, 1) == size(data.curves, 1)
        @test size(processed.curves, 2) <= size(data.curves, 2)
    end

    @testset "Smoothing (boxcar) preserves shape and time grid" begin
        opts = FitOptions(smooth=true, smooth_method=:boxcar, boxcar_window=3)
        processed = preprocess(data, opts)
        # boxcar is length-preserving: shape and times identical to input
        @test size(processed.curves) == size(data.curves)
        @test processed.times == data.times
        # smoothing must change at least some values
        @test processed.curves != data.curves
    end

    @testset "Smoothing (gaussian) keeps original times when no grid given" begin
        opts = FitOptions(smooth=true, smooth_method=:gaussian, gaussian_h_mult=2.0)
        processed = preprocess(data, opts)
        @test size(processed.curves) == size(data.curves)
        @test processed.times == data.times
    end

    @testset "Gaussian smoothing with custom time grid changes times and columns" begin
        new_times = collect(0.0:0.5:9.0)   # finer grid → 19 points vs original 10
        opts = FitOptions(
            smooth=true,
            smooth_method=:gaussian,
            gaussian_time_grid=new_times,
        )
        processed = preprocess(data, opts)
        @test processed.times == new_times
        @test size(processed.curves, 1) == size(data.curves, 1)   # same n_curves
        @test size(processed.curves, 2) == length(new_times)
    end
end

# 2.3 — negative_value_correction (old-API utility)
@testset "negative_value_correction" begin
    times = collect(0.0:1.0:9.0)
    od    = [-0.05, 0.0, 0.05, 0.1, 0.2, 0.4, 0.8, 1.0, 1.2, 1.3]
    data_mat = Matrix(transpose(hcat(times, od)))   # 2×10 old-API format

    corrected = negative_value_correction(data_mat, Float64[];
                                          method="thr_correction", thr_negative=0.01)
    @test all(corrected[2, :] .>= 0.0)
    @test size(corrected, 1) == 2
end

# 2.4 — specific_gr_evaluation on exponential data
@testset "specific_gr_evaluation" begin
    gr_true = 0.5
    t  = collect(0.0:0.5:10.0)
    od = 0.05 .* exp.(gr_true .* t)
    data_mat = Matrix(transpose(hcat(t, od)))   # 2×N old-API format

    sgr = specific_gr_evaluation(data_mat, 7)

    @test length(sgr) > 0
    @test all(isfinite.(sgr))
    # mean SGR should be close to the true growth rate
    @test abs(mean(sgr) - gr_true) < 0.15
end

# 2.6 — stationary phase cutting (exercised through kinbiont_fit)
@testset "Stationary phase cutting via kinbiont_fit" begin
    Random.seed!(42)
    # Dense logistic curve that clearly reaches a plateau by t=20
    sim   = ODE_sim("logistic", [0.05], 0.0, 30.0, 0.3, [0.5, 1.5])
    times = collect(0.0:0.3:30.0)
    curve = sim[1, :]
    data  = GrowthData(reshape(curve, 1, :), times, ["A"])

    spec      = ModelSpec([MODEL_REGISTRY["NL_logistic"]], [[1.5, 0.3, 0.0]])
    opts_cut  = FitOptions(cut_stationary_phase=true, loss="RE")
    opts_full = FitOptions(cut_stationary_phase=false, loss="RE")

    res_cut  = kinbiont_fit(data, spec, opts_cut)
    res_full = kinbiont_fit(data, spec, opts_full)

    # With cutting, the fitted curve should cover fewer time points
    @test length(res_cut[1].times) < length(res_full[1].times)
end
