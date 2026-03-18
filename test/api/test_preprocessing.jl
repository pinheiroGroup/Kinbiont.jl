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

    @testset "Clustering assigns cluster ids and centroids (z-normalised)" begin
        n_k  = 2
        opts = FitOptions(cluster=true, n_clusters=n_k, cluster_trend_test=false)
        processed = preprocess(data, opts)
        @test processed.clusters  isa Vector{Int}
        @test length(processed.clusters) == size(data.curves, 1)
        # centroid matrix: one row per cluster, one column per timepoint
        @test processed.centroids isa Matrix{Float64}
        @test size(processed.centroids) == (n_k, length(data.times))
        # centroids are in z-normalised space: each non-empty row has ~zero mean
        for k in 1:n_k
            if any(==(k), processed.clusters)
                @test abs(mean(processed.centroids[k, :])) < 0.1
            end
        end
    end

    @testset "Clustering labels always within 1..n_clusters (cluster_trend_test=true)" begin
        n_k = 3
        opts = FitOptions(cluster=true, n_clusters=n_k, cluster_trend_test=true)
        processed = preprocess(data, opts)
        @test processed.clusters isa Vector{Int}
        @test all(1 .<= processed.clusters .<= n_k)
        @test size(processed.centroids) == (n_k, length(data.times))
    end

    @testset "Constant pre-screening keeps labels within 1..n_clusters" begin
        # Mix flat and growing curves so pre-screening has something to detect
        flat_curves = hcat(fill(0.1, 3), fill(0.1, 3), fill(0.1, 3),
                           fill(0.1, 3), fill(0.1, 3), fill(0.1, 3),
                           fill(0.1, 3), fill(0.1, 3), fill(0.1, 3), fill(0.1, 3))
        mixed = vcat(data.curves, flat_curves)
        mixed_data = GrowthData(mixed, data.times, ["c$i" for i in 1:8])
        n_k  = 3
        opts = FitOptions(cluster=true, n_clusters=n_k,
                          cluster_prescreen_constant=true, cluster_trend_test=false)
        processed = preprocess(mixed_data, opts)
        @test all(1 .<= processed.clusters .<= n_k)
        @test size(processed.centroids) == (n_k, length(data.times))
        # flat curves should be assigned to the last label
        @test all(processed.clusters[6:8] .== n_k)
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
