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
        # output length matches input (pipeline pads if needed)
        @test size(processed.curves) == size(data.curves)
    end
end
