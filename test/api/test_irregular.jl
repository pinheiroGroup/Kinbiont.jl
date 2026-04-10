# test/api/test_irregular.jl
# Tests are added incrementally — one @testset per task below.

@testset "_normalize_times_01" begin
    # normal case: maps to [0, 1]
    t = [0.0, 5.0, 10.0, 20.0]
    τ = Kinbiont._normalize_times_01(t)
    @test τ[1]   ≈ 0.0
    @test τ[end] ≈ 1.0
    @test τ[2]   ≈ 0.25
    @test τ[3]   ≈ 0.5

    # irregular spacing: result still in [0, 1]
    t2 = [3.0, 3.5, 7.0, 9.0]
    τ2 = Kinbiont._normalize_times_01(t2)
    @test τ2[1] ≈ 0.0
    @test τ2[end] ≈ 1.0
    @test all(0.0 .<= τ2 .<= 1.0)

    # zero-span (constant time): returns all zeros, no error
    t3 = [5.0, 5.0, 5.0]
    τ3 = Kinbiont._normalize_times_01(t3)
    @test all(τ3 .== 0.0)
    @test length(τ3) == 3
end

@testset "_build_union_grid" begin
    # always includes 0.0 and 1.0
    t01_list = [[0.0, 0.33, 0.66, 1.0], [0.0, 0.5, 1.0]]
    g = Kinbiont._build_union_grid(t01_list; step=0.01)
    @test g[1]   == 0.0
    @test g[end] == 1.0
    @test issorted(g)
    @test allunique(g)

    # step respected: all entries are multiples of step
    @test all(v -> isapprox(v / 0.01, round(v / 0.01); atol=1e-10), g)

    # single-curve list with step=0.25: sources snap to 0.0, 0.25, 0.75, 1.0
    g2 = Kinbiont._build_union_grid([[0.0, 0.25, 0.75, 1.0]]; step=0.25)
    @test g2[1] == 0.0
    @test g2[end] == 1.0
    @test length(g2) == 4

    # with many curves covering dense range, all 101 bins appear
    dense = [[i * 0.01 for i in 0:100] for _ in 1:3]
    g3 = Kinbiont._build_union_grid(dense; step=0.01)
    @test length(g3) == 101

    # endpoint guard fires when source times don't snap to 0.0 or 1.0
    g4 = Kinbiont._build_union_grid([[0.15, 0.45, 0.85]]; step=0.1)
    @test g4[1]   == 0.0
    @test g4[end] == 1.0
end

@testset "_interp_linear" begin
    x = [0.0, 0.5, 1.0]
    y = [0.0, 1.0, 0.0]   # triangle

    # exact at known points
    @test Kinbiont._interp_linear(x, y, [0.0])[1] ≈ 0.0
    @test Kinbiont._interp_linear(x, y, [0.5])[1] ≈ 1.0
    @test Kinbiont._interp_linear(x, y, [1.0])[1] ≈ 0.0

    # midpoint interpolation
    @test Kinbiont._interp_linear(x, y, [0.25])[1] ≈ 0.5
    @test Kinbiont._interp_linear(x, y, [0.75])[1] ≈ 0.5

    # clamped below lower bound
    @test Kinbiont._interp_linear(x, y, [-0.1])[1] ≈ 0.0   # clamped to y[1]

    # clamped above upper bound
    @test Kinbiont._interp_linear(x, y, [1.5])[1] ≈ 0.0    # clamped to y[end]

    # multiple query points at once
    result = Kinbiont._interp_linear(x, y, [0.0, 0.25, 0.5, 0.75, 1.0])
    @test result ≈ [0.0, 0.5, 1.0, 0.5, 0.0]

    # single-interval source
    x2 = [0.0, 1.0]; y2 = [2.0, 4.0]
    @test Kinbiont._interp_linear(x2, y2, [0.5])[1] ≈ 3.0
end

@testset "_resample_to_union_grid" begin
    # Two curves sampled at different (already-normalized) time points
    # Curve 1: linear ramp on [0.0, 0.5, 1.0]
    # Curve 2: flat 0.5 on [0.0, 1.0]
    t1 = [0.0, 0.5, 1.0];  y1 = [0.0, 0.5, 1.0]
    t2 = [0.0, 1.0];        y2 = [0.5, 0.5]
    grid = [0.0, 0.25, 0.5, 0.75, 1.0]

    X = Kinbiont._resample_to_union_grid([t1, t2], [y1, y2], grid; step=0.25)
    @test size(X) == (2, 5)

    # curve 1 is a ramp — check sampled values
    @test X[1, :] ≈ [0.0, 0.25, 0.5, 0.75, 1.0]

    # curve 2 is flat — all values equal 0.5
    @test all(X[2, :] .≈ 0.5)

    # output matrix has no NaN or Inf
    @test all(isfinite.(X))
end

@testset "IrregularGrowthData constructor" begin
    # Three curves with irregular times (minutes, not normalized)
    rc = [[0.1, 0.2, 0.5, 0.9], [0.05, 0.15, 0.4, 0.8], [0.08, 0.18, 0.45, 0.85]]
    rt = [[0.0, 15.0, 45.0, 90.0],
          [5.0, 20.0, 50.0, 95.0],
          [0.0, 30.0, 60.0, 120.0]]
    labels = ["well_A", "well_B", "well_C"]

    data = IrregularGrowthData(rc, rt, labels; step=0.01)

    # raw fields preserved unchanged
    @test data.raw_curves === rc
    @test data.raw_times  === rt
    @test data.labels     === labels

    # resampled matrix has correct shape
    @test size(data.curves, 1) == 3
    @test data.times[1]   == 0.0
    @test data.times[end] == 1.0
    @test issorted(data.times)

    # clustering fields are nothing until preprocess is called
    @test data.clusters  === nothing
    @test data.centroids === nothing
    @test data.wcss      === nothing

    # resampled values are finite
    @test all(isfinite.(data.curves))
end

@testset "IrregularGrowthData constructor validation" begin
    good_rc = [[0.1, 0.5, 0.9], [0.1, 0.5, 0.9]]
    good_rt = [[0.0, 5.0, 10.0], [0.0, 5.0, 10.0]]
    good_lb = ["A", "B"]

    # raw_curves / raw_times length mismatch
    @test_throws ArgumentError IrregularGrowthData(good_rc, [[0.0, 5.0, 10.0]], good_lb)

    # labels length mismatch
    @test_throws ArgumentError IrregularGrowthData(good_rc, good_rt, ["A"])

    # curve/time length mismatch within a pair
    @test_throws ArgumentError IrregularGrowthData(
        [[0.1, 0.5], [0.1, 0.5, 0.9]],
        [[0.0, 5.0, 10.0], [0.0, 5.0, 10.0]],
        good_lb,
    )

    # time vector too short (< 2 points)
    @test_throws ArgumentError IrregularGrowthData(
        [[0.1], [0.1, 0.5]],
        [[0.0], [0.0, 5.0]],
        good_lb,
    )
end

@testset "IrregularGrowthData getindex subsetting" begin
    rc = [[0.1, 0.2, 0.5], [0.2, 0.4, 0.7], [0.05, 0.3, 0.6]]
    rt = [[0.0, 10.0, 20.0], [5.0, 15.0, 25.0], [0.0, 8.0, 16.0]]
    lb = ["A", "B", "C"]
    data = IrregularGrowthData(rc, rt, lb)

    # single label
    sub = data[["B"]]
    @test sub.labels == ["B"]
    @test sub.raw_curves == [rc[2]]
    @test sub.raw_times  == [rt[2]]
    @test size(sub.curves, 1) == 1
    @test sub.times == data.times       # same union grid
    @test sub.clusters  === nothing
    @test sub.centroids === nothing
    @test sub.wcss      === nothing

    # multiple labels in different order
    sub2 = data[["C", "A"]]
    @test sub2.labels     == ["C", "A"]
    @test sub2.raw_curves == [rc[3], rc[1]]

    # unknown label throws ArgumentError
    @test_throws ArgumentError data[["A", "Z"]]

    # clustering fields are dropped even when source has them
    data_cl = IrregularGrowthData(rc, rt, lb, data.curves, data.times, [1,2,3])
    sub3 = data_cl[["A"]]
    @test sub3.clusters === nothing
end

@testset "preprocess(IrregularGrowthData) clustering" begin
    Random.seed!(42)
    # 6 curves, 2 clear shape groups: fast-rising (logistic, midpoint=0.2) and
    # slow-rising (logistic, midpoint=0.8), 3 curves each, on irregular times
    function logistic_sample(tau, mid)
        0.05 .+ 0.9 ./ (1 .+ exp.(-12 .* (tau .- mid)))
    end

    raw_curves = Vector{Vector{Float64}}()
    raw_times  = Vector{Vector{Float64}}()
    for i in 1:3
        t  = sort(rand(20)) .* 100.0   # random times in [0, 100]
        τ  = (t .- minimum(t)) ./ (maximum(t) - minimum(t))
        push!(raw_times, t)
        push!(raw_curves, logistic_sample(τ, 0.2) .+ 0.005 .* randn(20))
    end
    for i in 1:3
        t  = sort(rand(20)) .* 100.0
        τ  = (t .- minimum(t)) ./ (maximum(t) - minimum(t))
        push!(raw_times, t)
        push!(raw_curves, logistic_sample(τ, 0.8) .+ 0.005 .* randn(20))
    end
    labels = vcat(["early_$i" for i in 1:3], ["late_$i" for i in 1:3])
    data   = IrregularGrowthData(raw_curves, raw_times, labels)

    opts = FitOptions(cluster=true, n_clusters=2,
                      cluster_trend_test=false, kmeans_seed=42, kmeans_n_init=10)
    result = preprocess(data, opts)

    # type preserved
    @test result isa IrregularGrowthData

    # raw fields unchanged
    @test result.raw_curves === data.raw_curves
    @test result.raw_times  === data.raw_times
    @test result.labels     === data.labels
    @test result.curves     === data.curves
    @test result.times      === data.times

    # clustering fields populated
    @test result.clusters  isa Vector{Int}
    @test length(result.clusters) == 6
    @test all(1 .<= result.clusters .<= 2)
    @test result.centroids isa Matrix{Float64}
    @test size(result.centroids) == (2, length(data.times))
    @test result.wcss isa Float64
    @test result.wcss >= 0.0

    # the two groups should be in different clusters
    early_cluster = result.clusters[1]
    late_cluster  = result.clusters[4]
    @test early_cluster != late_cluster
    @test all(result.clusters[1:3] .== early_cluster)
    @test all(result.clusters[4:6] .== late_cluster)
end

@testset "preprocess(IrregularGrowthData) unsupported opts warnings" begin
    rc = [[0.1, 0.5, 0.9], [0.1, 0.5, 0.9]]
    rt = [[0.0, 5.0, 10.0], [0.0, 5.0, 10.0]]
    data = IrregularGrowthData(rc, rt, ["A", "B"])

    @test_warn r"blank_subtraction" preprocess(data, FitOptions(blank_subtraction=true))
    @test_warn r"correct_negatives" preprocess(data, FitOptions(correct_negatives=true))
    @test_warn r"smooth"            preprocess(data, FitOptions(smooth=true))
end

@testset "preprocess(IrregularGrowthData) pure function" begin
    rc = [[0.1, 0.5, 0.9], [0.2, 0.6, 1.0]]
    rt = [[0.0, 5.0, 10.0], [0.0, 5.0, 10.0]]
    data = IrregularGrowthData(rc, rt, ["A", "B"])
    original_curves = copy(data.curves)
    _ = preprocess(data, FitOptions(cluster=true, n_clusters=2, cluster_trend_test=false))
    @test data.curves == original_curves
    @test data.clusters === nothing
end

@testset "preprocess(IrregularGrowthData) cluster=false no-op" begin
    rc = [[0.1, 0.5, 0.9], [0.2, 0.6, 1.0]]
    rt = [[0.0, 5.0, 10.0], [0.0, 5.0, 10.0]]
    data = IrregularGrowthData(rc, rt, ["A", "B"])
    result = preprocess(data, FitOptions())   # cluster=false by default
    @test result isa IrregularGrowthData
    @test result.clusters  === nothing
    @test result.centroids === nothing
    @test result.wcss      === nothing
end
