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
