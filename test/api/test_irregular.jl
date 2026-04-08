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
