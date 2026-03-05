@testset "GrowthData constructor" begin
    curves = rand(3, 10)
    times  = collect(1.0:10.0)
    labels = ["a", "b", "c"]

    gd = GrowthData(curves, times, labels)
    @test gd.curves === curves
    @test gd.times  === times
    @test gd.labels === labels
    @test gd.clusters === nothing

    # times length mismatch
    @test_throws ErrorException GrowthData(curves, collect(1.0:9.0), labels)

    # labels length mismatch
    @test_throws ErrorException GrowthData(curves, times, ["a", "b"])

    # clusters length mismatch
    @test_throws ErrorException GrowthData(curves, times, labels, [1, 2])
end

@testset "FitOptions defaults" begin
    opts = FitOptions()
    @test opts.loss   == "RE"
    @test opts.smooth == false
    @test opts.cluster == false
end

@testset "FitOptions keyword override" begin
    opts = FitOptions(loss="L2")
    @test opts.loss == "L2"
    @test opts.smooth == false
end

@testset "ModelSpec validation" begin
    m1 = MODEL_REGISTRY["NL_logistic"]
    m2 = MODEL_REGISTRY["logistic"]

    # valid construction
    spec = ModelSpec([m1], [[1.0, 0.3, 0.0]])
    @test length(spec.models) == 1

    # models/params length mismatch
    @test_throws ErrorException ModelSpec([m1, m2], [[1.0]])

    # lower bounds length mismatch
    @test_throws ErrorException ModelSpec([m1], [[1.0, 0.3, 0.0]]; lower=[[0.1], [0.1]])

    # upper bounds length mismatch
    @test_throws ErrorException ModelSpec([m1], [[1.0, 0.3, 0.0]]; upper=[[10.0], [10.0]])
end
