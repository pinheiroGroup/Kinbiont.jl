@testset "GrowthData constructor" begin
    curves = rand(3, 10)
    times  = collect(1.0:10.0)
    labels = ["a", "b", "c"]

    gd = GrowthData(curves, times, labels)
    @test gd.curves    === curves
    @test gd.times     === times
    @test gd.labels    === labels
    @test gd.clusters  === nothing
    @test gd.centroids === nothing

    # times length mismatch
    @test_throws ErrorException GrowthData(curves, collect(1.0:9.0), labels)

    # labels length mismatch
    @test_throws ErrorException GrowthData(curves, times, ["a", "b"])

    # clusters length mismatch
    @test_throws ErrorException GrowthData(curves, times, labels, [1, 2])
end

@testset "GrowthData from CSV path" begin
    # Write a minimal CSV in the Kinbiont column format:
    # first column = time, remaining columns = curves (header = label)
    tmp = tempname() * ".csv"
    open(tmp, "w") do io
        println(io, "time,well_A,well_B")
        println(io, "0.0,0.05,0.06")
        println(io, "1.0,0.10,0.12")
        println(io, "2.0,0.20,0.22")
    end

    gd = GrowthData(tmp)
    @test gd.times  == [0.0, 1.0, 2.0]
    @test gd.labels == ["well_A", "well_B"]
    @test size(gd.curves) == (2, 3)          # 2 curves × 3 timepoints
    @test gd.curves[1, :] ≈ [0.05, 0.10, 0.20]
    @test gd.curves[2, :] ≈ [0.06, 0.12, 0.22]
    @test gd.clusters === nothing
end

@testset "FitOptions defaults" begin
    opts = FitOptions()
    @test opts.loss       == "RE"
    @test opts.smooth     == false
    @test opts.cluster    == false
    @test opts.opt_params == (;)
end

@testset "FitOptions keyword override" begin
    opts = FitOptions(loss="L2")
    @test opts.loss == "L2"
    @test opts.smooth == false
end

@testset "FitOptions opt_params passthrough" begin
    opts = FitOptions(opt_params=(maxiters=100_000, abstol=1e-8))
    @test opts.opt_params.maxiters == 100_000
    @test opts.opt_params.abstol   == 1e-8
end

@testset "DDDEModel construction" begin
    m = DDDEModel()
    @test m isa AbstractGrowthModel
    @test m.max_degree  == 4
    @test m.lambda_min  == -5.0
    @test m.lambda_max  == -1.0
    @test m.lambda_step == 0.5

    m2 = DDDEModel(max_degree=3, lambda_min=-4.0, lambda_max=-2.0, lambda_step=1.0)
    @test m2.max_degree  == 3
    @test m2.lambda_min  == -4.0
    @test m2.lambda_max  == -2.0
    @test m2.lambda_step == 1.0
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

@testset "GrowthData subsetting via getindex" begin
    curves = Float64[1 2 3; 4 5 6; 7 8 9]   # 3 curves × 3 timepoints
    times  = [0.0, 1.0, 2.0]
    labels = ["A1", "B2", "C3"]
    gd = GrowthData(curves, times, labels)

    # single well
    sub = gd[["B2"]]
    @test sub.labels == ["B2"]
    @test size(sub.curves) == (1, 3)
    @test sub.curves[1, :] == [4.0, 5.0, 6.0]
    @test sub.times == times
    @test sub.clusters === nothing
    @test sub.centroids === nothing

    # multiple wells, different order than original
    sub2 = gd[["C3", "A1"]]
    @test sub2.labels == ["C3", "A1"]
    @test sub2.curves[1, :] == [7.0, 8.0, 9.0]
    @test sub2.curves[2, :] == [1.0, 2.0, 3.0]

    # unknown label throws ArgumentError
    @test_throws ArgumentError gd[["B2", "Z9"]]

    # clustering fields are dropped even if present on source
    gd_with_clusters = GrowthData(curves, times, labels, [1, 2, 3])
    sub3 = gd_with_clusters[["A1"]]
    @test sub3.clusters === nothing
end

# 1.5 — Custom model trait constructors
@testset "Custom model constructors (NLModel, ODEModel, LogLinModel)" begin

    @testset "NLModel without guess" begin
        f(p, t) = p[1] ./ (1 .+ exp.(-p[2] .* t))
        m = NLModel("my_nl", f, ["K", "gr"])
        @test m isa AbstractGrowthModel
        @test m.name == "my_nl"
        @test m.param_names == ["K", "gr"]
        @test m.guess === nothing
    end

    @testset "NLModel with guess function" begin
        f(p, t) = p[1] .* exp.(p[2] .* t)
        g(data) = [data[2, end], 0.3]
        m = NLModel("my_exp_nl", f, ["A", "gr"], g)
        @test m.guess === g
    end

    @testset "ODEModel" begin
        f!(du, u, p, t) = (du[1] = p[1] * u[1] * (1 - u[1] / p[2]))
        m = ODEModel("my_ode", f!, ["gr", "K"], 1)
        @test m isa AbstractGrowthModel
        @test m.name == "my_ode"
        @test m.n_eq == 1
        @test m.guess === nothing
    end

    @testset "LogLinModel" begin
        m = LogLinModel()
        @test m isa AbstractGrowthModel
        @test m isa LogLinModel
    end
end
