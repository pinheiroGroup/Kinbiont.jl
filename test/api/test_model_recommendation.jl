@testset "model_recommendation" begin

    # -----------------------------------------------------------------------
    # 1. Build a small DB (n_samples=10 for speed)
    # -----------------------------------------------------------------------
    db = build_model_fingerprint_db(n_samples=10, tmax=24.0, n_points=20)

    @test db isa ModelFingerprintDB
    @test !isempty(db.model_names)
    @test size(db.features, 2) == 15
    @test all(isfinite, db.features)

    # Every feature row should be (approximately) unit-norm
    for i in axes(db.features, 1)
        @test sqrt(sum(db.features[i, :].^2)) ≈ 1.0 atol=1e-6
    end

    # All model names in the DB should be in MODEL_REGISTRY
    @test all(n -> haskey(MODEL_REGISTRY, n), db.model_names)

    # aHPM (n_eq=2) should appear in the DB
    @test "aHPM" in db.model_names

    # -----------------------------------------------------------------------
    # 2. recommend_models
    # -----------------------------------------------------------------------
    sim  = ODE_sim("logistic", [0.05], 0.0, 24.0, 0.5, [0.3, 1.2])
    t    = sim.t
    y    = [u[1] for u in sim.u]

    recs = recommend_models(y, t, db; top_k=3)

    @test recs isa Vector{String}
    @test 1 <= length(recs) <= 3
    @test all(n -> haskey(MODEL_REGISTRY, n), recs)

    # Recommendations are unique
    @test length(recs) == length(unique(recs))

    # -----------------------------------------------------------------------
    # 3. smart_fit
    # -----------------------------------------------------------------------
    data = GrowthData(reshape(y, 1, :), t, ["curve1"])
    res  = smart_fit(data, db; top_k=2)

    @test res isa GrowthFitResults
    @test length(res) == 1
    @test res[1].label == "curve1"
    @test res[1].best_model isa AbstractGrowthModel

    # -----------------------------------------------------------------------
    # 4. Edge cases
    # -----------------------------------------------------------------------
    # Empty DB should return no recommendations
    empty_db = ModelFingerprintDB(String[], AbstractGrowthModel[],
                                  Matrix{Float64}(undef, 0, 15),
                                  db.feature_names, Vector{Float64}[])
    @test recommend_models(y, t, empty_db; top_k=3) == String[]

    # Feature extraction is stable on monotone, flat, and noisy curves
    t2 = collect(range(0.0, 24.0; length=30))
    @test all(isfinite, Kinbiont._extract_features(t2, ones(30)))
    @test all(isfinite, Kinbiont._extract_features(t2, collect(range(0.01, 1.2; length=30))))

end
