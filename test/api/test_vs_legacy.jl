@testset "New API == Legacy API" begin
    # Shared test data
    Random.seed!(42)
    sim   = ODE_sim("logistic", [0.05], 0.0, 24.0, 1.0, [0.3, 1.2])
    times = collect(0.0:1.0:24.0)
    curve = sim[1, :]

    # Legacy format: 2×n matrix [times; values]
    data_mat = Matrix(transpose(hcat(times, curve)))

    # New API data container (single curve)
    data_g = GrowthData(reshape(curve, 1, :), times, ["curve1"])

    # -----------------------------------------------------------------
    @testset "NL logistic equivalence" begin
        p0 = [1.5, 0.3, 0.0]

        # Set seed, run legacy
        Random.seed!(7)
        raw_nl = fit_NL_model(
            data_mat,
            "NL_logistic",   # name_well
            "curve1",        # label_exp
            NL_model_logistic,
            p0;
            type_of_loss = "RE",
            smoothing    = false,
        )
        old_params = raw_nl[2][4:end-3]
        old_loss   = Float64(raw_nl[2][end])
        old_curve  = Vector{Float64}(raw_nl[3])

        # Same seed, run new API
        Random.seed!(7)
        spec_nl = ModelSpec([MODEL_REGISTRY["NL_logistic"]], [p0])
        res_nl  = kinbiont_fit(data_g, spec_nl, FitOptions(loss="RE"))

        @test Float64.(res_nl[1].best_params) == Float64.(old_params)
        @test res_nl[1].loss                  == old_loss
        @test res_nl[1].fitted_curve          == old_curve
    end

    # -----------------------------------------------------------------
    @testset "ODE logistic equivalence" begin
        p0_ode = [0.3, 1.2]

        # Set seed, run legacy
        Random.seed!(7)
        raw_ode = fitting_one_well_custom_ODE(
            data_mat,
            "logistic",   # name_well
            "curve1",     # label_exp
            logistic,
            p0_ode,
            1;            # n_equation
            smoothing    = false,
            type_of_loss = "RE",
        )
        old_params_ode = raw_ode[2][3]        # fitted parameter vector
        old_loss_ode   = Float64(raw_ode[2][end])
        old_curve_ode  = Vector{Float64}(raw_ode[3])

        # Same seed, run new API
        Random.seed!(7)
        spec_ode = ModelSpec([MODEL_REGISTRY["logistic"]], [p0_ode])
        res_ode  = kinbiont_fit(data_g, spec_ode, FitOptions(loss="RE"))

        @test Float64.(res_ode[1].best_params) == Float64.(old_params_ode)
        @test res_ode[1].loss                  == old_loss_ode
        @test res_ode[1].fitted_curve          == old_curve_ode
    end

    # -----------------------------------------------------------------
    @testset "Log-linear equivalence" begin
        # Legacy call (same arguments the new API wrapper uses; no randomness in log-lin)
        raw_ll = fitting_one_well_Log_Lin(
            data_mat,
            "log_lin",
            "curve1";
            type_of_smoothing       = "NO",
            pt_smoothing_derivative = 7,
        )
        # raw_ll is a Tuple: (method, params, fit_matrix, smoothed_data, confidence_band)
        # params contains mixed types (strings + floats)
        old_params_ll = raw_ll[2]
        old_curve_ll  = raw_ll[3][:, 2]    # log-scale fitted values (2nd column)

        # New API call (no randomness — seed irrelevant)
        spec_ll = ModelSpec([LogLinModel()], [Float64[]])
        res_ll  = kinbiont_fit(data_g, spec_ll)

        @test res_ll[1].best_params == old_params_ll
        @test res_ll[1].fitted_curve == old_curve_ll
    end

    # -----------------------------------------------------------------
    @testset "AICc consistency" begin
        p0 = [1.5, 0.3, 0.0]

        Random.seed!(7)
        raw_nl = fit_NL_model(
            data_mat,
            "NL_logistic",
            "curve1",
            NL_model_logistic,
            p0;
            type_of_loss = "RE",
            smoothing    = false,
        )
        old_loss  = Float64(raw_nl[2][end])
        old_times = Vector{Float64}(raw_nl[4])
        old_aic   = Kinbiont.AICc_evaluation2(length(p0), 2.0, old_times, old_loss;
                                              correction=true)

        Random.seed!(7)
        spec_nl = ModelSpec([MODEL_REGISTRY["NL_logistic"]], [p0])
        res_nl  = kinbiont_fit(data_g, spec_nl, FitOptions(loss="RE"))

        @test res_nl[1].best_aic == Float64(old_aic)
    end
end
