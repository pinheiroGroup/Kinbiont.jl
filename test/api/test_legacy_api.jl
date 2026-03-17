# Group 4 — Legacy (old) single-curve fitting API
# Tests each legacy fitting function independently, verifying return structure
# and finite output rather than exact parameter recovery (optimizer is stochastic).

@testset "Legacy single-curve fitting" begin

    # Shared synthetic logistic data (same generation as test_fitting.jl / test_vs_legacy.jl)
    sim   = ODE_sim("logistic", [0.05], 0.0, 24.0, 1.0, [0.3, 1.2])
    times = collect(0.0:1.0:24.0)
    curve = sim[1, :]
    # Old API expects 2×N matrix: first row = times, second row = OD
    data_mat = Matrix(transpose(hcat(times, curve)))

    # ------------------------------------------------------------------
    @testset "4.1 fitting_one_well_ODE_constrained — return structure" begin
        Random.seed!(42)
        result = fitting_one_well_ODE_constrained(
            data_mat,
            "well_A",
            "exp1",
            "logistic",
            [0.3, 1.2];
            lb = [0.01, 0.1], ub = [2.0, 3.0],
        )

        @test result[1] == "ODE"
        # res_param layout: [label_exp, well, model, param..., th_gr, em_gr, loss]
        @test result[2][2] == "well_A"
        @test result[2][3] == "logistic"
        @test isfinite(Float64(result[2][end]))       # loss is finite
        @test length(result[3]) == length(result[4])  # fitted values match times
        @test all(isfinite.(Float64.(result[3])))
    end

    # ------------------------------------------------------------------
    @testset "4.2 fitting_one_well_Log_Lin — return structure and R²" begin
        result = fitting_one_well_Log_Lin(
            data_mat,
            "well_A",
            "exp1";
            type_of_smoothing = "NO",
        )

        @test result[1] == "Log-lin"
        # R² is stored at params[end]
        r2 = Float64(result[2][end])
        @test isfinite(r2)
        @test 0.0 <= r2 <= 1.0
        # fit_matrix is n×2: [times  log_fitted_values]
        @test size(result[3], 2) == 2
        @test length(result[3][:, 1]) > 0
    end

    # ------------------------------------------------------------------
    @testset "4.3 fit_NL_model (NL_logistic) — return structure" begin
        Random.seed!(42)
        result = fit_NL_model(
            data_mat,
            "well_A",
            "exp1",
            NL_model_logistic,
            [1.5, 0.3, 0.0];
            lb = [0.1, 0.01, 0.0], ub = [3.0, 1.0, 2.0],
        )

        @test result[1] == "NL"
        # res_param: [label_exp, well, model_name, params..., th_gr, em_gr, loss]
        @test isfinite(Float64(result[2][end]))        # loss
        @test length(result[3]) == length(times)       # fitted curve length
        @test all(isfinite.(Float64.(result[3])))
    end

    # ------------------------------------------------------------------
    @testset "4.4 ODE_Model_selection — selects one of the candidates" begin
        Random.seed!(42)
        result = ODE_Model_selection(
            data_mat,
            "well_A",
            "exp1",
            ["logistic", "gompertz"],
            [[0.3, 1.2], [0.3, 1.2]];
            lb_param_array = [[0.01, 0.1], [0.01, 0.1]],
            ub_param_array = [[2.0, 3.0],  [2.0, 3.0]],
            maxiters = 2_000,
        )

        @test result[1] == "ODE_model_selection"
        @test result[9] in ["logistic", "gompertz"]    # result[9] = selected model name
        @test all(isfinite.(Float64.(result[3])))       # fitted values finite
        @test length(result[4]) > 0                    # times non-empty
    end

    # ------------------------------------------------------------------
    @testset "4.5 NL_model_selection — selects one of the candidates" begin
        Random.seed!(42)
        result = NL_model_selection(
            data_mat,
            "well_A",
            "exp1",
            [NL_model_logistic, NL_model_Gompertz],
            [[1.5, 0.3, 0.0], [1.5, 0.3, 20.0]];
            lb_param_array = [[0.1, 0.0, 0.0],  [0.1, 0.0, 0.0]],
            ub_param_array = [[3.0, 2.0, 2.0],  [3.0, 2.0, 50.0]],
        )

        @test result[1] == "NL_model_selection"
        @test all(isfinite.(Float64.(result[3])))       # fitted values
        @test isfinite(Float64(result[6]))              # top_loss
    end

    # ------------------------------------------------------------------
    @testset "4.6 fitting_one_well_custom_ODE — user-supplied ODE" begin
        # Re-implement logistic as a custom function to verify the custom ODE path
        function custom_logistic!(du, u, p, t)
            du[1] = (p[1] / p[2]) * u[1] * (p[2] - u[1])
        end

        Random.seed!(42)
        result = fitting_one_well_custom_ODE(
            data_mat,
            "well_A",
            "exp1",
            custom_logistic!,
            [0.3, 1.2],
            1;                                   # n_equation
            lb = [0.01, 0.1], ub = [2.0, 3.0],
        )

        @test result[1] == "custom_ODE"
        # res_param: [name_well, "custom_model", fitted_params_vec, th_gr, em_gr, loss]
        @test isfinite(Float64(result[2][end]))
        @test length(result[3]) == length(result[4])
        @test all(isfinite.(Float64.(result[3])))
    end
end
