# aHPM (adjusted Heterogeneous Population Model) tests
#
# Biology: cells start in a lag (dormant) pool u[1] and transition to an
# actively-growing pool u[2] at rate exit_lag_rate.
# The ODE is ODEs_adjusted_McKellar:
#   du[1] = -u[1] * exit_lag_rate
#   du[2] =  u[1] * exit_lag_rate
#            + gr * u[2] * (1 - ((u[1]+u[2])/N_max)^shape)
#
# Observable: total biomass = u[1] + u[2]  (Kinbiont sums states automatically)
# Parameters: [gr, exit_lag_rate, N_max, shape]  (4 params, n_eq = 2)
#
# Canonical parameterisation used throughout:
#   gr=0.5, exit_lag_rate=0.15, N_max=1.5, shape=1.0
# The slow exit_lag_rate (0.15) produces a clearly visible lag phase.

const AHPM_PARAMS  = [0.5, 0.15, 1.5, 1.0]
const AHPM_N0      = 0.05    # initial total OD (all in lag pool)
const AHPM_TMAX    = 30.0
const AHPM_DT      = 0.5
const AHPM_TIMES   = collect(0.0:AHPM_DT:AHPM_TMAX)

# Helper: simulate and return the 2×N data matrix for old-API calls
function _ahpm_data_mat()
    sim = ODE_sim("aHPM", [AHPM_N0, 0.0], 0.0, AHPM_TMAX, AHPM_DT, AHPM_PARAMS)
    total = sim[1, :] .+ sim[2, :]          # u[1] + u[2] = observable
    return Matrix(transpose(hcat(AHPM_TIMES, total)))  # 2×N
end

# ─────────────────────────────────────────────────────────────────────────────
@testset "aHPM model" begin

    # ------------------------------------------------------------------
    @testset "Simulation shape and biology" begin
        sim   = ODE_sim("aHPM", [AHPM_N0, 0.0], 0.0, AHPM_TMAX, AHPM_DT, AHPM_PARAMS)
        lag   = sim[1, :]   # lag-cell compartment
        grow  = sim[2, :]   # growing-cell compartment
        total = lag .+ grow

        # Output dimensions: 2 states × n_timepoints
        @test length(sim.t) == length(AHPM_TIMES)
        @test all(isfinite.(lag))
        @test all(isfinite.(grow))

        # Lag cells decay monotonically
        @test all(diff(lag) .<= 1e-6)

        # Total biomass converges toward N_max
        @test total[end] ≈ AHPM_PARAMS[3] atol=0.05

        # Lag phase is visible: total biomass grows slowly at first.
        # Slope in the first 20 % of the time course should be less than
        # slope in the middle 20 % (exponential phase).
        n = length(total)
        early_slope  = (total[div(n, 5)] - total[1]) / (AHPM_TIMES[div(n, 5)])
        middle_slope = (total[div(n, 2)] - total[div(2*n, 5)]) /
                       (AHPM_TIMES[div(n, 2)] - AHPM_TIMES[div(2*n, 5)])
        @test early_slope < middle_slope
    end

    # ------------------------------------------------------------------
    @testset "Old API: fitting_one_well_ODE_constrained" begin
        data_mat = _ahpm_data_mat()
        Random.seed!(42)

        result = fitting_one_well_ODE_constrained(
            data_mat,
            "well_A",
            "exp1",
            "aHPM",
            AHPM_PARAMS;                           # initial guess = true params
            lb = [0.01, 0.01, 0.5, 0.1],
            ub = [2.0,  2.0,  3.0, 5.0],
        )

        @test result[1] == "ODE"
        # res_param: [label_exp, well, "aHPM", gr, exit_lag_rate, N_max, shape, th_gr, em_gr, loss]
        @test result[2][3] == "aHPM"
        @test isfinite(Float64(result[2][end]))       # loss
        @test length(result[3]) == length(result[4])  # curve length matches times
        @test all(isfinite.(Float64.(result[3])))
        # Fitted curve should reach N_max at steady state
        @test Float64(result[3][end]) ≈ AHPM_PARAMS[3] atol=0.15
    end

    # ------------------------------------------------------------------
    @testset "New API: kinbiont_fit with MODEL_REGISTRY[\"aHPM\"]" begin
        sim    = ODE_sim("aHPM", [AHPM_N0, 0.0], 0.0, AHPM_TMAX, AHPM_DT, AHPM_PARAMS)
        total  = sim[1, :] .+ sim[2, :]
        data   = GrowthData(reshape(total, 1, :), AHPM_TIMES, ["well_A"])

        spec = ModelSpec(
            [MODEL_REGISTRY["aHPM"]],
            [AHPM_PARAMS];                         # initial guess = true params
            lower = [[0.01, 0.01, 0.5, 0.1]],
            upper = [[2.0,  2.0,  3.0, 5.0]],
        )

        Random.seed!(42)
        res = kinbiont_fit(data, spec, FitOptions(loss="RE"))

        @test res isa GrowthFitResults
        @test length(res) == 1
        @test res[1].best_model isa ODEModel
        @test res[1].best_model.name == "aHPM"
        @test length(res[1].best_params) == 4          # gr, exit_lag_rate, N_max, shape
        @test isfinite(res[1].best_aic)
        @test res[1].loss >= 0
        # Fitted curve converges toward N_max
        @test res[1].fitted_curve[end] ≈ AHPM_PARAMS[3] atol=0.15
    end

    # ------------------------------------------------------------------
    @testset "New API: multi-curve GrowthData" begin
        sim1 = ODE_sim("aHPM", [AHPM_N0, 0.0], 0.0, AHPM_TMAX, AHPM_DT, AHPM_PARAMS)
        sim2 = ODE_sim("aHPM", [AHPM_N0, 0.0], 0.0, AHPM_TMAX, AHPM_DT,
                       [0.4, 0.2, 1.2, 1.0])   # different params for second curve
        c1 = sim1[1, :] .+ sim1[2, :]
        c2 = sim2[1, :] .+ sim2[2, :]
        data = GrowthData(vcat(reshape(c1, 1, :), reshape(c2, 1, :)),
                          AHPM_TIMES, ["A", "B"])

        spec = ModelSpec(
            [MODEL_REGISTRY["aHPM"]],
            [AHPM_PARAMS];
            lower = [[0.01, 0.01, 0.5, 0.1]],
            upper = [[2.0,  2.0,  3.0, 5.0]],
        )

        Random.seed!(42)
        res = kinbiont_fit(data, spec, FitOptions(loss="RE"))

        @test length(res) == 2
        for r in res
            @test r.best_model.name == "aHPM"
            @test isfinite(r.best_aic)
            @test all(isfinite.(r.fitted_curve))
        end
    end

    # ------------------------------------------------------------------
    @testset "Model selection: aHPM wins over logistic on lag-phase data" begin
        # aHPM-generated data has a clear lag phase that logistic cannot model.
        # With correct bounds and enough data, aHPM should have lower AICc.
        data_mat = _ahpm_data_mat()
        total    = data_mat[2, :]
        data     = GrowthData(reshape(total, 1, :), AHPM_TIMES, ["well_A"])

        spec = ModelSpec(
            [MODEL_REGISTRY["aHPM"], MODEL_REGISTRY["logistic"]],
            [AHPM_PARAMS,            [0.3, 1.2]];
            lower = [[0.01, 0.01, 0.5, 0.1], [0.01, 0.1]],
            upper = [[2.0,  2.0,  3.0, 5.0],  [2.0,  3.0]],
        )

        Random.seed!(42)
        res = kinbiont_fit(data, spec, FitOptions(loss="RE"))

        @test res[1].best_model.name == "aHPM"

        ahpm_aic    = first(r.aic for r in res[1].all_results if r.model_name == "aHPM")
        logistic_aic = first(r.aic for r in res[1].all_results if r.model_name == "logistic")
        @test ahpm_aic < logistic_aic
    end
end
