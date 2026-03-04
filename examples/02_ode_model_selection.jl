# =============================================================================
# Example 2 — ODE Model Selection across Multiple Curves
# =============================================================================
# Demonstrates:
#   - Mixing NLModel and ODEModel candidates in one ModelSpec
#   - Using custom ODE functions alongside registry models
#   - Fitting a batch of curves with a single fit() call
#   - Accessing per-model details via all_results
# =============================================================================

using Kinbiont
using OrdinaryDiffEq
using Distributions
using Random

Random.seed!(7)

# ---------------------------------------------------------------------------
# Generate synthetic data: 4 curves, each from a different true model
# ---------------------------------------------------------------------------

tstart, tmax, delta_t = 0.0, 72.0, 2.0
times = collect(tstart:delta_t:tmax)
n_tp  = length(times)

# True models and parameters
scenarios = [
    ("logistic",                     [0.20, 1.0]),
    ("alogistic",                    [0.18, 1.0, 1.5]),
    ("baranyi_roberts",              [0.22, 1.0, 8.0, 2.0, 1.5]),
    ("triple_piecewise_adjusted_logistic", [0.18, 1.0, 12.0, 1.5, 0.002, 60.0, -0.001]),
]

n_curves = length(scenarios)
curves   = Matrix{Float64}(undef, n_curves, n_tp)
labels   = String[]

for (i, (model_name, params)) in enumerate(scenarios)
    sim = ODE_sim(model_name, [0.05], tstart, tmax, delta_t, params)
    noise = rand(Uniform(-0.015, 0.015), n_tp)
    curves[i, :] = clamp.(reduce(hcat, sim.u)[1, :] .+ noise, 0.0, Inf)
    push!(labels, "curve_$(i)_true=$(model_name)")
end

# ---------------------------------------------------------------------------
# New API: define a custom ODE and mix it with registry models
# ---------------------------------------------------------------------------

# A simple custom ODE (modified logistic with Allee effect)
function allee_logistic!(du, u, p, t)
    du[1] = p[1] * u[1] * (u[1] / p[3] - 1.0) * (1.0 - u[1] / p[2])
end

custom_allee = ODEModel(
    "allee_logistic",
    allee_logistic!,
    ["r", "K", "allee_threshold"],
    1,           # n_eq — explicit, no detection loop
    nothing,     # no guess function; we provide params in ModelSpec
)

# Candidate models: mix of NL, ODE from registry, and a custom ODE
candidate_models = [
    MODEL_REGISTRY["NL_logistic"],
    MODEL_REGISTRY["NL_Gompertz"],
    MODEL_REGISTRY["logistic"],         # ODE version
    MODEL_REGISTRY["alogistic"],
    MODEL_REGISTRY["baranyi_roberts"],
    custom_allee,
]

initial_params = [
    [1.0, 0.2, 0.05],               # NL_logistic:   N_max, gr, lag
    [1.0, 0.2, 30.0],               # NL_Gompertz:   N_max, gr, lag
    [0.2, 1.0],                      # logistic ODE:  gr, N_max
    [0.2, 1.0, 1.5],                 # alogistic:     gr, N_max, shape
    [0.2, 1.0, 8.0, 2.0, 1.5],      # baranyi_roberts: gr, N_max, lag, shape1, shape2
    [0.2, 1.0, 0.1],                 # allee_logistic: r, K, allee_threshold
]

spec = ModelSpec(candidate_models, initial_params)

# Light smoothing, no blank correction
opts = FitOptions(smooth = true, smooth_method = :lowess, lowess_frac = 0.1, loss = "L2")

data = GrowthData(curves, times, labels)

# ---------------------------------------------------------------------------
# Single call — all curves × all models × AICc selection
# ---------------------------------------------------------------------------

results = kinbiont_fit(data, spec, opts)

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

println("=== ODE model selection results ===\n")
for r in results
    println("  $(r.label)")
    println("    → selected: $(r.best_model.name)  AICc=$(round(r.best_aic; digits=2))")
    println("    → params:   ", round.(Float64.(r.best_params); digits=4))

    # Show AICc ranking of all candidates for the first curve
    if r === results[1]
        println("\n  AICc ranking for $(r.label):")
        sorted = sort(r.all_results; by = x -> x.aic)
        for (rank, candidate) in enumerate(sorted)
            marker = candidate.model_name == r.best_model.name ? " ← best" : ""
            println("    $(rank). $(candidate.model_name)  AICc=$(round(candidate.aic; digits=2))$marker")
        end
    end
    println()
end
