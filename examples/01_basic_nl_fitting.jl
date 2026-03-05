# =============================================================================
# Example 1 — Basic NL Model Fitting
# =============================================================================
# Shows the minimal usage of the new unified API:
#   1. Build a GrowthData from a matrix
#   2. Pick models from MODEL_REGISTRY
#   3. Call fit() and inspect results
#
# Compare the old API at the bottom to see the difference.
# =============================================================================

using Kinbiont
using OrdinaryDiffEq          # for ODE_sim (data generation)
using Distributions            # for noise generation
using Random

Random.seed!(42)

# ---------------------------------------------------------------------------
# Generate synthetic data with the existing ODE_sim helper
# (logistic growth + uniform noise on 3 replicate curves)
# ---------------------------------------------------------------------------

tstart, tmax, delta_t = 0.0, 48.0, 1.0
true_params = [0.25, 1.5]   # gr, N_max  (logistic ODE)

times = collect(tstart:delta_t:tmax)
n_tp  = length(times)

# Simulate 3 replicate curves with small noise
n_curves = 3
curves   = Matrix{Float64}(undef, n_curves, n_tp)

for i in 1:n_curves
    sim = ODE_sim("logistic", [0.05], tstart, tmax, delta_t, true_params)
    noise = rand(Uniform(-0.02, 0.02), n_tp)
    curves[i, :] = reduce(hcat, sim.u)[1, :] .+ noise
end

labels = ["replicate_$i" for i in 1:n_curves]

# ---------------------------------------------------------------------------
# New API — 3 structs + 1 function call
# ---------------------------------------------------------------------------

# 1. Wrap data
data = GrowthData(curves, times, labels)

# 2. Specify candidate models with initial guesses
#    MODEL_REGISTRY holds all built-in NL and ODE models
spec = ModelSpec(
    [MODEL_REGISTRY["NL_logistic"], MODEL_REGISTRY["NL_Gompertz"], MODEL_REGISTRY["NL_exponential"]],
    [
        [1.5,  0.1, 0.05],    # NL_logistic:   N_max, growth_rate, lag
        [1.5,  0.3, 20.0],    # NL_Gompertz:   N_max, growth_rate, lag
        [0.05, 0.25],         # NL_exponential: N0, growth_rate
    ],
)

# 3. Configure fitting (all defaults — just override loss)
opts = FitOptions(loss = "RE")

# 4. Fit — one call, all curves, AICc model selection built in
results = kinbiont_fit(data, spec, opts)

# ---------------------------------------------------------------------------
# Inspect results
# ---------------------------------------------------------------------------

println("=== New API results ===")
for r in results
    println(
        "  $(r.label):  best=$(r.best_model.name)  " *
        "AICc=$(round(r.best_aic; digits=2))  " *
        "loss=$(round(r.loss; digits=4))"
    )
    println("    params: ", round.(Float64.(r.best_params); digits=4))
end

# ---------------------------------------------------------------------------
# Old API — for comparison: same fit, one curve at a time, more boilerplate
# ---------------------------------------------------------------------------

println("\n=== Old API (equivalent for replicate_1) ===")

data_mat_old = Matrix(transpose(hcat(times, curves[1, :])))

nl_result = fit_NL_model(
    data_mat_old,
    "replicate_1",
    "example",
    NL_model_logistic,     # pass the function directly
    [1.5, 0.1, 0.05];
    type_of_loss = "RE",
)

println("  method: ", nl_result[1])
println("  loss:   ", nl_result[2][end])

# ---------------------------------------------------------------------------
# Save results to CSV
# ---------------------------------------------------------------------------

files = save_results(results, "output/example_01"; prefix = "basic_nl")
println("\n=== Saved files ===")
for (k, v) in pairs(files)
    println("  $k → $v")
end
