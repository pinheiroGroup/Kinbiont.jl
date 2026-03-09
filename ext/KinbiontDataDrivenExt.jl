"""
    KinbiontDataDrivenExt

Extension that upgrades the `DDDEModel` backend to the full
`DataDrivenDiffEq` pipeline when all three optional packages are loaded:

    using DataDrivenDiffEq, DataDrivenSparse, ModelingToolkit

Advantages over the built-in pure-Julia fallback:
- Uses DataDrivenDiffEq's validated STLSQ implementation.
- Returns the discovered symbolic equation (e.g. `[φ₁ ~ 0.3·u - 0.25·u²]`)
  as a human-readable string in `CurveFitResult.best_params[1]`.
- Uses DataDrivenDiffEq's native AICc (StatsBase.aicc) as the model score.
"""
module KinbiontDataDrivenExt

using Kinbiont
using DataDrivenDiffEq
using DataDrivenSparse
using ModelingToolkit

# Called by Kinbiont._ddde_backend via Base.get_extension when this
# module is loaded.  Must NOT override any method from the base package.
function _ddde_backend_dde(
    t::Vector{Float64},
    y::Vector{Float64},
    model::Kinbiont.DDDEModel,
)
    N = length(t)

    # 1) Numerical derivative (forward / central / backward differences)
    dydt = Kinbiont._ddde_finite_differences(t, y)

    # 2) DataDrivenDiffEq Direct problem: state X → derivative Y
    X       = reshape(y,    1, :)
    Y       = reshape(dydt, 1, :)
    problem = DirectDataDrivenProblem(X, Y)

    # 3) Polynomial basis in a single symbolic variable u
    u     = ModelingToolkit.variable(:u)
    basis = Basis(polynomial_basis([u], model.max_degree), [u])

    # 4) STLSQ sparse regression
    lambdas = exp10.(model.lambda_min:model.lambda_step:model.lambda_max)
    opt     = STLSQ(lambdas)
    ddsol   = solve(problem, basis, opt;
                    options = DataDrivenCommonOptions(digits = 3))

    # 5) Extract results
    sys  = get_basis(ddsol)
    pmap = get_parameter_map(sys)                  # (param => value) pairs
    n_active = max(length(pmap), 1)

    # Build coefficient vector aligned with the polynomial basis [1, u, u², …]
    # DataDrivenDiffEq keeps all basis terms as parameters (non-zero or zero).
    xi = zeros(Float64, model.max_degree + 1)
    for (k, (_, v)) in enumerate(pmap)
        k <= length(xi) && (xi[k] = v)
    end

    # Human-readable symbolic equation — the main advantage of using DDE
    eq_str = string(equations(sys))

    return (
        xi           = xi,
        rss          = DataDrivenDiffEq.StatsBase.rss(ddsol),
        n_active     = n_active,
        params_out   = Vector{Any}([eq_str; pmap]),
        aic_override = Float64(DataDrivenDiffEq.StatsBase.aicc(ddsol)),
    )
end

end # module KinbiontDataDrivenExt
