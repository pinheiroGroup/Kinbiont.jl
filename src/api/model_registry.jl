# =============================================================================
# Unified Model Registry
# =============================================================================
# Wraps the existing JModel / NL_Model definitions into the new NLModel /
# ODEModel trait types used for dispatch.
#
# Single source of truth: MODEL_REGISTRY["name"] -> AbstractGrowthModel
#
# ODE models have n_eq explicitly declared here so there is no need for the
# error-prone detection-loop that probes up to 100 state dimensions.
#
# NL models carry their existing guess functions; ODE models currently have
# no guess (pass initial params explicitly via ModelSpec).
# =============================================================================

# ---------------------------------------------------------------------------
# Helpers to strip result-column metadata from the legacy param vectors
# ---------------------------------------------------------------------------

# ODE params format: ["label_exp", "well", "model", p1, p2, ..., "th_max_gr", "emp_max_gr", "loss"]
_ode_param_names(jmodel::JModel) = jmodel.params[4:end-3]

# NL params format: ["model", "well", p1, p2, ..., "th_max_gr", "emp_max_gr", "loss"]
_nl_param_names(nlm::NL_Model) = nlm.params[3:end-3]

# ---------------------------------------------------------------------------
# Wrap each ODE JModel into an ODEModel
# (n_eq determined from the ODE functions' state variable usage)
# ---------------------------------------------------------------------------

const _ODE_N_EQ = Dict{String,Int}(
    # 3-state HPM variants
    "HPM_3_inhibition"                   => 3,
    "HPM_3_death"                        => 3,
    "aHPM_3_death_resistance"            => 3,
    # 2-state HPM variants
    "HPM"                                => 2,
    "HPM_exp"                            => 2,
    "HPM_inhibition"                     => 2,
    "aHPM_inhibition"                    => 2,
    "ODEs_HPM_SR"                        => 2,
    "aHPM"                               => 2,
    # 1-state models (all others)
    "Diauxic_replicator_1"               => 1,
    "Diauxic_replicator_2"               => 1,
    "Diauxic_piecewise_adjusted_logistic"=> 1,
    "exponential"                        => 1,
    "piecewise_adjusted_logistic"        => 1,
    "triple_piecewise_adjusted_logistic" => 1,
    "triple_piecewise"                   => 1,
    "triple_piecewise_sub_linear"        => 1,
    "gbsm_piecewise"                     => 1,
    "ODE_four_piecewise"                 => 1,
    "hyper_gompertz"                     => 1,
    "hyper_logistic"                     => 1,
    "bertalanffy_richards"               => 1,
    "ode_von_bertalanffy"                => 1,
    "triple_piecewise_bertalanffy_richards" => 1,
    "logistic"                           => 1,
    "alogistic"                          => 1,
    "gompertz"                           => 1,
    "baranyi_richards"                   => 1,
    "baranyi_exp"                        => 1,
    "baranyi_roberts"                    => 1,
)

# Build ODE entries from the existing models_list
_ode_entries = [
    ODEModel(
        jm.name,
        jm.func,
        _ode_param_names(jm),
        _ODE_N_EQ[jm.name],
    )
    for jm in models_list
]

# Build NL entries from the existing NL_models_list
_nl_entries = [
    NLModel(
        nlm.name,
        nlm.func,
        _nl_param_names(nlm),
        nlm.guess,
    )
    for nlm in NL_models_list
]

# ---------------------------------------------------------------------------
# Public registry — Dict{String, AbstractGrowthModel}
# ---------------------------------------------------------------------------

"""
    MODEL_REGISTRY

Dictionary mapping model name strings to their [`AbstractGrowthModel`](@ref) instances.

Use this to build a [`ModelSpec`](@ref):

```julia
spec = ModelSpec(
    [MODEL_REGISTRY["NL_logistic"], MODEL_REGISTRY["logistic"]],
    [[1.2, 0.5, 0.01],              [0.5, 1.2]],
)
```

Keys with the `"NL_"` prefix are [`NLModel`](@ref)s; all others are [`ODEModel`](@ref)s.
The sentinel [`LogLinModel`](@ref) is not in this registry — construct it directly.
"""
const MODEL_REGISTRY = Dict{String, AbstractGrowthModel}(
    m.name => m
    for m in Iterators.flatten((_nl_entries, _ode_entries))
)

export MODEL_REGISTRY
