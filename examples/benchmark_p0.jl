# benchmark_p0.jl
#
# Benchmarks three fitting strategies for a chosen plate + well:
#
#   A) BBO   + default_p0  (bounds midpoint, wide search) ← baseline
#   B) BBO   + narrow bounds around suggest_p0            ← option B
#   C) LBFGS + suggest_p0  (gradient-based, warm start)  ← option A
#
# Each strategy is run N_RUNS times with fixed seeds so that timing and
# AICc variability reflect the optimiser behaviour, not randomness in p0.
#
# Usage (interactive):
#   julia --project examples/benchmark_p0.jl
#
# Usage (non-interactive):
#   julia --project examples/benchmark_p0.jl <folder> <well> [model1,model2,…]
#   julia --project examples/benchmark_p0.jl /data/LG177 A1
#   julia --project examples/benchmark_p0.jl /data/LG177 A1 aHPM,logistic

using Kinbiont
using Random
using Serialization
using Statistics: mean, std
using Logging
using OptimizationBBO: BBO_adaptive_de_rand_1_bin_radiuslimited
using OptimizationOptimJL: LBFGS

global_logger(ConsoleLogger(stderr, Logging.Error))

const N_RUNS = 3
const SEEDS  = [42, 43, 44]

# ── narrow-bounds factor around suggest_p0 ────────────────────────────────────
# BBO searches within [p * BOUND_LO, p * BOUND_HI] per parameter.
# A tighter window speeds up BBO by reducing wasted exploration.
const BOUND_LO = 0.2    # lower = 20% of suggested value
const BOUND_HI = 5.0    # upper = 5× suggested value

# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────
function find_stationary_phase(
    times::Vector{Float64}, od::Vector{Float64};
    percentile_thr::Float64 = 0.05, pt_smooth_derivative::Int = 7,
    win_size::Int = 5, thr_od::Float64 = 0.001,
)::Union{Float64, Nothing}
    data_mat = Matrix(transpose(hcat(times, od)))
    index_od = findall(data_mat[2, :] .> thr_od)
    isempty(index_od) && return nothing
    data_t = data_mat[:, index_od]
    length(data_t[1, :]) <= pt_smooth_derivative && return nothing
    sgr = Kinbiont.specific_gr_evaluation(data_t, pt_smooth_derivative)
    isempty(sgr) && return nothing
    max_sgr = maximum(sgr)
    max_sgr <= 0 && return nothing
    thr = max_sgr * percentile_thr
    max_i = argmax(sgr)
    for i in max_i:(length(sgr) - win_size)
        all(sgr[i:(i+win_size-1)] .< thr) && return times[min(i + index_od[1] - 1, length(times))]
    end
    return nothing
end

function ask(prompt::String)::String
    print(prompt); strip(readline())
end

# ─────────────────────────────────────────────────────────────────────────────
# Resolve inputs
# ─────────────────────────────────────────────────────────────────────────────
folder = if length(ARGS) >= 1
    strip(ARGS[1])
else
    f = ""
    while !isdir(f); f = ask("Experiment folder: "); isdir(f) || println("  Not found: '$f'"); end
    f
end
isdir(folder) || error("Folder not found: $folder")

data_path  = joinpath(folder, "data_channel_1.csv")
annot_path = joinpath(folder, "annotation_clean.csv")
isfile(data_path) || error("data_channel_1.csv not found in $folder")

println("\n━━━ p0 Benchmark — $(basename(folder)) ━━━\n")
println("  Strategies compared:")
println("  [A] BBO   + default_p0  (bounds midpoint, full search space) ← baseline")
println("  [B] BBO   + narrow bounds ($(Int(BOUND_LO*100))%–$(Int(BOUND_HI*100))% of suggest_p0)")
println("  [C] LBFGS + suggest_p0  (gradient-based, warm start)\n")

# ─────────────────────────────────────────────────────────────────────────────
# Load data + blank
# ─────────────────────────────────────────────────────────────────────────────
data = GrowthData(data_path)
println("  Data: $(size(data.curves,1)) curves × $(size(data.curves,2)) timepoints")

blank_value = 0.0
if isfile(annot_path)
    annot       = [(split(l,",")[1], strip(split(l,",")[2]))
                   for l in readlines(annot_path) if !isempty(strip(l))]
    blank_wells = [w for (w,t) in annot if lowercase(t) == "b"]
    if !isempty(blank_wells)
        blank_idx   = filter(!isnothing, [findfirst(==(w), data.labels) for w in blank_wells])
        blank_value = mean(data.curves[blank_idx, :])
        println("  Blank OD: $(round(blank_value; digits=5))")
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Choose well
# ─────────────────────────────────────────────────────────────────────────────
well = if length(ARGS) >= 2
    strip(ARGS[2])
else
    println("  Wells: $(join(data.labels, "  "))")
    ask("Well: ")
end
well_idx = findfirst(==(well), data.labels)
isnothing(well_idx) && error("Well '$well' not found.")
println("  Well: $well")

raw_od    = max.(data.curves[well_idx, :] .- blank_value, 0.001)
t_stat    = find_stationary_phase(data.times, raw_od)
if t_stat !== nothing
    println("  Stationary phase at t=$(round(t_stat;digits=2))h  " *
            "($(count(data.times .≤ t_stat))/$(length(data.times)) points)")
else
    println("  No stationary phase detected — full time range used")
end

mask      = t_stat !== nothing ? data.times .≤ t_stat : trues(length(data.times))
t_trunc   = data.times[mask]
y_trunc   = raw_od[mask]
well_data = GrowthData(reshape(y_trunc, 1, :), t_trunc, [well])

# ─────────────────────────────────────────────────────────────────────────────
# Load / build fingerprint DB
# ─────────────────────────────────────────────────────────────────────────────
db_path = joinpath(tempdir(), "kinbiont_db_$(basename(folder)).jls")
db = if isfile(db_path)
    candidate = deserialize(db_path)
    if !(candidate isa ModelFingerprintDB) || !isdefined(candidate, :params) || isempty(candidate.params) || !isdefined(candidate, :curve_matrix)
        println("\n  Cached DB outdated — rebuilding …"); nothing
    else
        println("  DB: $(length(unique(candidate.model_names))) models, " *
                "$(length(candidate.model_names)) fingerprints (from cache)")
        candidate
    end
else
    nothing
end

if isnothing(db)
    print("  Building fingerprint DB … ")
    Random.seed!(42)
    db = build_model_fingerprint_db(n_samples=200, tmax=data.times[end],
                                    n_points=60, od_range=(0.001, 0.02))
    serialize(db_path, db)
    println("$(length(unique(db.model_names))) models, $(length(db.model_names)) fingerprints")
end

# ─────────────────────────────────────────────────────────────────────────────
# Choose models
# ─────────────────────────────────────────────────────────────────────────────
model_names = if length(ARGS) >= 3
    String.(split(strip(ARGS[3]), r"[\s,]+"; keepempty=false))
else
    println("\nModels:")
    println("  [1] Recommended top-5 + aHPM")
    println("  [2] Manual list")
    ch = ask("Choice [1/2]: ")
    if ch == "2"
        print("  Names (space/comma separated): ")
        String.(split(strip(readline()), r"[\s,]+"; keepempty=false))
    else
        recs = recommend_models(y_trunc, t_trunc, db; top_k=5)
        "aHPM" in recs || pushfirst!(recs, "aHPM")
        recs
    end
end

valid   = filter(n -> haskey(MODEL_REGISTRY, n), model_names)
invalid = filter(n -> !haskey(MODEL_REGISTRY, n), model_names)
isempty(invalid) || println("  Ignored unknown: $(join(invalid, ", "))")
isempty(valid)   && error("No valid models.")
println("\n  Models : $(join(valid, ", "))")
println("  Runs   : $N_RUNS per strategy (seeds: $(join(SEEDS, ", ")))\n")

# ─────────────────────────────────────────────────────────────────────────────
# Fit options per strategy
# ─────────────────────────────────────────────────────────────────────────────
base_opts = (
    blank_subtraction  = false,
    smooth             = false,
    correct_negatives  = true,
    negative_method    = :thr_correction,
    negative_threshold = 0.001,
    loss               = "RE",
)

opts_A = FitOptions(; base_opts...,
    optimizer  = BBO_adaptive_de_rand_1_bin_radiuslimited(),
    opt_params = (maxiters = 10_000, MaxTime = 30.0),
)

opts_B = FitOptions(; base_opts...,
    optimizer  = BBO_adaptive_de_rand_1_bin_radiuslimited(),
    opt_params = (maxiters = 10_000, MaxTime = 30.0),
)

opts_C = FitOptions(; base_opts...,
    optimizer  = LBFGS(),
    opt_params = (maxiters = 2_000,),
)

# ─────────────────────────────────────────────────────────────────────────────
# Benchmark
# ─────────────────────────────────────────────────────────────────────────────
struct BenchRow
    model    :: String
    strategy :: String
    p0       :: Vector{Float64}
    lb       :: Union{Nothing, Vector{Float64}}
    ub       :: Union{Nothing, Vector{Float64}}
    aicc     :: Float64
    loss     :: Float64
    mean_t   :: Float64
    std_t    :: Float64
end

all_rows = BenchRow[]

function run_strategy(label, model, p0, lb, ub, opts)
    aiccs = Float64[]; losses = Float64[]; times = Float64[]
    for seed in SEEDS
        Random.seed!(seed)
        t0 = time()
        lb_arg = isnothing(lb) ? nothing : [lb]
        ub_arg = isnothing(ub) ? nothing : [ub]
        spec = ModelSpec([model], [p0], lb_arg, ub_arg)
        r = try
            kinbiont_fit(well_data, spec, opts)
        catch
            nothing
        end
        el = time() - t0
        if !isnothing(r)
            push!(aiccs,  r.results[1].best_aic)
            push!(losses, r.results[1].loss)
            push!(times,  el)
        end
    end
    isempty(times) && return nothing
    return (mean(aiccs), mean(losses), mean(times), std(times))
end

for name in valid
    model   = MODEL_REGISTRY[name]
    p0_nn   = suggest_p0(y_trunc, t_trunc, db, name)
    p0_def  = Kinbiont._default_p0(model, well_data)
    p0_used = isempty(p0_nn) ? p0_def : p0_nn

    println("── $name " * "─"^max(1, 55 - length(name)))

    # Strategy A: BBO + default_p0, no custom bounds
    print("  [A] BBO   + default_p0  …\r")
    res = run_strategy("A", model, p0_def, nothing, nothing, opts_A)
    if !isnothing(res)
        aicc, loss, mt, st = res
        println("  [A] BBO   + default_p0  " *
                "  AICc $(rpad(round(aicc;digits=1),9))" *
                "  loss $(rpad(round(loss;digits=5),10))" *
                "  $(round(mt;digits=1))s ± $(round(st;digits=1))s")
        push!(all_rows, BenchRow(name,"A: BBO default",p0_def,nothing,nothing,aicc,loss,mt,st))
    else
        println("  [A] BBO   + default_p0     FAILED")
    end

    # Strategy B: BBO + narrow bounds around suggest_p0
    lb_b = max.(p0_used .* BOUND_LO, 1e-6)
    ub_b = p0_used .* BOUND_HI
    # ensure lb < ub
    for i in eachindex(lb_b); lb_b[i] = min(lb_b[i], ub_b[i] * 0.5); end
    print("  [B] BBO   + narrow bounds …\r")
    res = run_strategy("B", model, p0_used, lb_b, ub_b, opts_B)
    if !isnothing(res)
        aicc, loss, mt, st = res
        println("  [B] BBO   + narrow bounds " *
                "  AICc $(rpad(round(aicc;digits=1),9))" *
                "  loss $(rpad(round(loss;digits=5),10))" *
                "  $(round(mt;digits=1))s ± $(round(st;digits=1))s")
        push!(all_rows, BenchRow(name,"B: BBO narrow",p0_used,lb_b,ub_b,aicc,loss,mt,st))
    else
        println("  [B] BBO   + narrow bounds  FAILED")
    end

    # Strategy C: LBFGS + suggest_p0
    print("  [C] LBFGS + suggest_p0   …\r")
    res = run_strategy("C", model, p0_used, nothing, nothing, opts_C)
    if !isnothing(res)
        aicc, loss, mt, st = res
        println("  [C] LBFGS + suggest_p0   " *
                "  AICc $(rpad(round(aicc;digits=1),9))" *
                "  loss $(rpad(round(loss;digits=5),10))" *
                "  $(round(mt;digits=1))s ± $(round(st;digits=1))s")
        push!(all_rows, BenchRow(name,"C: LBFGS",p0_used,nothing,nothing,aicc,loss,mt,st))
    else
        println("  [C] LBFGS + suggest_p0    FAILED (try BBO if model ODE is non-smooth)")
    end

    # Δ vs baseline A
    row_a = findfirst(r -> r.model == name && r.strategy == "A: BBO default", all_rows)
    for (tag, strat) in [("B", "B: BBO narrow"), ("C", "C: LBFGS")]
        row_x = findfirst(r -> r.model == name && r.strategy == strat, all_rows)
        if !isnothing(row_a) && !isnothing(row_x)
            Δa = all_rows[row_x].aicc    - all_rows[row_a].aicc
            Δt = all_rows[row_x].mean_t  - all_rows[row_a].mean_t
            sa = Δa < -0.5 ? "▼ better" : (Δa > 0.5 ? "▲ worse " : "≈ same  ")
            st = Δt < 0    ? "▼"        : "▲"
            println("      Δ vs A  [$tag]: AICc $sa $(rpad(round(Δa;digits=1),7))  " *
                    "time $st $(round(Δt;digits=1))s")
        end
    end
    println()
end

# ─────────────────────────────────────────────────────────────────────────────
# Summary
# ─────────────────────────────────────────────────────────────────────────────
println("━━━ Summary ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
println(rpad("Model",32) * rpad("Strategy",26) *
        rpad("AICc",10) * rpad("Loss",11) * "Time (s)")
println("─"^88)
for name in valid
    rows = filter(r -> r.model == name, all_rows)
    for r in rows
        println(rpad(r.model,32) * rpad(r.strategy,26) *
                rpad(round(r.aicc;digits=1),10) *
                rpad(round(r.loss;digits=5),11) *
                "$(round(r.mean_t;digits=1)) ± $(round(r.std_t;digits=1))")
    end
    println()
end

# ─────────────────────────────────────────────────────────────────────────────
# Save CSV
# ─────────────────────────────────────────────────────────────────────────────
out_dir  = joinpath(folder, "kinbiont_results"); mkpath(out_dir)
csv_path = joinpath(out_dir, "benchmark_p0_$(well).csv")
open(csv_path, "w") do f
    println(f, "well,model,strategy,aicc,loss,mean_time_s,std_time_s,p0")
    for r in all_rows
        println(f, "$well,$(r.model),$(r.strategy),$(r.aicc),$(r.loss)," *
                   "$(r.mean_t),$(r.std_t),\"$(join(round.(r.p0;digits=6), ";"))\"")
    end
end
println("Results saved → $csv_path")
println("\n━━━ Done ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
