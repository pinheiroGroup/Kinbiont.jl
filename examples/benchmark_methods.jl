# benchmark_methods.jl
#
# Compares two model-recommendation + p0-guessing backends:
#
#   :features — 15-element shape-feature vector (cosine similarity)
#   :curves   — min-max normalised full curve (Pearson correlation, scale-invariant)
#
# For a chosen plate + well the script runs the same set of models under both
# backends and reports, for each model:
#   • p0 chosen by each method
#   • fitting time (N_RUNS seeds for reliability)
#   • best AICc reached
#   • Δ AICc  and  Δ time  vs :features baseline
#
# Usage:
#   julia --project examples/benchmark_methods.jl
#   julia --project examples/benchmark_methods.jl <folder>
#   julia --project examples/benchmark_methods.jl <folder> <well>
#   julia --project examples/benchmark_methods.jl <folder> <well> <top_k>
#   julia --project examples/benchmark_methods.jl <folder> <well> <top_k> model1,model2,…

using Kinbiont
using Random
using Serialization
using Statistics: mean, std
using Logging

global_logger(ConsoleLogger(stderr, Logging.Error))

const N_RUNS = 3
const SEEDS  = [42, 43, 44]
const TOP_K  = 5

# ── helpers ───────────────────────────────────────────────────────────────────
function find_stationary_phase(times, od;
        percentile_thr=0.05, pt_smooth_derivative=7, win_size=5, thr_od=0.001)
    data_mat = Matrix(transpose(hcat(times, od)))
    idx_od   = findall(data_mat[2, :] .> thr_od)
    isempty(idx_od) && return nothing
    data_t = data_mat[:, idx_od]
    length(data_t[1, :]) <= pt_smooth_derivative && return nothing
    sgr = Kinbiont.specific_gr_evaluation(data_t, pt_smooth_derivative)
    isempty(sgr) && return nothing
    max_sgr = maximum(sgr); max_sgr <= 0 && return nothing
    thr = max_sgr * percentile_thr; max_i = argmax(sgr)
    for i in max_i:(length(sgr)-win_size)
        all(sgr[i:(i+win_size-1)] .< thr) &&
            return times[min(i + idx_od[1] - 1, length(times))]
    end
    nothing
end

function run_fit(well_data, model, p0, opts)
    aiccs = Float64[]; times = Float64[]
    for seed in SEEDS
        Random.seed!(seed)
        t0 = time()
        r  = try kinbiont_fit(well_data, ModelSpec([model], [p0]), opts) catch; nothing end
        el = time() - t0
        isnothing(r) || (push!(aiccs, r.results[1].best_aic); push!(times, el))
    end
    isempty(times) && return nothing
    (mean(aiccs), mean(times), std(times))
end

# ── resolve inputs ────────────────────────────────────────────────────────────
folder = if !isempty(ARGS)
    f = strip(ARGS[1]); isdir(f) || error("Folder not found: $f"); f
else
    f = ""
    while !isdir(f); print("Experiment folder: "); f = strip(readline()); isdir(f) || println("  Not found."); end; f
end

data_path  = joinpath(folder, "data_channel_1.csv")
annot_path = joinpath(folder, "annotation_clean.csv")
isfile(data_path) || error("data_channel_1.csv not found")

println("\n━━━ Method Benchmark — $(basename(folder)) ━━━\n")
println("  Comparing p0 source:")
println("  [:features] 15-element shape-feature cosine similarity")
println("  [:curves]   min-max normalised Pearson correlation (scale-invariant)\n")

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

# ── choose well ───────────────────────────────────────────────────────────────
well = if length(ARGS) >= 2
    strip(ARGS[2])
else
    println("  Wells: $(join(data.labels, "  "))")
    print("Well: "); strip(readline())
end
well_idx = findfirst(==(well), data.labels)
isnothing(well_idx) && error("Well '$well' not found.")

raw_od    = max.(data.curves[well_idx, :] .- blank_value, 0.001)
t_stat    = find_stationary_phase(data.times, raw_od)
mask      = t_stat !== nothing ? data.times .≤ t_stat : trues(length(data.times))
t_trunc   = data.times[mask]
y_trunc   = raw_od[mask]
well_data = GrowthData(reshape(y_trunc, 1, :), t_trunc, [well])

println("  Well: $well  " *
        (t_stat !== nothing ? "stationary at $(round(t_stat;digits=2))h" : "full time range"))

# ── load / build DB ───────────────────────────────────────────────────────────
db_path = joinpath(tempdir(), "kinbiont_db_$(basename(folder)).jls")
db = if isfile(db_path)
    candidate = deserialize(db_path)
    if !(candidate isa ModelFingerprintDB) ||
       !isdefined(candidate, :params)       || isempty(candidate.params) ||
       !isdefined(candidate, :curve_matrix)
        println("  Cached DB outdated — rebuilding …"); nothing
    else
        println("  DB: $(length(unique(candidate.model_names))) models, " *
                "$(length(candidate.model_names)) fingerprints  (from cache)")
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

# ── score_models diagnostic ──────────────────────────────────────────────────
println("\n  ── Similarity scores for this well ──────────────────────────────")
for method in [:features, :curves]
    scores = score_models(y_trunc, t_trunc, db; method)
    top10  = scores[1:min(10, length(scores))]
    ahpm_rank = findfirst(p -> p.first == "aHPM", scores)
    println("\n  [$method]  top-10:")
    for (rank, (name, score)) in enumerate(top10)
        marker = name == "aHPM" ? " ◄" : ""
        println("    $(lpad(rank,2)). $(rpad(name,32))  $(round(score;digits=4))$marker")
    end
    if !isnothing(ahpm_rank) && ahpm_rank > 10
        println("    …  aHPM at rank $ahpm_rank / $(length(scores))  " *
                "score=$(round(scores[ahpm_rank].second;digits=4))")
    end
end

# ── override top_k from ARGS[3] if provided ──────────────────────────────────
top_k = if length(ARGS) >= 3
    v = tryparse(Int, strip(ARGS[3]))
    isnothing(v) ? TOP_K : v
else
    TOP_K
end

# ── choose models ─────────────────────────────────────────────────────────────
model_names = if length(ARGS) >= 4
    # ARGS[4] = comma/space separated model names
    String.(split(strip(ARGS[4]), r"[\s,]+"; keepempty=false))
else
    println("\n  Which models to benchmark?")
    println("  [1] Top-$top_k recommendations (union of both methods) + aHPM")
    println("  [2] Manual list")
    print("Choice [1/2]: "); ch = strip(readline())
    if ch == "2"
        print("  Names (space/comma): "); String.(split(strip(readline()), r"[\s,]+"; keepempty=false))
    else
        recs_f = recommend_models(y_trunc, t_trunc, db; top_k, method=:features)
        recs_c = recommend_models(y_trunc, t_trunc, db; top_k, method=:curves)
        recs   = unique(vcat(recs_f, recs_c))
        "aHPM" in recs || push!(recs, "aHPM")
        recs
    end
end

valid   = filter(n -> haskey(MODEL_REGISTRY, n), model_names)
invalid = filter(n -> !haskey(MODEL_REGISTRY, n), model_names)
isempty(invalid) || println("  Ignored unknown: $(join(invalid, ", "))")
isempty(valid)   && error("No valid models.")
println("\n  Models  : $(join(valid, ", "))")
println("  Runs    : $N_RUNS per method (seeds: $(join(SEEDS, ", ")))\n")

# ── fit options ───────────────────────────────────────────────────────────────
fit_opts = FitOptions(
    blank_subtraction  = false,
    smooth             = false,
    correct_negatives  = true,
    negative_method    = :thr_correction,
    negative_threshold = 0.001,
    loss               = "RE",
    opt_params         = (maxiters = 10_000, MaxTime = 30.0),
)

# ── benchmark ─────────────────────────────────────────────────────────────────
struct BenchRow
    model     :: String
    method    :: Symbol
    aicc      :: Float64
    mean_t    :: Float64
    std_t     :: Float64
end

all_rows = BenchRow[]

for name in valid
    model = MODEL_REGISTRY[name]
    println("── $name " * "─"^max(1, 55 - length(name)))

    for method in [:features, :curves]
        p0 = suggest_p0(y_trunc, t_trunc, db, name; method)
        if isempty(p0)
            p0 = Kinbiont._default_p0(model, well_data)
        end

        tag = "[:$method]"
        print("  $tag $(rpad("",12)) fitting …\r")
        res = run_fit(well_data, model, p0, fit_opts)

        if !isnothing(res)
            aicc, mt, st = res
            println("  $tag $(rpad("",12)) AICc $(rpad(round(aicc;digits=1),9))  " *
                    "$(round(mt;digits=1))s ± $(round(st;digits=1))s  " *
                    "p0=[$(join(round.(p0;digits=3), ", "))]")
            push!(all_rows, BenchRow(name, method, aicc, mt, st))
        else
            println("  $tag  FAILED")
        end
    end

    # Δ :curves vs :features
    row_f = findfirst(r -> r.model == name && r.method == :features, all_rows)
    row_c = findfirst(r -> r.model == name && r.method == :curves,   all_rows)
    if !isnothing(row_f) && !isnothing(row_c)
        Δa = all_rows[row_c].aicc   - all_rows[row_f].aicc
        Δt = all_rows[row_c].mean_t - all_rows[row_f].mean_t
        sa = abs(Δa) < 2 ? "≈ same  " : (Δa < 0 ? "▼ better" : "▲ worse ")
        println("  Δ [:curves − :features]  AICc $sa $(rpad(round(Δa;digits=1),7))  " *
                "time $(Δt < 0 ? "▼" : "▲") $(round(Δt;digits=1))s")
    end
    println()
end

# ── summary ───────────────────────────────────────────────────────────────────
println("━━━ Summary ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
println(rpad("Model",32) * rpad("Method",12) * rpad("AICc",10) * "Time (s)")
println("─"^68)
for name in valid
    rows = filter(r -> r.model == name, all_rows)
    for r in rows
        println(rpad(r.model,32) * rpad(":$(r.method)",12) *
                rpad(round(r.aicc;digits=1),10) *
                "$(round(r.mean_t;digits=1)) ± $(round(r.std_t;digits=1))")
    end
    println()
end

curves_wins  = count(n -> begin
    f = findfirst(r -> r.model==n && r.method==:features, all_rows)
    c = findfirst(r -> r.model==n && r.method==:curves,   all_rows)
    !isnothing(f) && !isnothing(c) && all_rows[c].aicc < all_rows[f].aicc - 2
end, valid)
feat_wins    = count(n -> begin
    f = findfirst(r -> r.model==n && r.method==:features, all_rows)
    c = findfirst(r -> r.model==n && r.method==:curves,   all_rows)
    !isnothing(f) && !isnothing(c) && all_rows[f].aicc < all_rows[c].aicc - 2
end, valid)
ties = length(valid) - curves_wins - feat_wins
println("  :curves  better (ΔAICC > 2) for $curves_wins / $(length(valid)) models")
println("  :features better (ΔAICC > 2) for $feat_wins / $(length(valid)) models")
println("  equivalent (|ΔAICC| ≤ 2)    for $ties / $(length(valid)) models")

# ── save CSV ──────────────────────────────────────────────────────────────────
out_dir  = joinpath(folder, "kinbiont_results"); mkpath(out_dir)
csv_path = joinpath(out_dir, "benchmark_methods_$(well).csv")
open(csv_path, "w") do f
    println(f, "well,model,method,aicc,mean_time_s,std_time_s")
    for r in all_rows
        println(f, "$well,$(r.model),:$(r.method),$(r.aicc),$(r.mean_t),$(r.std_t)")
    end
end
println("\nResults saved → $csv_path")
println("\n━━━ Done ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
