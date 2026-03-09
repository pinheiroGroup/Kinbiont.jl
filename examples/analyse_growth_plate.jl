using Kinbiont
using Random
using Serialization
using Statistics: mean
using Logging

global_logger(ConsoleLogger(stderr, Logging.Error))

const DONE = "  [✓]"
const RUN  = "  [ ]"
const WARN = "  [!]"

# ── Stationary phase detection ────────────────────────────────────────────────
# Port of find_stationary_phase from aHPM-web-interface/function_for_fitting.jl.
# Returns the cutoff time (Float64) before stationary phase, or nothing if not found.
#
# Algorithm:
#   1. Ignore points with OD < thr_od (below noise floor)
#   2. Compute specific growth rate via rolling OLS on log(OD) (pt_smooth_derivative points)
#   3. Threshold = percentile_thr × max(growth rate)
#   4. Starting from the growth-rate peak, scan forward for the first window of
#      win_size consecutive points all below the threshold → stationary phase start
function find_stationary_phase(
    times::Vector{Float64},
    od::Vector{Float64};
    percentile_thr::Float64 = 0.05,
    pt_smooth_derivative::Int = 7,
    win_size::Int = 5,
    thr_od::Float64 = 0.001,
)::Union{Float64, Nothing}

    data_mat = Matrix(transpose(hcat(times, od)))   # 2 × n (Kinbiont legacy format)

    index_od = findall(data_mat[2, :] .> thr_od)
    isempty(index_od) && return nothing

    data_t = data_mat[:, index_od]
    length(data_t[1, :]) <= pt_smooth_derivative && return nothing

    sgr = Kinbiont.specific_gr_evaluation(data_t, pt_smooth_derivative)
    isempty(sgr) && return nothing

    max_sgr = maximum(sgr)
    max_sgr <= 0 && return nothing

    thr           = max_sgr * percentile_thr
    maximum_index = argmax(sgr)

    for i in maximum_index:(length(sgr) - win_size)
        if all(sgr[i:(i + win_size - 1)] .< thr)
            cutoff_idx = i + index_od[1] - 1   # map back to original indices
            return times[min(cutoff_idx, length(times))]
        end
    end
    return nothing
end

# ── 0. Resolve folder path ────────────────────────────────────────────────────
folder = if !isempty(ARGS)
    f = strip(ARGS[1])
    isdir(f) || error("Folder not found: $f")
    f
else
    f = ""
    while !isdir(f)
        print("Enter experiment folder path: ")
        f = strip(readline())
        isdir(f) || println("  Folder not found: '$f' — please try again.")
    end
    f
end

data_path  = joinpath(folder, "data_channel_1.csv")
annot_path = joinpath(folder, "annotation_clean.csv")
isfile(data_path)  || error("data_channel_1.csv not found in $folder")
isfile(annot_path) || error("annotation_clean.csv not found in $folder")

println("\n━━━ KinBiont Analysis — $(basename(folder)) ━━━\n")

# ── 1. Parse annotation ───────────────────────────────────────────────────────
print("$RUN Parsing annotation …\r")
annot = [(split(line, ",")[1], split(line, ",")[2])
         for line in readlines(annot_path) if !isempty(strip(line))]

blank_wells    = [well for (well, type) in annot if lowercase(type) == "b"]
excluded_wells = [well for (well, type) in annot if lowercase(type) == "x"]
skip_wells     = union(blank_wells, excluded_wells)

println("$DONE Annotation parsed")
println("       Blank wells   : ", isempty(blank_wells)    ? "none" : join(blank_wells, ", "))
println("       Excluded wells: ", isempty(excluded_wells) ? "none" : join(excluded_wells, ", "))

# ── 2. Load data ──────────────────────────────────────────────────────────────
print("$RUN Loading data …\r")
data = GrowthData(data_path)
println("$DONE Data loaded  — $(size(data.curves,1)) curves × $(size(data.curves,2)) timepoints, ",
        "$(round(data.times[1];digits=2))–$(round(data.times[end];digits=2)) h")

if isempty(blank_wells)
    blank_value = 0.0
    do_blank    = false
    println("       No blank wells found — blank subtraction disabled.")
else
    blank_idx   = filter(!isnothing, [findfirst(==(w), data.labels) for w in blank_wells])
    blank_value = mean(data.curves[blank_idx, :])
    do_blank    = true
    println("       Blank OD (mean): $(round(blank_value; digits=5))")
end

# ── 3. Select wells ───────────────────────────────────────────────────────────
sample_wells = [l for l in data.labels if l ∉ skip_wells]
println("\n       Available sample wells ($(length(sample_wells))): ", join(sample_wells, "  "))

println("""
\nFit which curves?
  [1] All sample wells ($(length(sample_wells)))
  [2] Select wells manually
""")
print("Choice [1/2]: ")
choice = strip(readline())

wells_to_fit = if choice == "2"
    println("Enter well names separated by spaces or commas (e.g. A1 A2 B3):")
    print("> ")
    raw      = strip(readline())
    selected = String.(split(raw, r"[\s,]+"; keepempty=false))
    valid    = filter(w -> w in sample_wells, selected)
    invalid  = filter(w -> w ∉ sample_wells, selected)
    isempty(invalid) || println("       Ignored unknown/excluded wells: $(join(invalid, ", "))")
    isempty(valid)   && error("No valid wells selected.")
    println("       Selected: ", join(valid, ", "))
    valid
else
    println("       Fitting all $(length(sample_wells)) sample wells.")
    sample_wells
end

fit_indices = [findfirst(==(w), data.labels) for w in wells_to_fit]
fit_data    = GrowthData(data.curves[fit_indices, :], data.times, wells_to_fit)

# ── 4. Build (or load) fingerprint DB ─────────────────────────────────────────
db_path = joinpath(tempdir(), "kinbiont_db_$(basename(folder)).jls")

db = if isfile(db_path)
    print("$RUN Loading cached fingerprint DB …\r")
    candidate = deserialize(db_path)
    if !(candidate isa ModelFingerprintDB) || !isdefined(candidate, :params) || isempty(candidate.params)
        println("$WARN Cached DB is outdated (no param storage) — rebuilding …")
        nothing
    else
        println("$DONE Fingerprint DB loaded  — $(length(unique(candidate.model_names))) models, $(length(candidate.model_names)) fingerprints  (from cache)")
        candidate
    end
else
    nothing
end

db = if isnothing(db)
    print("$RUN Building fingerprint DB  (sampling all models — this takes ~1–2 min) …\r")
    Random.seed!(42)
    db = build_model_fingerprint_db(
        n_samples = 200,
        tmax      = data.times[end],
        n_points  = 60,
        od_range  = (0.001, 0.02),
    )
    serialize(db_path, db)
    println("$DONE Fingerprint DB built   — $(length(unique(db.model_names))) models, $(length(db.model_names)) fingerprints  (cached to $db_path)")
    db
end

# ── 5. Select models ───────────────────────────────────────────────────────────
println("""
\nWhich models to fit?
  [1] Recommended  — top-5 per curve via fingerprint DB + aHPM always included
  [2] Manual list  — enter model names (DB still used for initial parameter guessing)
""")
print("Choice [1/2]: ")
model_choice = strip(readline())

all_recs = if model_choice == "2"
    println("Enter model name(s) separated by spaces or commas")
    println("  (e.g.  aHPM   or   aHPM, logistic, baranyi_roberts)")
    println("  Available: $(join(sort(collect(keys(MODEL_REGISTRY))), ", "))\n")
    print("> ")
    raw_models   = strip(readline())
    entered      = String.(split(raw_models, r"[\s,]+"; keepempty=false))
    valid_models = filter(n -> haskey(MODEL_REGISTRY, n), entered)
    bad_models   = filter(n -> !haskey(MODEL_REGISTRY, n), entered)
    isempty(bad_models) || println("$WARN  Unknown model(s) ignored: $(join(bad_models, ", "))")
    isempty(valid_models) && error("No valid model names entered.")
    println("$DONE Using manual model list: $(join(valid_models, ", "))")
    valid_models
else
    print("$RUN Recommending models …\r")
    recs = unique(vcat([
        recommend_models(fit_data.curves[i, :], fit_data.times, db; top_k=5)
        for i in axes(fit_data.curves, 1)
    ]...))
    # aHPM is always included regardless of recommendation rank.
    "aHPM" in recs || pushfirst!(recs, "aHPM")
    println("$DONE Model recommendation done — $(length(recs)) model(s): $(join(recs, ", "))")
    recs
end

# ── 6. Fit options ─────────────────────────────────────────────────────────────
# Used when fitting the truncated (already blank-subtracted) data.
# Blank subtraction + smoothing happen before truncation (see below), so disabled here.
fit_opts = FitOptions(
    blank_subtraction  = false,
    smooth             = false,
    correct_negatives  = true,
    negative_method    = :thr_correction,
    negative_threshold = 0.001,
    loss               = "RE",
    opt_params         = (maxiters = 10_000, MaxTime = 60.0),
)

# ── 7. Per-well fitting with stationary phase detection ───────────────────────
println("\n       $(length(wells_to_fit)) well(s) × $(length(all_recs)) model(s) — p0 from fingerprint DB nearest neighbour\n")

# Collected per-well results
best_curve_results = CurveFitResult[]
stationary_times   = Union{Float64,Nothing}[]
preprocessed_obs   = Tuple{Vector{Float64},Vector{Float64}}[]  # (times, od) per well, full range

for (wi, well) in enumerate(wells_to_fit)
    println("  ── $well ──────────────────────────────────────────────────────")

    raw_times = fit_data.times
    raw_od    = fit_data.curves[wi, :]

    # --- Manual blank subtraction + negative clipping for detection/observation ---
    od_corrected = max.(raw_od .- blank_value, 0.001)

    push!(preprocessed_obs, (raw_times, od_corrected))

    # --- Stationary phase detection ---
    t_stat = find_stationary_phase(raw_times, od_corrected;
                 percentile_thr    = 0.05,
                 pt_smooth_derivative = 7,
                 win_size          = 5,
                 thr_od            = 0.001)
    push!(stationary_times, t_stat)

    if t_stat !== nothing
        println("       Stationary phase detected at t = $(round(t_stat; digits=2)) h  ",
                "(using $(count(raw_times .≤ t_stat))/$(length(raw_times)) points)")
    else
        println("$WARN  No stationary phase found — using full time range")
    end

    # --- Truncate to growth phase only ---
    mask    = t_stat !== nothing ? raw_times .≤ t_stat : trues(length(raw_times))
    t_trunc = raw_times[mask]
    y_trunc = od_corrected[mask]
    well_data = GrowthData(reshape(y_trunc, 1, :), t_trunc, [well])

    # --- Fit each recommended model ---
    well_model_results = CurveFitResult[]
    for (mi, name) in enumerate(all_recs)
        pad = rpad("  ($mi/$(length(all_recs))) $name", 55)
        print("$RUN $pad\r")
        model  = MODEL_REGISTRY[name]
        p0_nn  = suggest_p0(y_trunc, t_trunc, db, name)
        p0     = isempty(p0_nn) ? Kinbiont._default_p0(model, well_data) : p0_nn
        spec_i = ModelSpec([model], [p0])
        t0     = time()
        r      = kinbiont_fit(well_data, spec_i, fit_opts)
        elapsed = round(time() - t0; digits=1)
        println("$DONE $pad  $(elapsed)s  AICc=$(round(r.results[1].best_aic; digits=1))")
        push!(well_model_results, r.results[1])
    end

    best = well_model_results[argmin([r.best_aic for r in well_model_results])]
    push!(best_curve_results, best)
    println("       → Best: $(best.best_model.name)  AICc=$(round(best.best_aic; digits=1))")
    println()
end

# ── 8. Summary ────────────────────────────────────────────────────────────────
println("\n━━━ Results ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
println(rpad("Well", 6), rpad("Stat. phase (h)", 18), rpad("Best model", 38),
        rpad("AICc", 9), "Loss")
println("─"^80)
for (i, well) in enumerate(wells_to_fit)
    r      = best_curve_results[i]
    t_stat = stationary_times[i]
    println(
        rpad(well, 6),
        rpad(t_stat !== nothing ? string(round(t_stat; digits=1)) : "not found", 18),
        rpad(r.best_model.name, 38),
        rpad(round(r.best_aic; digits=1), 9),
        round(r.loss; digits=5),
    )
end

freq = Dict{String,Int}()
for r in best_curve_results
    freq[r.best_model.name] = get(freq, r.best_model.name, 0) + 1
end
println("\n─── Model selection frequency ──────────────────────────────────────")
for (name, count) in sort(collect(freq); by = x -> -x[2])
    println("  ", rpad(name, 38), "× $count")
end

# ── 9. Save tabular results ───────────────────────────────────────────────────
print("$RUN Saving tabular results …\r")
out_dir = joinpath(folder, "kinbiont_results")
mkpath(out_dir)

# Write summary CSV manually (GrowthFitResults expected; build a minimal one per-well)
open(joinpath(out_dir, "summary.csv"), "w") do f
    println(f, "well,stationary_phase_h,best_model,best_aicc,best_loss")
    for (i, well) in enumerate(wells_to_fit)
        r      = best_curve_results[i]
        t_stat = stationary_times[i]
        println(f, "$well,$(t_stat !== nothing ? t_stat : ""),$(r.best_model.name),$(r.best_aic),$(r.loss)")
    end
end
println("$DONE Tabular results saved to $out_dir/summary.csv")

# ── 10. Save plot data CSVs ───────────────────────────────────────────────────
print("$RUN Saving plot data …\r")
plot_data_dir = joinpath(out_dir, "plot_data")
mkpath(plot_data_dir)

for (i, well) in enumerate(wells_to_fit)
    r      = best_curve_results[i]
    t_stat = stationary_times[i]
    obs_t, obs_y = preprocessed_obs[i]

    # Full observed curve
    open(joinpath(plot_data_dir, "$(well)_obs.csv"), "w") do f
        println(f, "time,od")
        for (t, y) in zip(obs_t, obs_y); println(f, "$t,$y"); end
    end

    # Best model fit (on truncated range) + stationary phase time
    open(joinpath(plot_data_dir, "$(well)_best.csv"), "w") do f
        println(f, "time,od,model,aicc,loss,t_stat")
        t_s = t_stat !== nothing ? t_stat : ""
        for (t, y) in zip(r.times, r.fitted_curve)
            println(f, "$t,$y,$(r.best_model.name),$(r.best_aic),$(r.loss),$t_s")
        end
    end
end
println("$DONE Plot data saved to $plot_data_dir/")

# ── 11. Write and run plot script in a subprocess ─────────────────────────────
plot_script = joinpath(out_dir, "make_plots.jl")
open(plot_script, "w") do f
    print(f, raw"""
using Plots
gr()

out_dir       = ARGS[1]
plot_data_dir = joinpath(out_dir, "plot_data")

obs_files = filter(f -> endswith(f, "_obs.csv"), readdir(plot_data_dir))
wells = [replace(basename(f), "_obs.csv" => "") for f in obs_files]

function read_csv(path)
    lines  = readlines(path)
    header = split(lines[1], ",")
    cols   = Dict(h => Float64[] for h in header)
    strs   = Dict(h => String[]  for h in header)
    for line in lines[2:end]
        vals = split(line, ",")
        for (h, v) in zip(header, vals)
            n = tryparse(Float64, strip(v))
            if n !== nothing
                push!(cols[h], n)
            else
                push!(strs[h], strip(v))
            end
        end
    end
    return cols, strs
end

for well in wells
    obs,  _      = read_csv(joinpath(plot_data_dir, well * "_obs.csv"))
    best, best_s = read_csv(joinpath(plot_data_dir, well * "_best.csv"))

    best_name = isempty(best_s["model"]) ? "best" : best_s["model"][1]
    best_aicc = isempty(best["aicc"])    ? NaN    : best["aicc"][1]
    t_stat    = isempty(best["t_stat"])  ? nothing : best["t_stat"][1]

    fig = plot(obs["time"], obs["od"];
        seriestype = :scatter, label = "Data (full range)",
        color = :gray, markersize = 2, markerstrokewidth = 0, alpha = 0.6,
        xlabel = "Time (h)", ylabel = "OD",
        title = "$well — $best_name", legend = :bottomright, titlefontsize = 10,
        size = (700, 420), dpi = 150,
        left_margin = 5Plots.mm, bottom_margin = 5Plots.mm)
    plot!(fig, best["time"], best["od"];
        label = "$best_name  AICc=$(round(best_aicc; digits=1))",
        color = :steelblue, linewidth = 2)
    if t_stat !== nothing
        vline!(fig, [t_stat];
            label = "Stationary phase (t=$(round(t_stat; digits=1)) h)",
            color = :black, linestyle = :dot, linewidth = 1.5)
    end

    out_path = joinpath(out_dir, "fit_$(well).png")
    savefig(fig, out_path)
    println("  [✓] $well → $out_path")
end
""")
end

println()
print("$RUN Generating plots (spawning clean Julia process) …\r")
plot_proc = run(ignorestatus(`julia $plot_script $out_dir`))
if plot_proc.exitcode == 0
    println("$DONE Plots saved to $out_dir/fit_<well>.png")
else
    println("$WARN Plot generation failed (exit $(plot_proc.exitcode)) — run manually:")
    println("      julia $plot_script $out_dir")
end

println("\n━━━ Done ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
