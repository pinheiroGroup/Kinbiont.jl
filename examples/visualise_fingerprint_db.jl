# visualise_fingerprint_db.jl
#
# Re-simulates a handful of growth curves per model using the same parameter
# bounds as build_model_fingerprint_db, then saves one PNG per model into
# an output folder.
#
# Usage:
#   julia --project examples/visualise_fingerprint_db.jl \
#         /tmp/kinbiont_db_LG177.jls \
#         /tmp/db_curves/
#
# If the DB path is omitted, looks for /tmp/kinbiont_db_*.jls and picks the
# most recently modified one.

using Kinbiont
using Serialization
using Random

# ── resolve paths ─────────────────────────────────────────────────────────────
db_path = if length(ARGS) >= 1
    ARGS[1]
else
    candidates = filter(f -> startswith(basename(f), "kinbiont_db_") &&
                              endswith(f, ".jls"),
                        readdir(tempdir(); join=true))
    isempty(candidates) && error("No kinbiont_db_*.jls found in /tmp/. Build one first.")
    sort(candidates; by=mtime)[end]
end

out_dir = length(ARGS) >= 2 ? ARGS[2] : joinpath(dirname(db_path), "db_curves")
mkpath(out_dir)

println("Loading DB: $db_path")
db = deserialize(db_path)
println("  $(length(unique(db.model_names))) models, $(length(db.model_names)) fingerprint rows\n")

# ── simulation settings (match what the DB was built with) ───────────────────
const TMAX     = 24.0
const N_POINTS = 60
const N_CURVES = 8       # curves to draw per model
const OD_RANGE = (0.001, 0.02)

Random.seed!(1)
t_grid = collect(range(0.0, TMAX; length=N_POINTS))

# ── write CSV data per model, then plot in a subprocess ──────────────────────
csv_dir = joinpath(out_dir, "csv")
mkpath(csv_dir)

unique_models = unique(db.model_names)
println("Simulating $N_CURVES curves for each of $(length(unique_models)) models …\n")

for model_name in unique_models
    # find the model object from the DB
    idx   = findfirst(==(model_name), db.model_names)
    model = db.model_objects[idx]

    lb, ub = try
        Kinbiont._model_param_bounds(model, TMAX, N_POINTS)
    catch
        println("  [!] $model_name — could not get param bounds, skipping")
        continue
    end
    isempty(lb) && continue

    samples = Kinbiont._lhs_sample(lb, ub, N_CURVES; rng=Random.GLOBAL_RNG)

    csv_path = joinpath(csv_dir, "$(model_name).csv")
    open(csv_path, "w") do f
        println(f, "curve,time,od")
        for s in 1:N_CURVES
            p   = samples[s, :]
            od0 = OD_RANGE[1] + rand() * (OD_RANGE[2] - OD_RANGE[1])

            result = if model isa Kinbiont.ODEModel && model.n_eq >= 2
                ic_alpha = 0.05 + rand() * 0.90
                Kinbiont._simulate_model(model, p, t_grid, od0, TMAX, N_POINTS; ic_alpha)
            else
                Kinbiont._simulate_model(model, p, t_grid, od0, TMAX, N_POINTS)
            end

            result === nothing && continue
            t_sim, y = result
            for (t, y_val) in zip(t_sim, y)
                println(f, "$s,$t,$y_val")
            end
        end
    end
    println("  [✓] $model_name")
end

println("\nCSVs written to $csv_dir/")

# ── plot script (runs in a clean subprocess — no --project) ──────────────────
plot_script = joinpath(out_dir, "make_plots.jl")
open(plot_script, "w") do f
    print(f, """
using Plots
gr()

csv_dir = ARGS[1]
out_dir = ARGS[2]
mkpath(out_dir)

for csv_file in readdir(csv_dir; join=true)
    endswith(csv_file, ".csv") || continue
    model_name = replace(basename(csv_file), ".csv" => "")

    lines  = readlines(csv_file)
    header = split(lines[1], ",")
    curves = Dict{String, Tuple{Vector{Float64}, Vector{Float64}}}()

    for line in lines[2:end]
        isempty(strip(line)) && continue
        parts = split(line, ",")
        id    = parts[1]
        t     = parse(Float64, parts[2])
        y     = parse(Float64, parts[3])
        if !haskey(curves, id)
            curves[id] = (Float64[], Float64[])
        end
        push!(curves[id][1], t)
        push!(curves[id][2], y)
    end

    isempty(curves) && continue

    p = plot(title=model_name, xlabel="Time (h)", ylabel="OD",
             legend=false, titlefontsize=10, size=(600,380), dpi=130)
    for (_, (t, y)) in curves
        plot!(p, t, y; linewidth=1.5, alpha=0.75)
    end

    savefig(p, joinpath(out_dir, "\$(model_name).png"))
    println("  [✓] \$(model_name)")
end
""")
end

print("Generating plots …")
proc = run(ignorestatus(`julia $plot_script $csv_dir $out_dir`))
if proc.exitcode == 0
    println(" done.")
    println("Plots saved to $out_dir/")
else
    println(" failed (exit $(proc.exitcode)). Run manually:")
    println("  julia $plot_script $csv_dir $out_dir")
end
