# examples/04_irregular_time_clustering.jl
#
# Validates IrregularGrowthData by running the standalone approach
# (ported from newClusteringIdea.jl) and the KinBiont approach on the
# same synthetic data and comparing results.
#
# Run from KinBiont.jl root:
#   julia --project examples/04_irregular_time_clustering.jl

using Random
using Statistics
using Clustering
using Kinbiont

# ─────────────────────────────────────────────────────────────────────────────
# Shared parameters
# ─────────────────────────────────────────────────────────────────────────────
const N_CURVES    = 240
const K           = 6
const STEP        = 0.01
const N_INIT      = 15
const KMEANS_SEED = 2026
const DATA_SEED   = 2026

# ─────────────────────────────────────────────────────────────────────────────
# Data generation helpers (from newClusteringIdea.jl)
# ─────────────────────────────────────────────────────────────────────────────
function logistic_curve(tau, midpoint, slope, amp, baseline)
    baseline .+ amp ./ (1 .+ exp.(-slope .* (tau .- midpoint)))
end

function exponential_curve(tau, rate, amp, baseline)
    baseline .+ amp .* ((exp.(rate .* tau) .- 1.0) ./ (exp(rate) - 1.0))
end

function make_shape_curve(shape_name::String, tau::Vector{Float64})
    baseline = 0.03 + 0.02 * rand()
    amp      = 0.6  + 0.5  * rand()
    y = if shape_name == "early"
        logistic_curve(tau, 0.18 + 0.04 * randn(), 14.0 + 2.0 * randn(), amp, baseline)
    elseif shape_name == "mild-early"
        logistic_curve(tau, 0.28 + 0.05 * randn(),  9.0 + 1.5 * randn(), amp, baseline)
    elseif shape_name == "mild"
        logistic_curve(tau, 0.48 + 0.05 * randn(),  8.0 + 1.0 * randn(), amp, baseline)
    elseif shape_name == "mild-late"
        logistic_curve(tau, 0.62 + 0.05 * randn(),  7.0 + 1.0 * randn(), amp, baseline)
    elseif shape_name == "late-growth"
        logistic_curve(tau, 0.78 + 0.04 * randn(), 12.0 + 2.0 * randn(), amp, baseline)
    elseif shape_name == "exponential"
        exponential_curve(tau, 4.5 + 0.7 * randn(), amp, baseline)
    else
        error("Unknown shape: $shape_name")
    end
    noise_sd = 0.01 + 0.005 * rand()
    return max.(y .+ noise_sd .* randn(length(tau)), 0.0)
end

function generate_irregular_times(n_points::Int)
    start_time = rand() < 0.5 ? 0.0 : rand(5.0:1.0:20.0)
    steps      = rand([15.0, 30.0, 45.0, 60.0], n_points - 1)
    times      = Vector{Float64}(undef, n_points)
    times[1]   = start_time
    for i in 2:n_points
        times[i] = times[i-1] + steps[i-1]
    end
    return times
end

function normalize_times_01(times::Vector{Float64})
    tmin = minimum(times); tmax = maximum(times)
    span = tmax - tmin
    span <= 0.0 && return zeros(Float64, length(times))
    return (times .- tmin) ./ span
end

function generate_dataset(n_curves::Int; min_pts::Int=50, max_pts::Int=200, seed::Int=2026)
    Random.seed!(seed)
    shapes       = ["early","mild-early","mild","mild-late","late-growth","exponential"]
    times_list   = Vector{Vector{Float64}}(undef, n_curves)
    times01_list = Vector{Vector{Float64}}(undef, n_curves)
    values_list  = Vector{Vector{Float64}}(undef, n_curves)
    true_labels  = Vector{String}(undef, n_curves)
    for i in 1:n_curves
        shape      = shapes[mod1(i, length(shapes))]
        n_pts      = rand(min_pts:max_pts)
        raw_t      = generate_irregular_times(n_pts)
        tau        = normalize_times_01(raw_t)
        vals       = make_shape_curve(shape, tau)
        times_list[i]   = raw_t
        times01_list[i] = tau
        values_list[i]  = vals
        true_labels[i]  = shape
    end
    return times_list, times01_list, values_list, true_labels
end

# ─────────────────────────────────────────────────────────────────────────────
# Standalone helpers (from newClusteringIdea.jl)
# ─────────────────────────────────────────────────────────────────────────────
function standalone_build_union_grid(t01_list; step=0.01)
    points = Set{Float64}()
    for t in t01_list
        for v in t
            push!(points, clamp(round(v / step) * step, 0.0, 1.0))
        end
    end
    grid = sort!(collect(points))
    grid[1] > 0.0  && pushfirst!(grid, 0.0)
    grid[end] < 1.0 && push!(grid, 1.0)
    return grid
end

function standalone_interp(x, y, x_new)
    n   = length(x)
    out = Vector{Float64}(undef, length(x_new))
    j   = 1
    for (i, xi) in enumerate(x_new)
        if xi <= x[1];      out[i] = y[1];   continue; end
        if xi >= x[end];    out[i] = y[end]; continue; end
        while j < n-1 && x[j+1] < xi; j += 1; end
        w      = (xi - x[j]) / (x[j+1] - x[j])
        out[i] = (1.0 - w) * y[j] + w * y[j+1]
    end
    return out
end

function standalone_resample(t01_list, values_list, grid; step=0.01)
    n = length(t01_list)
    X = Matrix{Float64}(undef, n, length(grid))
    for i in 1:n
        t_s  = [clamp(round(v / step) * step, 0.0, 1.0) for v in t01_list[i]]
        y    = values_list[i]
        seen = Dict{Float64,Bool}()
        keep = [begin !haskey(seen, tv) ? (seen[tv]=true; true) : false end for tv in t_s]
        t_u  = t_s[keep];  y_u = y[keep]
        perm = sortperm(t_u);  t_u = t_u[perm];  y_u = y_u[perm]
        X[i, :] = standalone_interp(t_u, y_u, grid)
    end
    return X
end

function standalone_zscore_rows(X)
    Z = similar(X)
    for i in axes(X, 1)
        row = X[i, :];  mu = mean(row);  s = std(row)
        Z[i, :] = s < 1e-12 ? zeros(size(X, 2)) : (row .- mu) ./ s
    end
    return Z
end

function standalone_kmeans_best(X, k; n_init=15, maxiter=1000, seed=2026)
    Xt   = permutedims(X)
    best = nothing
    best_cost = Inf
    for rep in 1:n_init
        rng = MersenneTwister(seed + rep)
        r   = kmeans(Xt, k; maxiter=maxiter, rng=rng)
        if r.totalcost < best_cost
            best      = r
            best_cost = r.totalcost
        end
    end
    return best
end

# ─────────────────────────────────────────────────────────────────────────────
# Generate data (shared between both approaches)
# ─────────────────────────────────────────────────────────────────────────────
println("Generating $N_CURVES synthetic curves (seed=$DATA_SEED)...")
times_list, times01_list, values_list, true_labels =
    generate_dataset(N_CURVES; min_pts=50, max_pts=200, seed=DATA_SEED)

# ─────────────────────────────────────────────────────────────────────────────
# Section A — Standalone approach
# ─────────────────────────────────────────────────────────────────────────────
println("\n=== Section A: Standalone approach ===")
union_grid_sa  = standalone_build_union_grid(times01_list; step=STEP)
X_sa           = standalone_resample(times01_list, values_list, union_grid_sa; step=STEP)
Z_sa           = standalone_zscore_rows(X_sa)
km_sa          = standalone_kmeans_best(Z_sa, K; n_init=N_INIT, seed=KMEANS_SEED)
assign_sa      = assignments(km_sa)
println("Union grid length: $(length(union_grid_sa))")
println("K-means WCSS: $(km_sa.totalcost)")
println("Cluster sizes: ", [sum(assign_sa .== c) for c in 1:K])

# ─────────────────────────────────────────────────────────────────────────────
# Section B — KinBiont approach
# ─────────────────────────────────────────────────────────────────────────────
println("\n=== Section B: KinBiont IrregularGrowthData approach ===")
kb_data = IrregularGrowthData(values_list, times_list, true_labels; step=STEP)
kb_opts = FitOptions(
    cluster              = true,
    n_clusters           = K,
    cluster_trend_test   = false,
    kmeans_seed          = KMEANS_SEED,
    kmeans_n_init        = N_INIT,
    kmeans_max_iters     = 1000,
)
kb_result  = preprocess(kb_data, kb_opts)
assign_kb  = kb_result.clusters
println("Union grid length: $(length(kb_data.times))")
println("K-means WCSS: $(kb_result.wcss)")
println("Cluster sizes: ", [sum(assign_kb .== c) for c in 1:K])

# ─────────────────────────────────────────────────────────────────────────────
# Section C — Resampled matrix comparison
# ─────────────────────────────────────────────────────────────────────────────
println("\n=== Section C: Resampled matrix comparison ===")
grids_match = union_grid_sa == kb_data.times
println("Union grids identical: $grids_match")
if grids_match
    max_diff = maximum(abs.(X_sa .- kb_data.curves))
    println("Max element-wise difference in resampled matrices: $max_diff")
    println("Matrices identical (tol=1e-12): $(max_diff < 1e-12)")
end

# ─────────────────────────────────────────────────────────────────────────────
# Section D — Cluster assignment comparison
# ─────────────────────────────────────────────────────────────────────────────
println("\n=== Section D: Cluster assignment comparison ===")

# Contingency table
table = zeros(Int, K, K)
for i in eachindex(assign_sa)
    table[assign_sa[i], assign_kb[i]] += 1
end
println("Contingency table (standalone rows × KinBiont cols):")
for r in 1:K
    println("  SA cluster $r: ", table[r, :])
end

# Greedy label matching (works well for well-separated clusters)
function greedy_label_match(tbl::Matrix{Int})
    k    = size(tbl, 1)
    perm = zeros(Int, k)
    used = Set{Int}()
    row_order = sortperm([maximum(tbl[r, :]) for r in 1:k]; rev=true)
    for r in row_order
        avail   = [j for j in 1:k if j ∉ used]
        best_j  = avail[argmax(tbl[r, avail])]
        perm[r] = best_j
        push!(used, best_j)
    end
    return perm
end

perm      = greedy_label_match(table)
n_agree   = sum(table[r, perm[r]] for r in 1:K)
agreement = n_agree / N_CURVES
println("\nBest label permutation (SA→KB): $perm")
println("Curves in agreement: $n_agree / $N_CURVES  ($(round(100*agreement; digits=1))%)")

# Per-true-shape breakdown
println("\nPer-shape cluster agreement:")
shapes = sort(unique(true_labels))
for shape in shapes
    idx_shape   = findall(==(shape), true_labels)
    sa_dominant = argmax([count(==(c), assign_sa[idx_shape]) for c in 1:K])
    kb_dominant = perm[sa_dominant]
    kb_agree    = count(==(kb_dominant), assign_kb[idx_shape])
    println("  $shape: SA dominant=cluster $sa_dominant, KB mapped=$kb_dominant, " *
            "agreement=$(kb_agree)/$(length(idx_shape))")
end
