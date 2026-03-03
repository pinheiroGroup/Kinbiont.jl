module timeSeriesKmeansLib

using LinearAlgebra, Statistics, Random, CSV, DataFrames, Base.Threads

export  zscore_vector, euclid_sqdist,
        build_exp_protos, km_random_init, assign_step!, update_step!, kmeans_parallel,
        zscore_rows, best_const_distance_with_level,
        best_exp_distance, best_km_distance, best_nonconst,
        curves_clustering

include("mannKendallTest.jl")

# ============================================================
# Normalization & Euclidean Distance Functions
# ============================================================

"Z-score normalize a real vector. If std ~ 0, return zeros."
@inline function zscore_vector(x::AbstractVector{<:Real})
    x64 = Float64.(x)
    mu = mean(x64)
    sigma = std(x64; corrected=false)  # avoids NaN for length(x)==1
    return (sigma < 1e-12) ? zeros(length(x64)) : ((x64 .- mu) ./ sigma)
end

"Squared Euclidean distance."
@inline function euclid_sqdist(u::AbstractVector{<:Real}, v::AbstractVector{<:Real})
    s = 0.0
    @inbounds @simd for j in eachindex(u, v)
        d = float(u[j]) - float(v[j])
        s += d * d
    end
    return s
end

# ==================================================================
# Build Exponential Prototypes for the optional Exponential Cluster
# ==================================================================

"""
Build exponential-like prototypes and normalize them with z-score.
Each prototype is (b^tau - 1)/(b-1) for tau in [0,1], then standardized.
"""
function build_exp_protos(
    time::AbstractVector{<:Real};
    bases::Tuple{Vararg{Int}}=(2^6,2^7,2^8,2^9,2^10,2^11,2^12,
                               2^13,2^14,2^15,2^16)
)
    t64 = Float64.(time)
    tmin, tmax = minimum(t64), maximum(t64)
    delta_t = max(tmax - tmin, 1e-12)
    tau = (t64 .- tmin) ./ delta_t

    protos = Vector{Vector{Float64}}()
    for b in bases
        bf = float(b)
        y = [(bf^tau_j - 1.0) / (bf - 1.0) for tau_j in tau]
        yz = zscore_vector(y)
        push!(protos, yz)
    end
    return protos
end

# ============================================================
# K-means
# ============================================================

"Randomly initialize K centroids from rows of X (Matrix{Float64})."
function km_random_init(X::Matrix{Float64}, k::Int; rng=Random.default_rng())
    n, _ = size(X)
    idxs = rand(rng, 1:n, k)
    return [(@views copy(X[idxs[j], :])) for j in 1:k]
end

"Assign points to nearest centroid; returns total SSE."
function assign_step!(labels::Vector{Int},
                      X::Matrix{Float64},
                      centroids::Vector{Vector{Float64}})
    n = size(X, 1)
    k = length(centroids)
    inertia_local = zeros(Float64, nthreads())

    @threads for i in 1:n
        @views xi = X[i, :]
        best_k = 1
        best_d2 = euclid_sqdist(xi, centroids[1])

        @inbounds for j in 2:k
            d2 = euclid_sqdist(xi, centroids[j])
            if d2 < best_d2
                best_d2 = d2
                best_k = j
            end
        end

        labels[i] = best_k
        inertia_local[threadid()] += best_d2
    end

    return sum(inertia_local)
end

"Update centroid means from labels."
function update_step!(centroids::Vector{Vector{Float64}},
                      X::Matrix{Float64},
                      labels::Vector{Int})
    n, m = size(X)
    k = length(centroids)
    T = nthreads()

    sums_tls = [zeros(m, k) for _ in 1:T]
    cnts_tls = [zeros(Int, k) for _ in 1:T]

    @threads for i in 1:n
        t = threadid()
        lbl = labels[i]  # 1..k
        @views sums_tls[t][:, lbl] .+= X[i, :]
        cnts_tls[t][lbl] += 1
    end

    sums = zeros(size(sums_tls[1]))
    cnts = zeros(Int, k)
    for t in 1:T
        sums .+= sums_tls[t]
        cnts .+= cnts_tls[t]
    end

    for j in 1:k
        if cnts[j] == 0
            @views centroids[j] .= X[rand(1:n), :]
        else
            @views centroids[j] .= view(sums, :, j) ./ cnts[j]
        end
    end

    # Re-normalize centroids (keep shape-based consistency)
    for j in 1:k
        c = centroids[j]
        mu_c = mean(c)
        sigma_c = std(c; corrected=false)
        if sigma_c < 1e-12
            fill!(c, 0.0)
        else
            c .= (c .- mu_c) ./ sigma_c
        end
    end

    return nothing
end

"""
Run parallel K-means with multiple initializations.

X can be any AbstractMatrix{<:Real}; it is converted once to Matrix{Float64}.

Returns:
  - labels (0-based: 0..k-1)
  - centroids (z-normalized)
  - best_sse
  - wcss (== best_sse; returned for API compatibility)
"""
function kmeans_parallel(
    X::AbstractMatrix{<:Real},
    k::Int;
    max_iter::Int=300,
    n_init::Int=3,
    tol::Float64=1e-8,
    rng=Random.default_rng(),
    verbose::Bool=true
)
    Xmat = Array{Float64}(X)
    n, _ = size(Xmat)

    best_sse = Inf
    best_labels_1 = zeros(Int, n)  # internal 1-based labels
    best_centroids = Vector{Vector{Float64}}()

    for rep in 1:n_init
        verbose && println("Initialization $rep / $n_init")

        centroids = km_random_init(Xmat, k; rng=rng)
        labels_1 = zeros(Int, n)   # 1..k
        prev_sse = Inf

        for _ in 1:max_iter
            sse = assign_step!(labels_1, Xmat, centroids)
            update_step!(centroids, Xmat, labels_1)
            if abs(prev_sse - sse) < tol
                break
            end
            prev_sse = sse
        end

        sse_final = assign_step!(labels_1, Xmat, centroids)
        if sse_final < best_sse
            best_sse = sse_final
            best_labels_1 .= labels_1
            best_centroids = deepcopy(centroids)
        end
    end

    verbose && println("Best SSE: $(round(best_sse, digits=6))")

    # Convert labels to 0-based for the public API
    best_labels_0 = copy(best_labels_1)
    best_labels_0 .-= 1  # now 0..k-1

    wcss = best_sse  # WCSS == SSE for k-means with squared Euclidean distance
    return best_labels_0, best_centroids, best_sse, wcss
end

# ============================================================
# Distances and classifiers
# ============================================================

"Row-wise z-score (name kept for API compatibility)."
function zscore_rows(mat::AbstractMatrix{<:Real})
    n, m = size(mat)
    out = Array{Float64}(undef, n, m)
    @inbounds @views for i in 1:n
        out[i, :] = zscore_vector(mat[i, :])
    end
    return out
end

"Compute distance to the best constant level (L = mean(x)), with bias."
function best_const_distance_with_level(
    xnorm::AbstractVector{<:Real};
    use_bias::Bool=true,
    p::Float64=0.70,
    trend_factor::Float64=0.50,
    bias_cap::Float64=0.99
)
    x64 = Float64.(xnorm)
    N = length(x64)
    mu = mean(x64)
    d2 = sum((x64 .- mu) .^ 2) + 1e-18

    if !use_bias
        return d2, mu
    end

    qval = quantile(x64, p)
    Ninband = count(x -> x <= qval, x64)
    frac_inband = Ninband / N

    Lmax = 0
    cur = 0
    @inbounds @simd for t in 1:N
        if x64[t] > qval
            cur += 1
            if cur > Lmax
                Lmax = cur
            end
        else
            cur = 0
        end
    end
    run_max = Lmax / N

    bias = frac_inband * (1 - run_max)

    res = mk_original_test(x64)
    if res.trend == "increasing"
        bias *= trend_factor
    end

    if all(x -> x > qval, x64)
        bias = 0.0
    end

    bias = min(bias, bias_cap)
    comp_bias = 1 - bias
    comp_bias = comp_bias ^ 2  # keep consistent with squared distances

    return d2 * comp_bias, mu
end

"Compute distance to nearest exponential prototype."
@inline function best_exp_distance(
    xnorm::AbstractVector{<:Real},
    exp_protos::Vector{Vector{Float64}}
)
    x64 = Float64.(xnorm)
    bestd = Inf
    @inbounds for prot in exp_protos
        d2 = euclid_sqdist(x64, prot)
        if d2 < bestd
            bestd = d2
        end
    end
    return bestd
end

"Compute distance to nearest K-means centroid. Returns (distance, 0-indexed label)."
function best_km_distance(
    xnorm::AbstractVector{<:Real},
    centroids_norm::Vector{Vector{Float64}}
)
    x64 = Float64.(xnorm)
    bestd = Inf
    best_lbl0 = -1
    @inbounds for c in 1:length(centroids_norm)
        d2 = euclid_sqdist(x64, centroids_norm[c])
        if d2 < bestd
            bestd = d2
            best_lbl0 = c - 1
        end
    end
    return bestd, best_lbl0
end

"Compare constant vs nonconstant distances and assign label."
function best_nonconst(
    xnorm::AbstractVector{<:Real},
    centroids_norm::Vector{Vector{Float64}},
    exp_protos::Vector{Vector{Float64}},
    use_exp_cluster::Bool
)
    best_km, best_lbl0 = isempty(centroids_norm) ? (Inf, -1) :
                         best_km_distance(xnorm, centroids_norm)

    if use_exp_cluster
        best_exp = best_exp_distance(xnorm, exp_protos)

        if isempty(centroids_norm) || best_exp < best_km
            return best_exp, "exponential", -1
        else
            return best_km, "kmeans", best_lbl0
        end
    else
        # exponential cluster disabled: only K-means is used
        return best_km, "kmeans", best_lbl0
    end
end

# ============================================================
# Top-level classifier (no I/O)
# ============================================================

"""
Cluster a full plate matrix X (rows=wells, cols=time).

X and time_vec can be any AbstractMatrix/Vector{<:Real}; internally we
convert once to Float64.

Returns:
  - labels_final (1..K for kmeans, K+1=exponential, K+2=constant)
  - wcss_value   (== SSE from kmeans)
  - centroids_original (length K+2):
        index L is the prototype in original space for label L (L=1..K+2)
  - centroids_norm (length K+2):
        index L is the z-score prototype for label L (L=1..K+2)
"""
function curves_clustering(
    X::AbstractMatrix{<:Real},
    time_vec::AbstractVector{<:Real},
    K::Int;
    q_low::Float64=0.05,
    q_high::Float64=0.95,
    tol_const::Float64=1.5,
    const_qbias::Float64=0.70,
    trend_factor::Float64=0.50,
    use_exp_cluster::Bool=true,
    rng=Random.MersenneTwister(42)
)
    Xmat = Array{Float64}(X)
    tvec = Float64.(time_vec)

    # Build exponential prototypes only if used
    exp_protos = use_exp_cluster ? build_exp_protos(tvec) : Vector{Vector{Float64}}()
    W = size(Xmat, 1)

    # --- quantiles for equal/near-const/non-const
    qlowv  = [quantile(Xmat[i, :], q_low)  for i in 1:W]
    qhighv = [quantile(Xmat[i, :], q_high) for i in 1:W]

    # --- robust medians on low and high tails
    low_med  = Vector{Float64}(undef, W)
    high_med = Vector{Float64}(undef, W)

    for i in 1:W
        row = Xmat[i, :]
        lows  = row[row .<= qlowv[i]]
        highs = row[row .>= qhighv[i]]

        low_med[i]  = isempty(lows)  ? qlowv[i]  : median(lows)
        high_med[i] = isempty(highs) ? qhighv[i] : median(highs)
    end

    # --- equal / near-constant / non-constant
    eps_rel = 1e-6
    equal_mask     = falses(W)
    nearconst_mask = falses(W)

    for i in 1:W
        lm = low_med[i]
        hm = high_med[i]

        # 1) Curva praticamente piatta
        equal_mask[i] = abs(hm - lm) <= eps_rel * (abs(lm) + 1e-6)

        # 2) Near-constant: non ha superato "raddoppio" rispetto al low
        if !equal_mask[i]
            if lm > 0
                nearconst_mask[i] = hm <= tol_const * lm
            else
                nearconst_mask[i] = hm / tol_const <= lm
            end
        end
    end

    nonconst_mask = .!(equal_mask .| nearconst_mask)

    equal_idxs     = findall(equal_mask)
    nearconst_idxs = findall(nearconst_mask)
    nonconst_idxs  = findall(nonconst_mask)

    # --- K-means ONLY on non-const
    centroids_norm_km = Vector{Vector{Float64}}()
    labels_km0 = Int[]       # 0-based labels for nonconst subset
    wcss_value = 0.0

    if !isempty(nonconst_idxs)
        data_km = zscore_rows(Xmat[nonconst_idxs, :])
        labels_km0, centroids_norm_km, _, wcss_value =
            kmeans_parallel(data_km, K; rng=rng, verbose=false)
    else
        # no non-constant curves: no k-means signal, keep wcss_value = 0.0
        wcss_value = 0.0
    end

    # --- final labels
    LABEL_EXP0   = K
    LABEL_CONST0 = K + 1
    labels_final = fill(Int(-1), W)

    # equal -> constant
    for idx in equal_idxs
        labels_final[idx] = LABEL_CONST0
    end

    # near-const: constant (with bias) vs non-const (kmeans/exp)
    for idx in nearconst_idxs
        xi_norm = zscore_vector(@view Xmat[idx, :])
        best_const_adj, _ = best_const_distance_with_level(
            xi_norm; use_bias=true, p=const_qbias, trend_factor=trend_factor
        )
        best_nc, nc_type, nc_lbl0 =
            best_nonconst(xi_norm, centroids_norm_km, exp_protos, use_exp_cluster)
        labels_final[idx] = (best_const_adj <= best_nc) ? LABEL_CONST0 :
                            ((nc_type == "exponential") ? LABEL_EXP0 : nc_lbl0)
    end

    # non-const: kmeans vs exponential
    for idx in nonconst_idxs
        xi_norm = zscore_vector(@view Xmat[idx, :])
        _, nc_type, nc_lbl0 =
            best_nonconst(xi_norm, centroids_norm_km, exp_protos, use_exp_cluster)
        labels_final[idx] = (nc_type == "exponential") ? LABEL_EXP0 : nc_lbl0
    end

    # defensive fallback: any still -1 becomes constant
    for i in eachindex(labels_final)
        if labels_final[i] == -1
            labels_final[i] = LABEL_CONST0
        end
    end

    # --- build prototypes aligned with labels 0..K+1
    # index L+1 is the prototype for label L
    n_time = size(Xmat, 2)
    centroids_original = Vector{Vector{Float64}}(undef, K + 2)
    centroids_norm     = Vector{Vector{Float64}}(undef, K + 2)

    # 1) k-means clusters 0..K-1:
    #    use the actual k-means centroids for normalized prototypes,
    #    and the mean of original curves for raw prototypes.
    for j0 in 0:K-1
        if !isempty(nonconst_idxs) && !isempty(labels_km0)
            mask = labels_km0 .== j0
            if any(mask)
                cluster_rows  = nonconst_idxs[mask]
                cluster_data  = Xmat[cluster_rows, :]
                centroid_orig = vec(mean(cluster_data, dims=1))
                centroids_original[j0+1] = centroid_orig

                # for normalized prototype, we keep the k-means centroid as is
                centroids_norm[j0+1] = centroids_norm_km[j0+1]
            else
                centroids_original[j0+1] = zeros(n_time)
                centroids_norm[j0+1]     = zeros(n_time)
            end
        else
            centroids_original[j0+1] = zeros(n_time)
            centroids_norm[j0+1]     = zeros(n_time)
        end
    end

    # 2) exponential cluster (label K) and constant cluster (label K+1):
    #    we build prototypes directly from the raw curves assigned to them.
    for (label0, idx_out) in ((LABEL_EXP0,   K+1),   # exponential
                            (LABEL_CONST0, K+2))   # constant
        idxs = findall(==(label0), labels_final)
        if !isempty(idxs)
            cluster_data  = Xmat[idxs, :]
            centroid_orig = vec(mean(cluster_data, dims=1))
            centroid_norm = zscore_vector(centroid_orig)
        else
            centroid_orig = zeros(n_time)
            centroid_norm = zeros(n_time)
        end
        centroids_original[idx_out] = centroid_orig
        centroids_norm[idx_out]     = centroid_norm
    end

    # convert labels from 0-based (0..K+1) to 1-based (1..K+2)
    labels_final_1based = similar(labels_final)
    @inbounds for i in eachindex(labels_final)
        labels_final_1based[i] = labels_final[i] + 1
    end

    return labels_final_1based, wcss_value, centroids_original, centroids_norm

end

end # module
