# =============================================================================
# Clustering Helpers
# =============================================================================
# Public helpers used by GUIbiont exports and by users who want to reproduce the
# clustering data-preparation path without copying GUI-side glue code.
# =============================================================================

using CSV, DataFrames
using Clustering: silhouettes, clustering_quality
using Distributions: Normal, cdf

function apply_blank_timeseries(
    curves::Matrix{Float64},
    blank_timeseries::Vector{Float64};
    method::Symbol = :pointbypoint,
)::Matrix{Float64}
    n_tp = size(curves, 2)
    length(blank_timeseries) >= n_tp ||
        throw(ArgumentError("blank_timeseries length $(length(blank_timeseries)) < n_timepoints $n_tp"))
    ts = blank_timeseries[1:n_tp]
    out = copy(curves)

    if method == :pointbypoint
        for i in axes(out, 1)
            out[i, :] = out[i, :] .- ts
        end
    elseif method == :shift
        finite = filter(isfinite, ts)
        isempty(finite) && return out
        blank_mean = mean(finite)
        for i in axes(out, 1)
            corrected = out[i, :] .- blank_mean
            finite_corrected = filter(isfinite, corrected)
            shift = isempty(finite_corrected) ? 0.0 : min(0.0, minimum(finite_corrected))
            out[i, :] = corrected .- shift
        end
    elseif method == :clip
        finite = filter(isfinite, ts)
        isempty(finite) && return out
        blank_mean = mean(finite)
        for i in axes(out, 1)
            out[i, :] = max.(out[i, :] .- blank_mean, 0.0)
        end
    else
        throw(ArgumentError("Unknown blank correction method: $method"))
    end
    return out
end

function detect_blank_indices(
    curves::Matrix{Float64},
    times::Vector{Float64};
    flat_p_thr::Float64 = 0.05,
    flat_range_thr::Float64 = 0.005,
    od_percentile::Float64 = 0.10,
)::Vector{Int}
    length(times) == size(curves, 2) ||
        throw(ArgumentError("times length $(length(times)) != n_timepoints $(size(curves, 2))"))
    length(times) < 3 && return Int[]

    t_c = times .- mean(times)
    mean_ods = Float64[]
    flat_flags = Bool[]
    for i in axes(curves, 1)
        y = curves[i, :]
        finite_mask = isfinite.(y)
        nf = sum(finite_mask)
        if nf < 3
            push!(mean_ods, NaN); push!(flat_flags, false); continue
        end
        yf = y[finite_mask]
        tcf = t_c[finite_mask]
        push!(mean_ods, mean(yf))

        if maximum(yf) - minimum(yf) < flat_range_thr
            push!(flat_flags, true); continue
        end
        ss_t = sum(tcf .^ 2)
        if ss_t < 1e-12
            push!(flat_flags, true); continue
        end
        slope = sum(tcf .* yf) / ss_t
        yhat = mean(yf) .+ slope .* tcf
        s2 = sum((yf .- yhat) .^ 2) / (nf - 2)
        se = sqrt(max(s2, 0.0) / ss_t)
        t_stat = se < 1e-12 ? Inf : abs(slope / se)
        p_value = 2 * (1 - cdf(Normal(), t_stat))
        push!(flat_flags, p_value >= flat_p_thr)
    end

    finite_means = filter(isfinite, mean_ods)
    isempty(finite_means) && return Int[]
    od_thr = quantile(finite_means, od_percentile)
    return [i for i in axes(curves, 1)
            if flat_flags[i] && isfinite(mean_ods[i]) && mean_ods[i] <= od_thr]
end

function fill_nonfinite_colmean(curves::Matrix{Float64})::Matrix{Float64}
    out = copy(curves)
    for j in axes(out, 2)
        finite = filter(isfinite, out[:, j])
        fill_value = isempty(finite) ? 0.0 : mean(finite)
        for i in axes(out, 1)
            isfinite(out[i, j]) || (out[i, j] = fill_value)
        end
    end
    return out
end

function common_time_grid(
    times_all::Vector{Vector{Float64}};
    n_grid::Int = 100,
    q_low::Float64 = 0.05,
    q_high::Float64 = 0.95,
)::Vector{Float64}
    n = max(2, n_grid)
    starts = [minimum(t) for t in times_all if !isempty(t)]
    ends = [maximum(t) for t in times_all if !isempty(t)]
    isempty(starts) && return Float64[]
    t0 = quantile(starts, q_low)
    t1 = quantile(ends, q_high)
    t0 >= t1 && ((t0, t1) = (minimum(starts), maximum(ends)))
    return collect(range(t0, t1; length=n))
end

function interpolate_curves_to_grid(
    curves_all::Vector{Vector{Float64}},
    times_all::Vector{Vector{Float64}},
    grid::Vector{Float64},
)::Matrix{Float64}
    length(curves_all) == length(times_all) ||
        throw(ArgumentError("curves_all and times_all must have the same length"))
    out = Matrix{Float64}(undef, length(curves_all), length(grid))
    for (i, (y, t)) in enumerate(zip(curves_all, times_all))
        length(y) == length(t) ||
            throw(ArgumentError("curve $i: values and times must have the same length"))
        if isempty(y)
            out[i, :] .= NaN
            continue
        end
        ord = sortperm(t)
        t_s = t[ord]
        y_s = y[ord]
        for (j, tg) in enumerate(grid)
            if tg <= t_s[1]
                out[i, j] = y_s[1]
            elseif tg >= t_s[end]
                out[i, j] = y_s[end]
            else
                k = clamp(searchsortedlast(t_s, tg), 1, length(t_s) - 1)
                frac = (tg - t_s[k]) / (t_s[k + 1] - t_s[k])
                out[i, j] = y_s[k] + frac * (y_s[k + 1] - y_s[k])
            end
        end
    end
    return out
end

function _centroids_for_quality(X::Matrix{Float64}, ids::Vector{Int})::Matrix{Float64}
    centers = Matrix{Float64}(undef, maximum(ids), size(X, 2))
    for cid in 1:maximum(ids)
        centers[cid, :] = vec(mean(X[ids .== cid, :], dims=1))
    end
    return centers
end

function cluster_quality_indices(
    curves::Matrix{Float64},
    ids::Vector{Int};
    zscore_curves::Bool = true,
    max_pairwise_n::Int = 5000,
)::Dict{String,Any}
    length(ids) == size(curves, 1) ||
        throw(ArgumentError("ids length $(length(ids)) != n_curves $(size(curves, 1))"))
    keep = findall(!=(0), ids)
    X = curves[keep, :]
    raw_labels = ids[keep]
    unique_labels = sort(unique(raw_labels))
    label_map = Dict(label => i for (i, label) in enumerate(unique_labels))
    labels = [label_map[label] for label in raw_labels]
    q = Dict{String,Any}()
    if length(labels) < 3 || length(unique(labels)) < 2
        for key in ("silhouette_mean", "silhouettes", "dunn", "davies_bouldin",
                    "calinski_harabasz", "xie_beni")
            q[key] = nothing
        end
        return q
    end

    Xq = zscore_curves ? _zscore_rows(X) : X
    if length(labels) <= max_pairwise_n
        dmat = _pairwise_euclidean(Xq)
        sil = try silhouettes(labels, dmat) catch; nothing end
        q["silhouettes"] = sil === nothing ? nothing : collect(Float64, sil)
        q["silhouette_mean"] = sil === nothing ? nothing : mean(sil)
        q["dunn"] = try Float64(clustering_quality(Xq', labels; quality_index=:dunn)) catch; nothing end
    else
        q["silhouettes"] = nothing
        q["silhouette_mean"] = nothing
        q["dunn"] = nothing
    end

    centers = _centroids_for_quality(Xq, labels)'
    q["davies_bouldin"] = try Float64(clustering_quality(Xq', centers, labels; quality_index=:davies_bouldin)) catch; nothing end
    q["calinski_harabasz"] = try Float64(clustering_quality(Xq', centers, labels; quality_index=:calinski_harabasz)) catch; nothing end
    q["xie_beni"] = try Float64(clustering_quality(Xq', centers, labels; quality_index=:xie_beni)) catch; nothing end
    return q
end

_to_f64_vector(col) = [try parse(Float64, string(v)) catch; NaN end for v in col]

function _load_clustering_csv(path::String)
    raw_df = CSV.read(path, DataFrame)
    col_names = names(raw_df)
    length(col_names) >= 2 || error("CSV must have at least time + one curve column")
    time_col_idx = 1
    first_name = col_names[1]
    first_type = nonmissingtype(eltype(raw_df[!, first_name]))
    if (first_name == "Column1" || startswith(first_name, "Unnamed:")) && first_type <: Number
        time_col_idx = 2
    end
    times_raw = _to_f64_vector(raw_df[!, col_names[time_col_idx]])
    any(isnan, times_raw) && (times_raw = Float64.(0:(nrow(raw_df) - 1)))
    labels = String.(col_names[(time_col_idx + 1):end])
    curves_all = [_to_f64_vector(raw_df[!, label]) for label in labels]
    return [times_raw for _ in labels], curves_all, labels, Vector{Vector{Float64}}(), String[]
end

function _annotation_blank_wells(path::String)
    blanks = Set{String}()
    isfile(path) || return blanks
    ann = CSV.read(path, DataFrame; header=false, silencewarnings=true, stringtype=String)
    for i in 1:nrow(ann)
        ncol(ann) >= 2 || continue
        marker = string(ann[i, 2])
        marker in ("b", "X", "x") && push!(blanks, string(ann[i, 1]))
    end
    return blanks
end

function _load_clustering_experiments(clean_data_path::String, experiments::Vector{String})
    times_all = Vector{Vector{Float64}}()
    curves_all = Vector{Vector{Float64}}()
    labels = String[]
    blank_curves_all = Vector{Vector{Float64}}()
    blank_labels = String[]
    for exp_name in experiments
        data_file = joinpath(clean_data_path, exp_name, "data_channel_1.csv")
        annotation_file = joinpath(clean_data_path, exp_name, "annotation_clean.csv")
        isfile(data_file) || continue
        gd_raw = CSV.read(data_file, DataFrame; header=1, silencewarnings=true)
        blank_wells = _annotation_blank_wells(annotation_file)
        time_col = names(gd_raw)[1]
        time_data = gd_raw[!, time_col]
        time_numeric = if eltype(time_data) <: AbstractString
            try [0.0; [parse(Float64, string(t)) for t in time_data[2:end]]]
            catch; Float64.(0:(nrow(gd_raw) - 1)) end
        else
            Float64.(time_data)
        end
        for well in names(gd_raw)[2:end]
            od = _to_f64_vector(gd_raw[!, well])
            if well in blank_wells
                push!(blank_curves_all, od)
                push!(blank_labels, "$(exp_name)/$(well)")
            else
                push!(times_all, time_numeric)
                push!(curves_all, od)
                push!(labels, "$(exp_name)/$(well)")
            end
        end
    end
    return times_all, curves_all, labels, blank_curves_all, blank_labels
end

function _strip_nan_tail_for_clustering(t::Vector{Float64}, y::Vector{Float64})
    last_valid = findlast(!isnan, y)
    last_valid === nothing && return t[1:min(2, end)], fill(0.0, min(2, length(y)))
    last_valid = max(last_valid, 2)
    return t[1:last_valid], y[1:last_valid]
end

function prepare_clustering_data(;
    csv_path::String = "",
    clean_data_path::String = "",
    experiments::Vector{String} = String[],
    interpolate::Bool = false,
    interp_n::Int = 100,
    interp_quantile_lo::Float64 = 0.05,
    interp_quantile_hi::Float64 = 0.95,
    auto_detect_blanks::Bool = true,
    subtract_blank::Bool = false,
    blank_method::Symbol = :pointbypoint,
    blank_range_thr::Float64 = 0.005,
    blank_od_percentile::Float64 = 0.10,
)::GrowthData
    times_all, curves_all, labels, blank_curves_all, _ = if !isempty(csv_path)
        _load_clustering_csv(csv_path)
    else
        isempty(clean_data_path) && error("clean_data_path is required when csv_path is empty")
        isempty(experiments) && error("experiments must not be empty when csv_path is empty")
        _load_clustering_experiments(clean_data_path, experiments)
    end
    isempty(curves_all) && error("No data loaded")

    if interpolate
        times = common_time_grid(times_all; n_grid=max(10, interp_n),
                                 q_low=interp_quantile_lo, q_high=interp_quantile_hi)
        isempty(times) && error("Could not build interpolation grid")
        curves = interpolate_curves_to_grid(curves_all, times_all, times)
    elseif any(c -> any(isnan, c), curves_all)
        clean_times = Vector{Vector{Float64}}()
        clean_curves = Vector{Vector{Float64}}()
        clean_labels = String[]
        for i in eachindex(curves_all)
            ct, cy = _strip_nan_tail_for_clustering(times_all[i], curves_all[i])
            length(ct) < 2 && continue
            push!(clean_times, ct); push!(clean_curves, cy); push!(clean_labels, labels[i])
        end
        isempty(clean_curves) && error("No valid curves after NaN removal")
        igd = IrregularGrowthData(clean_curves, clean_times, clean_labels; step=0.01)
        times, curves, labels = igd.times, igd.curves, clean_labels
        blank_curves_all = Vector{Vector{Float64}}()
    else
        min_len = minimum(length.(curves_all))
        times = times_all[1][1:min_len]
        curves = Matrix{Float64}(undef, length(curves_all), min_len)
        for i in eachindex(curves_all)
            curves[i, :] = curves_all[i][1:min_len]
        end
    end

    if isempty(blank_curves_all) && auto_detect_blanks
        blank_idxs = detect_blank_indices(curves, times;
            flat_range_thr=blank_range_thr, od_percentile=blank_od_percentile)
        if !isempty(blank_idxs)
            blank_curves_all = [Vector(curves[i, :]) for i in blank_idxs]
            keep = setdiff(1:size(curves, 1), blank_idxs)
            curves, labels = curves[keep, :], labels[keep]
        end
    end

    if subtract_blank && !isempty(blank_curves_all)
        ncols = size(curves, 2)
        blen = min(ncols, minimum(length.(blank_curves_all)))
        blank_mat = Matrix{Float64}(undef, length(blank_curves_all), blen)
        for (i, bc) in enumerate(blank_curves_all)
            blank_mat[i, :] = bc[1:blen]
        end
        blank_ts = [mean(filter(isfinite, blank_mat[:, t])) for t in 1:blen]
        blank_ts_full = length(blank_ts) >= ncols ? blank_ts[1:ncols] :
                        vcat(blank_ts, fill(blank_ts[end], ncols - length(blank_ts)))
        curves = apply_blank_timeseries(curves, blank_ts_full; method=blank_method)
    end

    return GrowthData(fill_nonfinite_colmean(curves), times, labels)
end

export detect_blank_indices
export apply_blank_timeseries
export fill_nonfinite_colmean
export common_time_grid
export interpolate_curves_to_grid
export cluster_quality_indices
export prepare_clustering_data
