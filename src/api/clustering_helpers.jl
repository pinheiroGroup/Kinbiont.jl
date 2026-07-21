# =============================================================================
# Clustering Helpers
# =============================================================================
# Public helpers used by GUIbiont exports and by users who want to reproduce the
# clustering data-preparation path without copying GUI-side glue code.
# =============================================================================

using CSV, DataFrames
using Clustering: silhouettes, clustering_quality

"""
    apply_blank_timeseries(curves, blank_timeseries; method=:pointbypoint) -> Matrix{Float64}

Subtract a per-timepoint blank trace from every curve in `curves` (rows = curves,
columns = timepoints) and return a new matrix.

`method` selects the correction style:
- `:pointbypoint` — subtract `blank_timeseries[t]` from each timepoint. Non-finite
  blank values are treated as zero (no correction at that timepoint).
- `:shift` — subtract the mean of the (finite) blank trace, then shift each curve
  so its minimum is ≥ 0.
- `:clip` — subtract the mean of the (finite) blank trace, then clamp negatives
  to 0.

`blank_timeseries` must be at least as long as the number of timepoints.
"""
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
        # Replace non-finite blank entries with 0 so they leave the curve untouched
        # at that timepoint (consistent with :shift / :clip which filter to finite).
        ts_safe = [isfinite(v) ? v : 0.0 for v in ts]
        for i in axes(out, 1)
            out[i, :] = out[i, :] .- ts_safe
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

"""
    detect_blank_indices(curves, times; flat_p_thr=0.05, flat_range_thr=0.005, od_percentile=0.10) -> Vector{Int}

Return the row indices in `curves` that look like blank wells: curves that are
both *flat* (no significant linear trend at p ≥ `flat_p_thr` **or** OD range
< `flat_range_thr`) and have a mean OD in the lowest `od_percentile` quantile of
all finite curve means.
"""
function detect_blank_indices(
    curves::Matrix{Float64},
    times::Vector{Float64};
    flat_p_thr::Float64 = 0.05,
    flat_range_thr::Float64 = 0.005,
    od_percentile::Float64 = 0.10,
)::Vector{Int}
    flat_flags = _flat_curve_mask(curves, times;
        p_threshold=flat_p_thr,
        range_threshold=flat_range_thr,
    )

    mean_ods = Float64[]
    for i in axes(curves, 1)
        finite = filter(isfinite, curves[i, :])
        push!(mean_ods, isempty(finite) ? NaN : mean(finite))
    end

    finite_means = filter(isfinite, mean_ods)
    isempty(finite_means) && return Int[]
    od_thr = quantile(finite_means, od_percentile)
    return [i for i in axes(curves, 1)
            if flat_flags[i] && isfinite(mean_ods[i]) && mean_ods[i] <= od_thr]
end

"""
    fill_nonfinite_colmean(curves) -> Matrix{Float64}

Return a copy of `curves` (rows = curves, columns = timepoints) where each
non-finite entry is replaced by the mean of the finite values in its column.
Columns with no finite values are filled with 0.
"""
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

"""
    common_time_grid(times_all; n_grid=100, q_low=0.05, q_high=0.95) -> Vector{Float64}

Build a uniform interpolation grid covering the bulk of the time ranges in
`times_all`. The grid starts at the `q_low` quantile of all curve start times
and ends at the `q_high` quantile of all curve end times. Falls back to the
full min/max range if the quantile-based interval is degenerate.
"""
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

"""
    interpolate_curves_to_grid(curves_all, times_all, grid) -> Matrix{Float64}

Linearly interpolate each curve in `curves_all` (paired with its own time vector
in `times_all`) onto the common `grid`. Values outside a curve's own time range
are clamped to the nearest endpoint. Empty curves produce a row of `NaN`.
"""
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

"""
    cluster_quality_indices(curves, ids; zscore_curves=true, max_pairwise_n=5000) -> Dict{String,Any}

Compute clustering quality indices for `curves` (rows = curves) given cluster
assignments `ids`. Curves with `ids[i] == 0` (DBSCAN noise) are excluded.

Returns a dictionary with keys `silhouette_mean`, `silhouettes`, `dunn`,
`davies_bouldin`, `calinski_harabasz`, `xie_beni`. Entries that cannot be
computed (too few labels, too many samples for the pairwise step, or a
Clustering.jl failure) are returned as `nothing`; the underlying exception is
emitted on the `@debug` log level.
"""
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
        sil = try
            silhouettes(labels, dmat)
        catch e
            @debug "silhouettes computation failed" exception=e
            nothing
        end
        q["silhouettes"] = sil === nothing ? nothing : collect(Float64, sil)
        q["silhouette_mean"] = sil === nothing ? nothing : mean(sil)
        q["dunn"] = try
            Float64(clustering_quality(Xq', labels; quality_index=:dunn))
        catch e
            @debug "dunn computation failed" exception=e
            nothing
        end
    else
        q["silhouettes"] = nothing
        q["silhouette_mean"] = nothing
        q["dunn"] = nothing
    end

    centers = _compute_centroids(Xq, labels, maximum(labels))'
    q["davies_bouldin"] = try
        Float64(clustering_quality(Xq', centers, labels; quality_index=:davies_bouldin))
    catch e
        @debug "davies_bouldin computation failed" exception=e
        nothing
    end
    q["calinski_harabasz"] = try
        Float64(clustering_quality(Xq', centers, labels; quality_index=:calinski_harabasz))
    catch e
        @debug "calinski_harabasz computation failed" exception=e
        nothing
    end
    q["xie_beni"] = try
        Float64(clustering_quality(Xq', centers, labels; quality_index=:xie_beni))
    catch e
        @debug "xie_beni computation failed" exception=e
        nothing
    end
    return q
end

# Parse a column of mixed numeric/string values into a Float64 vector.
# Non-parseable entries become `NaN`.
function _to_f64_vector(col)
    out = Vector{Float64}(undef, length(col))
    for (i, v) in enumerate(col)
        if v isa Real
            out[i] = Float64(v)
        else
            parsed = tryparse(Float64, string(v))
            out[i] = parsed === nothing ? NaN : parsed
        end
    end
    return out
end

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
            parsed = [tryparse(Float64, string(t)) for t in time_data[2:end]]
            if any(==(nothing), parsed)
                Float64.(0:(nrow(gd_raw) - 1))
            else
                [0.0; Float64.(parsed)]
            end
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

# Strip a trailing NaN tail from a curve so it can be resampled by
# IrregularGrowthData. Returns `nothing` if the curve has fewer than 2 finite
# leading values (callers must drop such curves rather than fabricate data).
function _strip_nan_tail_for_clustering(t::Vector{Float64}, y::Vector{Float64})
    last_valid = findlast(!isnan, y)
    last_valid === nothing && return nothing
    last_valid < 2 && return nothing
    return t[1:last_valid], y[1:last_valid]
end

"""
    prepare_clustering_data(; kwargs...) -> GrowthData

Load growth curve data and assemble a `GrowthData` ready for the clustering
pipeline, mirroring GUIbiont's data-preparation steps.

Source (one of):
- `csv_path` — single CSV file: first column is time (or `"Column1"`/`"Unnamed:"`
  sentinel skipped), remaining columns are wells.
- `clean_data_path` + `experiments` — GUIbiont layout where each experiment is a
  subdirectory containing `data_channel_1.csv` and `annotation_clean.csv`
  (annotation markers `"b"`, `"X"`, `"x"` flag blanks).

Pipeline options:
- `interpolate` — resample all curves onto a common quantile-trimmed grid.
- `auto_detect_blanks` — flag flat, low-OD curves as blanks via
  [`detect_blank_indices`](@ref) when no explicit blanks were provided.
- `subtract_blank` — subtract the per-timepoint mean of the blank curves using
  [`apply_blank_timeseries`](@ref) with `blank_method`.

Curves that become invalid (e.g. all-NaN after trimming) are dropped silently
**from the curve list, not fabricated as zeros**.
"""
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
    times_all, curves_all, labels, blank_curves_all, blank_labels = if !isempty(csv_path)
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
            stripped = _strip_nan_tail_for_clustering(times_all[i], curves_all[i])
            stripped === nothing && continue
            ct, cy = stripped
            push!(clean_times, ct); push!(clean_curves, cy); push!(clean_labels, labels[i])
        end
        isempty(clean_curves) && error("No valid curves after NaN removal")
        igd = IrregularGrowthData(clean_curves, clean_times, clean_labels; step=0.01)
        times, curves, labels = igd.times, igd.curves, clean_labels
        blank_curves_all = Vector{Vector{Float64}}()
        blank_labels     = String[]
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
            blank_labels     = labels[blank_idxs]
            keep = setdiff(1:size(curves, 1), blank_idxs)
            curves, labels = curves[keep, :], labels[keep]
        end
    end

    # Blank correction is performed separately within each loaded experiment,
    # matching the GUIbiont clustering route: each experiment's wells are
    # corrected using only that experiment's own blanks.
    if subtract_blank && !isempty(blank_curves_all)
        curves = _subtract_blanks_per_experiment(curves, labels,
            blank_curves_all, blank_labels, blank_method, !isempty(csv_path))
    end

    return GrowthData(fill_nonfinite_colmean(curves), times, labels)
end

# Experiment key for a series label (mirrors GUIbiont's clustering route). In CSV
# mode all series share one group; in experiment mode the key is the substring
# before the first '/'.
function _clustering_experiment_key(label::AbstractString, csv_mode::Bool)::String
    csv_mode && return "__all__"
    idx = findfirst('/', label)
    return idx === nothing ? String(label) : String(label[1:prevind(label, idx)])
end

# Per-experiment pooled blank subtraction: each experiment's wells are corrected
# with a per-timepoint mean blank timeseries built only from that experiment's
# blanks. Series whose experiment has no blanks are returned unchanged. Returns a
# new curves matrix (never mutates the input).
function _subtract_blanks_per_experiment(
    curves::Matrix{Float64},
    labels::Vector{String},
    blank_curves::Vector{Vector{Float64}},
    blank_labels::Vector{String},
    blank_method::Symbol,
    csv_mode::Bool,
)::Matrix{Float64}
    out   = copy(curves)
    ncols = size(curves, 2)
    isempty(blank_curves) && return out
    blank_exps = [_clustering_experiment_key(l, csv_mode) for l in blank_labels]
    curve_exps = [_clustering_experiment_key(l, csv_mode) for l in labels]
    for exp in unique(blank_exps)
        bidx = findall(==(exp), blank_exps)
        isempty(bidx) && continue
        rows = findall(==(exp), curve_exps)
        isempty(rows) && continue
        blen = min(ncols, minimum(length.(blank_curves[bidx])))
        blen < 1 && continue
        bmat = Matrix{Float64}(undef, length(bidx), blen)
        for (i, bi) in enumerate(bidx)
            bmat[i, :] = blank_curves[bi][1:blen]
        end
        blank_ts      = [mean(filter(isfinite, bmat[:, t])) for t in 1:blen]
        blank_ts_full = length(blank_ts) >= ncols ? blank_ts[1:ncols] :
                        vcat(blank_ts, fill(blank_ts[end], ncols - length(blank_ts)))
        out[rows, :] = apply_blank_timeseries(out[rows, :], blank_ts_full; method=blank_method)
    end
    return out
end

export detect_blank_indices
export apply_blank_timeseries
export fill_nonfinite_colmean
export common_time_grid
export interpolate_curves_to_grid
export cluster_quality_indices
export prepare_clustering_data
