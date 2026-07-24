# =============================================================================
# Clustering Helpers
# =============================================================================
# Public helpers used by GUIbiont exports and by users who want to reproduce the
# clustering data-preparation path without copying GUI-side glue code.
# =============================================================================

using CSV, DataFrames
using Clustering: silhouettes, clustering_quality

"""
    detect_blank_indices(curves, times; flat_p_thr=0.05, flat_range_thr=0.005, od_percentile=0.10) -> Vector{Int}

Return the row indices in `curves` that look like blank wells: curves that are
both *flat* (no significant linear trend under a two-sided Normal approximation
at p ≥ `flat_p_thr` **or** OD range
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
    flat_flags = _normal_approx_flat_mask(curves, times;
        p_threshold=flat_p_thr, range_threshold=flat_range_thr)

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
in `times_all`) onto the common `grid`. Non-finite time/value pairs are ignored;
values outside a curve's finite time range are clamped to the nearest endpoint.
Curves without a finite measurement produce a row of `NaN`.
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
        finite = isfinite.(t) .& isfinite.(y)
        if !any(finite)
            out[i, :] .= NaN
            continue
        end
        t_f, y_f = t[finite], y[finite]
        ord = sortperm(t_f)
        t_s, y_s = t_f[ord], y_f[ord]
        if length(t_s) == 1
            out[i, :] .= y_s[1]
            continue
        end
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

function _normal_approx_flat_mask(
    curves::Matrix{Float64},
    times::Vector{Float64};
    p_threshold::Float64=0.05,
    range_threshold::Union{Nothing,Float64}=nothing,
)::BitVector
    length(times) == size(curves, 2) || throw(ArgumentError(
        "times length $(length(times)) != n_timepoints $(size(curves, 2))"
    ))
    mask = falses(size(curves, 1))
    length(times) < 3 && return mask

    for i in axes(curves, 1)
        finite = isfinite.(curves[i, :]) .& isfinite.(times)
        sum(finite) < 3 && continue
        y, t = curves[i, finite], times[finite]
        if range_threshold !== nothing && maximum(y) - minimum(y) < range_threshold
            mask[i] = true
            continue
        end
        tc = t .- mean(t)
        ss_t = sum(tc .^ 2)
        if ss_t <= 1e-12
            mask[i] = true
            continue
        end
        slope = sum(tc .* y) / ss_t
        residuals = y .- (mean(y) .+ slope .* tc)
        se = sqrt(max(sum(residuals .^ 2) / (length(y) - 2), 0.0) / ss_t)
        if se <= 0 || !isfinite(se)
            mask[i] = abs(slope) < 1e-12
            continue
        end
        p_value = 2 * (1 - cdf(Normal(), abs(slope / se)))
        mask[i] = p_value >= p_threshold
    end
    return mask
end

"""
    apply_grouped_blank_subtraction(curves, times, groups,
                                    blank_curves, blank_times, blank_groups;
                                    method=:pointbypoint)

Correct every sample curve using only blank curves from the same group, with
the same semantics used by single and batch fitting. Subtraction is enabled
only when the mean of all finite blank readings is positive. `:shift` and
`:clip` use that scalar mean; `:pointbypoint` averages blanks by measurement
index. Samples whose group has no usable positive blank are returned unchanged.
Corrected curves are kept at or above `floor` (default `1e-4`).

Returns `(corrected_curves, corrected_mask)`, where `corrected_mask[i]` records
whether sample `i` had a same-group blank available.
"""
function apply_grouped_blank_subtraction(
    curves::Vector{Vector{Float64}},
    times::Vector{Vector{Float64}},
    groups::Vector{String},
    blank_curves::Vector{Vector{Float64}},
    blank_times::Vector{Vector{Float64}},
    blank_groups::Vector{String};
    method::Symbol = :pointbypoint,
    floor::Union{Nothing,Float64} = 1e-4,
)::Tuple{Vector{Vector{Float64}},BitVector}
    length(curves) == length(times) == length(groups) || throw(ArgumentError(
        "curves, times, and groups must have the same length"
    ))
    length(blank_curves) == length(blank_times) == length(blank_groups) ||
        throw(ArgumentError("blank_curves, blank_times, and blank_groups must have the same length"))

    blanks_by_group = Dict{String,Vector{Int}}()
    for (i, group) in enumerate(blank_groups)
        push!(get!(blanks_by_group, group, Int[]), i)
    end

    corrected = deepcopy(curves)
    corrected_mask = falses(length(curves))
    for i in eachindex(curves)
        length(curves[i]) == length(times[i]) || throw(ArgumentError(
            "sample $i: values and times must have the same length"
        ))
        blank_idx = get(blanks_by_group, groups[i], Int[])
        isempty(blank_idx) && continue

        blank_value, blank_trace = _blank_summary(
            blank_curves[blank_idx],
            length(curves[i]),
        )
        blank_value > 0.0 || continue
        subtraction_trace = method == :pointbypoint ?
            blank_trace :
            fill(blank_value, length(curves[i]))
        corrected_row = apply_blank_timeseries(
            reshape(curves[i], 1, :), subtraction_trace; method=method, floor=floor
        )
        corrected[i] = vec(corrected_row)
        corrected_mask[i] = true
    end
    return corrected, corrected_mask
end

"""
    derive_blank_from_non_growing(curves, times, labels, groups,
                                  detection_curves, detection_times; kwargs...)

Use the selected non-growing criteria to derive blank wells. Matching rows are
removed from the sample matrix and combined with any annotated blanks, then the
remaining samples are corrected within their own group. `detection_curves` may
be a smoothed version of `curves`; the unsmoothed rows are used for the blank
trace itself.
"""
function derive_blank_from_non_growing(
    curves::Matrix{Float64},
    times::Vector{Float64},
    labels::Vector{String},
    groups::Vector{String},
    detection_curves::Matrix{Float64},
    detection_times::Vector{Float64};
    annotated_blank_curves::Vector{Vector{Float64}}=Vector{Vector{Float64}}(),
    annotated_blank_times::Vector{Vector{Float64}}=Vector{Vector{Float64}}(),
    annotated_blank_groups::Vector{String}=String[],
    annotated_blank_labels::Vector{String}=String[],
    prescreen_constant::Bool=false,
    trend_test::Bool=false,
    prescreen_tol::Float64=1.5,
    prescreen_q_low::Float64=0.05,
    prescreen_q_high::Float64=0.95,
    trend_p_threshold::Float64=0.05,
    blank_method::Symbol=:pointbypoint,
    blank_floor::Union{Nothing,Float64}=1e-4,
)
    size(curves, 1) == length(labels) == length(groups) || throw(ArgumentError(
        "curves, labels, and groups must contain the same number of samples"
    ))
    size(detection_curves, 1) == size(curves, 1) || throw(ArgumentError(
        "detection_curves must contain one row per sample"
    ))
    derived_idx = detect_non_growing_indices(
        detection_curves, detection_times;
        prescreen_constant, trend_test, prescreen_tol,
        prescreen_q_low, prescreen_q_high, trend_p_threshold,
    )
    keep = setdiff(1:size(curves, 1), derived_idx)

    derived_curves = [Vector(curves[i, :]) for i in derived_idx]
    all_blank_curves = vcat(annotated_blank_curves, derived_curves)
    all_blank_times = vcat(annotated_blank_times, [copy(times) for _ in derived_idx])
    all_blank_groups = vcat(annotated_blank_groups, groups[derived_idx])
    all_blank_labels = vcat(annotated_blank_labels, labels[derived_idx])

    remaining = [Vector(curves[i, :]) for i in keep]
    corrected, corrected_mask = apply_grouped_blank_subtraction(
        remaining, [copy(times) for _ in keep], groups[keep],
        all_blank_curves, all_blank_times, all_blank_groups;
        method=blank_method,
        floor=blank_floor,
    )
    corrected_matrix = isempty(corrected) ? zeros(Float64, 0, size(curves, 2)) :
                       reduce(vcat, permutedims.(corrected))
    return (
        curves=corrected_matrix,
        labels=labels[keep],
        groups=groups[keep],
        blank_labels=all_blank_labels,
        derived_indices=derived_idx,
        corrected_mask=corrected_mask,
    )
end

"""
    derive_blank_from_non_growing_sources(source_curves, source_times,
        source_labels, source_groups, detection_curves, detection_times,
        detection_labels; kwargs...)

Detect non-growing rows on a prepared clustering grid, then map the selected
labels back to their original source curves. Annotated and derived blanks are
combined and applied to the remaining source curves before interpolation or
irregular-grid construction. This keeps measurement-index blank subtraction
on one shared source schedule instead of mixing raw and resampled indices.
"""
function derive_blank_from_non_growing_sources(
    source_curves::Vector{Vector{Float64}},
    source_times::Vector{Vector{Float64}},
    source_labels::Vector{String},
    source_groups::Vector{String},
    detection_curves::Matrix{Float64},
    detection_times::Vector{Float64},
    detection_labels::Vector{String};
    annotated_blank_curves::Vector{Vector{Float64}}=Vector{Vector{Float64}}(),
    annotated_blank_times::Vector{Vector{Float64}}=Vector{Vector{Float64}}(),
    annotated_blank_groups::Vector{String}=String[],
    annotated_blank_labels::Vector{String}=String[],
    prescreen_constant::Bool=false,
    trend_test::Bool=false,
    prescreen_tol::Float64=1.5,
    prescreen_q_low::Float64=0.05,
    prescreen_q_high::Float64=0.95,
    trend_p_threshold::Float64=0.05,
    blank_method::Symbol=:pointbypoint,
    blank_floor::Union{Nothing,Float64}=1e-4,
)
    length(source_curves) == length(source_times) ==
        length(source_labels) == length(source_groups) || throw(ArgumentError(
        "source curves, times, labels, and groups must have the same length"
    ))
    size(detection_curves, 1) == length(detection_labels) || throw(ArgumentError(
        "detection_curves must contain one row per detection label"
    ))
    length(unique(source_labels)) == length(source_labels) || throw(ArgumentError(
        "source labels must be unique"
    ))
    length(unique(detection_labels)) == length(detection_labels) || throw(ArgumentError(
        "detection labels must be unique"
    ))
    length(annotated_blank_curves) == length(annotated_blank_times) ==
        length(annotated_blank_groups) == length(annotated_blank_labels) ||
        throw(ArgumentError(
            "annotated blank curves, times, groups, and labels must have the same length"
        ))

    source_index = Dict(label => i for (i, label) in enumerate(source_labels))
    missing_labels = filter(label -> !haskey(source_index, label), detection_labels)
    isempty(missing_labels) || throw(ArgumentError(
        "detection labels missing from source data: $(join(missing_labels, ", "))"
    ))

    derived_idx = detect_non_growing_indices(
        detection_curves, detection_times;
        prescreen_constant, trend_test, prescreen_tol,
        prescreen_q_low, prescreen_q_high, trend_p_threshold,
    )
    keep_idx = setdiff(eachindex(detection_labels), derived_idx)
    derived_labels = detection_labels[derived_idx]
    remaining_labels = detection_labels[keep_idx]
    derived_source_idx = [source_index[label] for label in derived_labels]
    remaining_source_idx = [source_index[label] for label in remaining_labels]

    all_blank_curves = vcat(
        annotated_blank_curves,
        source_curves[derived_source_idx],
    )
    all_blank_times = vcat(
        annotated_blank_times,
        source_times[derived_source_idx],
    )
    all_blank_groups = vcat(
        annotated_blank_groups,
        source_groups[derived_source_idx],
    )
    all_blank_labels = vcat(annotated_blank_labels, derived_labels)

    corrected, corrected_mask = apply_grouped_blank_subtraction(
        source_curves[remaining_source_idx],
        source_times[remaining_source_idx],
        source_groups[remaining_source_idx],
        all_blank_curves,
        all_blank_times,
        all_blank_groups;
        method=blank_method,
        floor=blank_floor,
    )
    return (
        curves=corrected,
        times=source_times[remaining_source_idx],
        labels=remaining_labels,
        groups=source_groups[remaining_source_idx],
        blank_labels=all_blank_labels,
        derived_indices=derived_idx,
        corrected_mask=corrected_mask,
    )
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
    groups = fill("uploaded_matrix", length(labels))
    return [times_raw for _ in labels], curves_all, labels, groups,
           Vector{Vector{Float64}}(), Vector{Vector{Float64}}(), String[], String[]
end

function _annotation_well_sets(path::String)
    excluded = Set{String}()
    blanks = Set{String}()
    isfile(path) || return excluded, blanks
    ann = CSV.read(path, DataFrame; header=false, silencewarnings=true, stringtype=String)
    for i in 1:nrow(ann)
        ncol(ann) >= 2 || continue
        marker_raw = ann[i, 2]
        marker = ismissing(marker_raw) ? "X" : string(marker_raw)
        marker in ("", "missing") && (marker = "X")
        well = string(ann[i, 1])
        marker in ("b", "X", "x") && push!(excluded, well)
        marker == "b" && push!(blanks, well)
    end
    return excluded, blanks
end

function _experiment_annotation_path(
    experiment_dir::String,
    channel::Int,
    channel_annotation::Bool,
)::String
    if !channel_annotation
        path = joinpath(experiment_dir, "annotation_clean.csv")
        isfile(path) || error("Annotation file not found: $path")
        return path
    end

    files = try
        readdir(experiment_dir)
    catch
        String[]
    end
    prefix = "annotation_channel_$(channel)_"
    candidates = sort(filter(
        f -> startswith(f, prefix) && endswith(f, ".csv"),
        files,
    ))
    !isempty(candidates) && return joinpath(experiment_dir, first(candidates))

    has_channel_annotations = any(
        f -> occursin(r"^annotation_channel_\d+_", f),
        files,
    )
    if !has_channel_annotations
        fallback = joinpath(experiment_dir, "annotation_clean.csv")
        isfile(fallback) && return fallback
    end
    error("No annotation file found for channel $channel in $experiment_dir")
end

function _blank_summary(
    blank_curves::Vector{Vector{Float64}},
    n_timepoints::Int,
)::Tuple{Float64,Vector{Float64}}
    isempty(blank_curves) && return 0.0, zeros(Float64, n_timepoints)

    values = filter(isfinite, reduce(vcat, blank_curves))
    blank_value = isempty(values) ? 0.0 : mean(values)
    blank_timeseries = Vector{Float64}(undef, n_timepoints)
    for j in 1:n_timepoints
        at_time = Float64[]
        for curve in blank_curves
            j <= length(curve) || continue
            isfinite(curve[j]) && push!(at_time, curve[j])
        end
        blank_timeseries[j] = isempty(at_time) ? NaN : mean(at_time)
    end
    valid = filter(isfinite, blank_timeseries)
    fallback = isempty(valid) ? 0.0 : mean(valid)
    replace!(blank_timeseries, NaN => fallback)
    return blank_value, blank_timeseries
end

"""
    load_experiment_data(clean_data_path, experiment;
                             channel=1, channel_annotation=false)

Load one experiment from GUIbiont's `Clean_data` layout without using any
browser-produced result. `clean_data_path` is the user-selected local path,
and `experiment` identifies its subdirectory.

The returned named tuple contains:
- `data`: raw nonblank, non-discarded wells as a [`GrowthData`](@ref);
- `blank_value`: the mean of all finite readings in wells annotated `"b"`;
- `blank_timeseries`: the pointwise mean of those annotated blank wells;
- `blank_labels` and `excluded_labels`.

Set `channel_annotation=true` for GUIbiont replicate selections: the loader
then uses `annotation_channel_N_*.csv` for the selected channel, falling back
to `annotation_clean.csv` only when no channel-specific annotations exist.
This helper is intended for exported workflows, which can recompute blank
correction and fitting directly from source files.
"""
function load_experiment_data(
    clean_data_path::String,
    experiment::String;
    channel::Int=1,
    channel_annotation::Bool=false,
)
    experiment_dir = joinpath(clean_data_path, experiment)
    data_file = joinpath(experiment_dir, "data_channel_$(channel).csv")
    isfile(data_file) || error("Data file not found: $data_file")
    annotation_file = _experiment_annotation_path(
        experiment_dir,
        channel,
        channel_annotation,
    )

    raw = CSV.read(data_file, DataFrame; header=1, silencewarnings=true)
    names_raw = String.(names(raw))
    length(names_raw) >= 2 || error("Data file must contain time and at least one well")
    time_data = raw[!, names(raw)[1]]
    times = if nonmissingtype(eltype(time_data)) <: AbstractString
        parsed = [tryparse(Float64, string(t)) for t in time_data[2:end]]
        any(==(nothing), parsed) ?
            Float64.(0:(nrow(raw) - 1)) :
            [0.0; Float64.(parsed)]
    else
        numeric = _to_f64_vector(time_data)
        any(isnan, numeric) ? Float64.(0:(nrow(raw) - 1)) : numeric
    end

    annotation = CSV.read(
        annotation_file,
        DataFrame;
        header=false,
        silencewarnings=true,
        stringtype=String,
    )
    excluded = Set{String}()
    blanks = Set{String}()
    for i in 1:nrow(annotation)
        ncol(annotation) >= 2 || continue
        well = string(annotation[i, 1])
        marker_raw = annotation[i, 2]
        marker = ismissing(marker_raw) ? "X" : string(marker_raw)
        marker in ("", "missing") && (marker = "X")
        marker in ("b", "X", "x") && push!(excluded, well)
        marker == "b" && push!(blanks, well)
    end

    sample_labels = String[]
    sample_curves = Vector{Vector{Float64}}()
    blank_labels = String[]
    blank_curves = Vector{Vector{Float64}}()
    for well in names_raw[2:end]
        curve = _to_f64_vector(raw[!, well])
        if well in blanks
            push!(blank_labels, well)
            push!(blank_curves, curve)
        elseif !(well in excluded)
            push!(sample_labels, well)
            push!(sample_curves, curve)
        end
    end
    isempty(sample_curves) && error("No nonblank sample wells found in $data_file")

    curves = reduce(vcat, permutedims.(sample_curves))
    blank_value, blank_timeseries = _blank_summary(
        blank_curves,
        length(times),
    )
    return (
        data=GrowthData(curves, times, sample_labels),
        blank_value=blank_value,
        blank_timeseries=blank_timeseries,
        blank_labels=sort(blank_labels),
        excluded_labels=sort(collect(excluded)),
    )
end

function _load_clustering_experiments(clean_data_path::String, experiments::Vector{String})
    times_all = Vector{Vector{Float64}}()
    curves_all = Vector{Vector{Float64}}()
    labels = String[]
    groups = String[]
    blank_curves_all = Vector{Vector{Float64}}()
    blank_times_all = Vector{Vector{Float64}}()
    blank_labels = String[]
    blank_groups = String[]
    for exp_name in experiments
        data_file = joinpath(clean_data_path, exp_name, "data_channel_1.csv")
        annotation_file = joinpath(clean_data_path, exp_name, "annotation_clean.csv")
        isfile(data_file) || continue
        gd_raw = CSV.read(data_file, DataFrame; header=1, silencewarnings=true)
        excluded_wells, blank_wells = _annotation_well_sets(annotation_file)
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
            if well in excluded_wells
                if well in blank_wells
                    push!(blank_curves_all, od)
                    push!(blank_times_all, time_numeric)
                    push!(blank_labels, "$(exp_name)/$(well)")
                    push!(blank_groups, exp_name)
                end
            else
                push!(times_all, time_numeric)
                push!(curves_all, od)
                push!(labels, "$(exp_name)/$(well)")
                push!(groups, exp_name)
            end
        end
    end
    return times_all, curves_all, labels, groups,
           blank_curves_all, blank_times_all, blank_labels, blank_groups
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
  (`"b"` marks blanks; `"X"`/`"x"` wells are excluded without being used as blanks).

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
    blank_floor::Union{Nothing,Float64} = 1e-4,
    blank_range_thr::Float64 = 0.005,
    blank_od_percentile::Float64 = 0.10,
    derive_non_growing_blanks::Bool = false,
    blank_prescreen_constant::Bool = false,
    blank_trend_test::Bool = false,
    blank_prescreen_tol::Float64 = 1.5,
    blank_prescreen_q_low::Float64 = 0.05,
    blank_prescreen_q_high::Float64 = 0.95,
    blank_trend_p_threshold::Float64 = 0.05,
    detection_smooth::Bool = false,
    detection_smooth_method::Symbol = :lowess,
    detection_smooth_pt_avg::Int = 7,
    detection_lowess_frac::Float64 = 0.05,
    detection_gaussian_h_mult::Float64 = 2.0,
)::GrowthData
    times_all, curves_all, labels, groups, blank_curves_all, blank_times_all,
    blank_labels, blank_groups = if !isempty(csv_path)
        _load_clustering_csv(csv_path)
    else
        isempty(clean_data_path) && error("clean_data_path is required when csv_path is empty")
        isempty(experiments) && error("experiments must not be empty when csv_path is empty")
        _load_clustering_experiments(clean_data_path, experiments)
    end
    isempty(curves_all) && error("No data loaded")
    source_times = times_all
    source_curves = curves_all
    source_labels = copy(labels)
    source_groups = copy(groups)
    used_irregular_grid = !interpolate && any(c -> any(isnan, c), curves_all)

    # Annotated blanks are paired within each experiment on the sample's real
    # time grid. This must happen before any cross-experiment resampling.
    annotated_blanks = !isempty(blank_curves_all)
    if subtract_blank && annotated_blanks && !derive_non_growing_blanks
        curves_all, _ = apply_grouped_blank_subtraction(
            curves_all, times_all, groups,
            blank_curves_all, blank_times_all, blank_groups;
            method=blank_method,
            floor=blank_floor,
        )
    end

    if interpolate
        times = common_time_grid(times_all; n_grid=max(10, interp_n),
                                 q_low=interp_quantile_lo, q_high=interp_quantile_hi)
        isempty(times) && error("Could not build interpolation grid")
        curves = interpolate_curves_to_grid(curves_all, times_all, times)
    elseif used_irregular_grid
        clean_times = Vector{Vector{Float64}}()
        clean_curves = Vector{Vector{Float64}}()
        clean_labels = String[]
        clean_groups = String[]
        for i in eachindex(curves_all)
            stripped = _strip_nan_tail_for_clustering(times_all[i], curves_all[i])
            stripped === nothing && continue
            ct, cy = stripped
            push!(clean_times, ct); push!(clean_curves, cy); push!(clean_labels, labels[i])
            push!(clean_groups, groups[i])
        end
        isempty(clean_curves) && error("No valid curves after NaN removal")
        igd = IrregularGrowthData(clean_curves, clean_times, clean_labels; step=0.01)
        times, curves, labels, groups = igd.times, igd.curves, clean_labels, clean_groups
    else
        min_len = minimum(length.(curves_all))
        times = times_all[1][1:min_len]
        curves = Matrix{Float64}(undef, length(curves_all), min_len)
        for i in eachindex(curves_all)
            curves[i, :] = curves_all[i][1:min_len]
        end
    end

    if derive_non_growing_blanks
        (blank_prescreen_constant || blank_trend_test) || error(
            "derive_non_growing_blanks requires blank_prescreen_constant or blank_trend_test"
        )
        detector_curves = fill_nonfinite_colmean(curves)
        detector_times = times
        if detection_smooth
            detector_data = preprocess(
                GrowthData(detector_curves, times, labels),
                FitOptions(
                    smooth=true,
                    smooth_method=detection_smooth_method,
                    smooth_pt_avg=detection_smooth_pt_avg,
                    lowess_frac=detection_lowess_frac,
                    gaussian_h_mult=detection_gaussian_h_mult,
                    cluster=false,
                ),
            )
            detector_curves = detector_data.curves
            detector_times = detector_data.times
        end
        derived = derive_blank_from_non_growing_sources(
            source_curves, source_times, source_labels, source_groups,
            detector_curves, detector_times, labels;
            annotated_blank_curves=blank_curves_all,
            annotated_blank_times=blank_times_all,
            annotated_blank_groups=blank_groups,
            annotated_blank_labels=blank_labels,
            prescreen_constant=blank_prescreen_constant,
            trend_test=blank_trend_test,
            prescreen_tol=blank_prescreen_tol,
            prescreen_q_low=blank_prescreen_q_low,
            prescreen_q_high=blank_prescreen_q_high,
            trend_p_threshold=blank_trend_p_threshold,
            blank_method,
            blank_floor,
        )
        isempty(derived.derived_indices) && error("No non-growing curves were detected for blank derivation")
        isempty(derived.labels) && error("No sample curves remain after blank derivation")
        if interpolate
            curves = interpolate_curves_to_grid(derived.curves, derived.times, times)
            labels, groups = derived.labels, derived.groups
        elseif used_irregular_grid
            clean_times = Vector{Vector{Float64}}()
            clean_curves = Vector{Vector{Float64}}()
            clean_labels = String[]
            clean_groups = String[]
            for i in eachindex(derived.curves)
                stripped = _strip_nan_tail_for_clustering(
                    derived.times[i], derived.curves[i],
                )
                stripped === nothing && continue
                ct, cy = stripped
                push!(clean_times, ct)
                push!(clean_curves, cy)
                push!(clean_labels, derived.labels[i])
                push!(clean_groups, derived.groups[i])
            end
            isempty(clean_curves) && error("No valid curves remain after blank derivation")
            igd = IrregularGrowthData(clean_curves, clean_times, clean_labels; step=0.01)
            times, curves = igd.times, igd.curves
            labels, groups = clean_labels, clean_groups
        else
            ncols = length(times)
            curves = Matrix{Float64}(undef, length(derived.curves), ncols)
            for i in eachindex(derived.curves)
                curves[i, :] = derived.curves[i][1:ncols]
            end
            labels, groups = derived.labels, derived.groups
        end
    end

    if !derive_non_growing_blanks && !annotated_blanks && auto_detect_blanks
        blank_idxs = Int[]
        for group in unique(groups)
            group_rows = findall(==(group), groups)
            local_idxs = detect_blank_indices(curves[group_rows, :], times;
                flat_range_thr=blank_range_thr, od_percentile=blank_od_percentile)
            append!(blank_idxs, group_rows[local_idxs])
        end
        if !isempty(blank_idxs)
            blank_curves_all = [Vector(curves[i, :]) for i in blank_idxs]
            blank_labels     = labels[blank_idxs]
            blank_times_all = [copy(times) for _ in blank_idxs]
            blank_groups = groups[blank_idxs]
            keep = setdiff(1:size(curves, 1), blank_idxs)
            isempty(keep) && error("No sample curves remain after blank removal")
            curves, labels, groups = curves[keep, :], labels[keep], groups[keep]
        end
    end

    # Auto-detected blanks already share the prepared grid, but remain grouped
    # so one experiment can never supply the correction for another.
    if !derive_non_growing_blanks && subtract_blank && !annotated_blanks && !isempty(blank_curves_all)
        curve_vectors = [Vector(curves[i, :]) for i in axes(curves, 1)]
        corrected, _ = apply_grouped_blank_subtraction(
            curve_vectors, [copy(times) for _ in curve_vectors], groups,
            blank_curves_all, blank_times_all, blank_groups;
            method=blank_method,
            floor=blank_floor,
        )
        curves = reduce(vcat, permutedims.(corrected))
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
export fill_nonfinite_colmean
export common_time_grid
export interpolate_curves_to_grid
export apply_grouped_blank_subtraction
export derive_blank_from_non_growing
export derive_blank_from_non_growing_sources
export cluster_quality_indices
export prepare_clustering_data
export load_experiment_data
