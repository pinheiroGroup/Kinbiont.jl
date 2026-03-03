include("gaussian_smoothing_utils.jl")
include("time_series_clustering_utils.jl")

using .timeSeriesKmeansLib

function Kinbiont_data_preprocessing(Kinbiont_data, Kinbiont_options_preprocessing)
  # data and options packages renaming
  opt_pack = Kinbiont_options_preprocessing
  data_pack = Kinbiont_data

  # unpack data
  data = data_pack.Data_matrix 
  data_time = data_pack.Times
  Labels = data_pack.Labels # time series labels

  # unpack preprocessing options
  smooth_data = opt_pack.smooth_data
  avg_replicate = opt_pack.avg_replicate
  negative_corrections = opt_pack.negative_corrections
  multiple_scattering_correction = opt_pack.multiple_scattering_correction
  calibration_OD_curve = opt_pack.calibration_OD_curve
  do_blank_subtraction = opt_pack.do_blank_subtraction

  # clustering parameters: 
  cluster_data = opt_pack.cluster_data # boolean (do clustering or not)
  if cluster_data
    K = opt_pack.K === nothing ? 2 : opt_pack.K # number of clusters
    # robust minimum:
    q_low = opt_pack.q_low === nothing ? 0.0 : opt_pack.q_low
    # robust maximum:  
    q_high = opt_pack.q_high === nothing ? 1.0 : opt_pack.q_high
    # constant curves threshold: ex.: maximum <= 2.0 * minimum
    tol_const = opt_pack.tol_const === nothing ? 2.0 : opt_pack.tol_const

    #= threshold for fraction of inbound / out-of-bound points for constant
      cluster bias for fair comparison against other more fitting clusters =#
    const_qbias = opt_pack.const_qbias === nothing ? 0.7 :
      opt_pack.const_qbias

    # trend factor for bias correction, used to obtain more growing curves
    trend_factor = opt_pack.trend_factor === nothing ? 0.5 :
      opt_pack.trend_factor
    # boolean: find exponential cluster or not
    use_exp_clust = opt_pack.use_exp_clust === nothing ? true : 
      opt_pack.use_exp_clust
    # clustering rng for initial random points
    clustering_rng = opt_pack.rng === nothing ? Random.MersenneTwister(42) :
      opt_pack.rng
  end

  # Step 1 multiple scattering correction
  if multiple_scattering_correction
    data = correction_OD_multiple_scattering_matrix(
        data,
        calibration_OD_curve;
    )
  end

  # Step 3.2 average replicates
  if avg_replicate
    data, Labels = average_replicate_matrix(data, Labels)
  end

  # Step 1 smoothing
  if smooth_data
    data, data_time = smoothing_data_matrix(data, data_time, 
      Kinbiont_options_preprocessing)
  end

  # === Blank Correction ===
  # Step 3.3 blank computation  ONLY MEAN VALUES!!!
  if do_blank_subtraction in ("yes", "auto")
    mean_blank = blank_mean_computation(data, Labels)
    # Step 3.4 blank subtraction SUBTRACT MEAN VALUES TO ALL THE DATA!!!
    data .-= mean_blank
  end
  # === End Blank Correction ===

  # Step 3.5 negative value removal FOR NOT blanks
  if negative_corrections !== nothing
    for (i, label) in enumerate(Labels)
        if label ∉ ["X", "b"]
            data[i, data[i, :] .<= 0.0] .= negative_corrections
        end
    end
  end
  
  # Step ? clustering 
  if cluster_data
    X = Array{Float64}(data)
    t = Float64.(data_time)
    clustering_labels, wcss, centroids_orig, centroids_norm =
        curves_clustering(
            X, # data
            t, # time points
            Int(K); # K is generic, put to Int for safety

            # robust minimum and maximum using quantiles
            q_low = q_low,
            q_high = q_high,
            tol_const = tol_const,

            # paramters for bias computation for constant cluster
            const_qbias = const_qbias,
            trend_factor = trend_factor,

            use_exp_cluster = use_exp_clust, # find exponential cluster
            rng = clustering_rng # rng value for random initialization
        )
    # TO DO: modify Kinbiont_data to add these parameters
    # save clustering_labels (vector of cluster id for each curve):
    Kinbiont_data.Clustering_labels = clustering_labels
    # save wcss result from clustering:
    Kinbiont_data.Clustering_wcss = wcss
    #= save z-norm centroids as computed by k-means (see this as a 
       K-length vector of vectors (K prototypes time points)) =#
    Kinbiont_data.Clustering_centroids_norm = centroids_norm
    # save original space centroids (see this as a vector of vectors)
    Kinbiont_data.Clustering_centroids_orig = centroids_orig
  end

  # wrap the data 
  Kinbiont_data.Data_matrix = data
  Kinbiont_data.Times = data_time 
  Kinbiont_data.Labels = Labels

  return Kinbiont_data
end

function blank_mean_computation(data, labels)
    index_blank = findall(==("b"), labels)
    if isempty(index_blank)
        @warn "No blanks ('b') found. Skipping blank correction (mean_blank = 0.0)."
        return 0.0
    end
    mean_blank = mean(mean(data[index_blank, :], dims=1))
    return mean_blank
end

# not tested old
#=
function correction_OD_multiple_scattering_matrix(
  data::Matrix{Float64},
  calibration_curve::String;
)
    od_calib = CSV.File(calibration_curve)
    names_of_cols = propertynames(od_calib)

    Od_real = od_calib[names_of_cols[1]]
    od_calib_array = Matrix(transpose(hcat(Od_real, od_calib[names_of_cols[2]])))
    soterd_calib = sort!(od_calib_array, rev=false, dims=2)
    itp = interpolate(soterd_calib[1, :], soterd_calib[2, :], SteffenMonotonicInterpolation())
    extrap_spline = extrapolate(itp, 0)

    # recumpute data matrix
    corrected_data = [extrap_spline(k) for k in data]
    return corrected_data
end
=# 

# not tested (For Fabrizio: remade accordingly to safe issues and CSV.File documentation)
function correction_OD_multiple_scattering_matrix(
  data::Matrix{Float64},
  # or data::AbstactMatrix{<:Real}?
  calibration_curve::String; 
  # or calibration_curve::AbstractString?
)
  # Read calibration curve (table)
  tbl = CSV.File(calibration_curve)
  colnames = propertynames(tbl)  # is it (:OD_measured, :OD_real)?

  # Extact first two columns as vectors
  x_measured = collect(getproperty(tbl, colnames[1]))  # OD misurata
  y_real = collect(getproperty(tbl, colnames[2]))  # OD "vera"

  # Construct an x_measured-ordered calibration...
  idx = sortperm(x_measured)
  # ... And apply it (maintains same order for safety)
  x_sorted = x_measured[idx]
  y_sorted = y_real[idx]

  # Monotonic interpolation via Steffen
  itp = interpolate(x_sorted, y_sorted, SteffenMonotonicInterpolation())
  # Extrapolate
  extrap_spline = extrapolate(itp, 0)

  # Apply element-wise correction
  corrected_data = extrap_spline.(data)

  return corrected_data
end

# tested
function average_replicate_matrix(data, labels)
  labels_str = String.(labels)

  # keep only real replicate labels (exclude blanks and controls)
  reps = filter(x -> x ∉ ["X", "b"], unique(labels_str))

  nrep  = length(reps)
  ncols = size(data, 2)
  out = Array{Float64}(undef, nrep, ncols)
  out_names = String[]

  for i in eachindex(reps)
    rep = reps[i]
    idx = findall(labels_str .== rep)  # rows that belong to this replicate
    out[i, :] = mean(data[idx, :], dims=1)[:]  # column-wise mean
    push!(out_names, rep)
  end

  return out, out_names
end

function smoothing_data_matrix(
    data::AbstractMatrix{<:Real},
    data_time::AbstractVector{<:Real},
    Kinbiont_options_preprocessing
  )::Tuple{Matrix{Float64}, Vector{Float64}}

  # Materialize inputs as Float64 for internal computations
  time_vec  = Float64.(data_time)
  data_mat  = Array{Float64}(data)

  # unpack options
  type_of_smoothing = Kinbiont_options_preprocessing.type_of_smoothing

  # default: same grid, copy of data
  smoothed_data_time   = copy(time_vec)
  smoothed_data_matrix = similar(data_mat)

  # === Boxcar Filter ===
  if type_of_smoothing == "boxcar"
    # Window length (fallback order: pt_boxcar -> 5)
    w = 5
    try
      if Kinbiont_options_preprocessing.pt_boxcar !== nothing
        w = Kinbiont_options_preprocessing.pt_boxcar
      end
    catch
      @warn "WARNING: Missing Kinbiont_options_preprocessing.pt_boxcar value, 
        fallback to w = 5"
    end
    w = Int(w)

    # Keeps original time grid
    smoothed_data_time = time_vec

    # If w is too small, just copy data
    if w < 2
      @warn "WARNING: user has put w < 2 (too small), returning the same data"
      smoothed_data_matrix = copy(data_mat)
    else
      half = w ÷ 2 # the symbol means "integer division"
      nrows, ncols = size(data_mat)
      smoothed_data_matrix = Matrix{Float64}(undef, nrows, ncols)
      for i in 1:nrows
        curve = @view data_mat[i, :]
        for j in 1:ncols
          left  = max(1, j - half)
          right = min(ncols, j + half)
          smoothed_data_matrix[i, j] = mean(@view curve[left:right])
        end
    end
  end
  # === End Boxcar Filter ===

  # === Rolling Average Filter ===
  elseif type_of_smoothing == "rolling_avg"
    # Get the window pts
    pt_avg = Kinbiont_options_preprocessing.pt_avg
    if pt_avg < 3
      @warn "WARNING: the number of points to do rolling average is " *
        "too low, falling back to lowess method..."
      type_of_smoothing = "lowess"
    else
      # perform rolling average smoothing of times
      smoothed_data_time = [
        sum(@view time_vec[i:(i+pt_avg-1)]) / pt_avg
            for i in 1:(eachindex(time_vec)[end] - (pt_avg-1))
      ]

      # perform rolling average smoothing of the matrix
      nrows, ncols = size(data_mat)
      smoothed_data_matrix = Matrix{Float64}(undef, nrows, ncols)
      for i in 1:nrows
        smoothed_data_matrix[i, :] = [
          sum(@view data_mat[i, j:(j+pt_avg-1)]) / pt_avg
            for j in 1:(eachindex(data_mat[i, :])[end] - (pt_avg-1))
        ]
      end
    end
  end
  # === End Rolling Average Filter ===

  # === Lowess Filter ===
  # === !!! REMEMBER TO ALWAYS LEAVE THIS AFTER ROLLING !!! ===
  if type_of_smoothing == "lowess"
    thr_lowess = Kinbiont_options_preprocessing.thr_lowess
    smoothed_data_time   = time_vec
    nrows, ncols         = size(data_mat)
    smoothed_data_matrix = Matrix{Float64}(undef, nrows, ncols)

    for i in axes(data_mat, 1)
      curve = @view data_mat[i, :]
      model_fit = lowess_model(time_vec, curve, thr_lowess)
      smoothed_data_matrix[i, :] = model_fit
    end
  # === End Lowess Filter ===

  # === Gaussian Filter ===
  # By specifying the target time_grid (smoothed_data_time)
  # the user can interpolate new data points (at the specified target times)
  elseif type_of_smoothing == "gaussian"
    # options for parameters
    h_mult    = Kinbiont_options_preprocessing.gaussian_h_mult
    time_grid = Kinbiont_options_preprocessing.gaussian_time_grid

    # new time grid for interpolation, if nothing --> same as time_vec
    smoothed_data_time = time_grid === nothing ? time_vec : Float64.(time_grid)

    nrows, _ = size(data_mat)
    smoothed_data_matrix = Matrix{Float64}(undef, nrows, length(smoothed_data_time))

    for i in 1:nrows
      y = @view data_mat[i, :]
      y_smoothed = gaussian_smoothing(
        time_vec,             # original time points (Float64)
        collect(y),           # raw curve as Vector{Float64}
        smoothed_data_time;   # target time grid
        h_mult = h_mult
      )
      smoothed_data_matrix[i, :] = y_smoothed
    end 
  # === End Gaussian Filter ===

  # === Safety Else ===
  else # check this else, might do that tomorrow
    @warn "Unknown type_of_smoothing = '$type_of_smoothing'. Returning original data."
    smoothed_data_matrix = copy(data_mat)
    smoothed_data_time   = time_vec
  end
  # === ---------- ====

  return smoothed_data_matrix, smoothed_data_time

end
