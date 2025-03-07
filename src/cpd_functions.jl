using Combinatorics
using Peaks
using ChangePointDetection
using Statistics
using LinearAlgebra

function getpoints_mod(profile; number_of_bin=100)
    points_list = ()
    # list of points
    min_profile = minimum(profile)
    max_profile = maximum(profile)
    delta_profile = (max_profile - min_profile) / number_of_bin
    seq_thr = min_profile:delta_profile:max_profile

    for thr_temp in seq_thr[2:end]
        points = Int[]
        exceeded = false

        for (index, value) in enumerate(profile)
            if (value > thr_temp) && exceeded == false
                push!(points, index)
                exceeded = true
            elseif (value < thr_temp) && exceeded == true
                exceeded = false
            end
        end

        if thr_temp == seq_thr[2]
            points_list = [points]
        else
            points_list = push!(points_list, points)
        end
    end

    return reverse(points_list[1:(end-1)])
end

"""
    fast_binary_segmentation(data::Matrix{Float64}, n_max_cp::Int; cost_func=:rss, min_segment_length=3)

Efficient change point detection using binary segmentation with dynamic programming.

# Arguments:
- `data::Matrix{Float64}`: A matrix containing the data for change point detection. First row is time, second row is values.
- `n_max_cp::Int`: Maximum number of change points to detect.

# Key Arguments:
- `cost_func=:rss`: Cost function to use for segmentation (:rss for residual sum of squares, :l1 for L1 norm)
- `min_segment_length=3`: Minimum length of segments to consider

# Output:
- A tuple containing the indices of detected change points and their corresponding times and values.
"""
function fast_binary_segmentation(
    data::Matrix{Float64}, 
    n_max_cp::Int; 
    cost_func=:rss,
    min_segment_length=3
)
    # Extract time and value vectors for cleaner code
    times = data[1, :]
    values = data[2, :]
    n = length(values)
    
    # Initialize change point storage
    change_points = Int[]
    change_point_scores = Float64[]
    
    # Precompute cumulative sums for efficient cost calculation
    sums = cumsum(values)
    sq_sums = cumsum(values.^2)
    
    # Function to calculate cost of a segment
    function segment_cost(start_idx, end_idx)
        if cost_func == :rss
            # RSS cost - efficient implementation using precomputed sums
            segment_length = end_idx - start_idx + 1
            
            if segment_length < min_segment_length
                return Inf
            end
            
            # Calculate segment mean
            segment_sum = start_idx > 1 ? sums[end_idx] - sums[start_idx-1] : sums[end_idx]
            segment_mean = segment_sum / segment_length
            
            # Calculate RSS
            segment_sq_sum = start_idx > 1 ? sq_sums[end_idx] - sq_sums[start_idx-1] : sq_sums[end_idx]
            return segment_sq_sum - 2 * segment_mean * segment_sum + segment_length * segment_mean^2
        elseif cost_func == :l1
            # L1 cost - median-based
            if end_idx - start_idx + 1 < min_segment_length
                return Inf
            end
            
            segment_median = median(values[start_idx:end_idx])
            return sum(abs.(values[start_idx:end_idx] .- segment_median))
        end
    end
    
    # Initial split evaluation over the whole sequence
    function evaluate_split(start_idx, end_idx)
        best_score = -Inf
        best_point = -1
        
        # Try all possible split points
        for split_idx in (start_idx+min_segment_length-1):(end_idx-min_segment_length)
            # Cost without split
            cost_whole = segment_cost(start_idx, end_idx)
            
            # Cost with split
            cost_left = segment_cost(start_idx, split_idx)
            cost_right = segment_cost(split_idx+1, end_idx)
            cost_split = cost_left + cost_right
            
            # Improvement score from splitting
            score = cost_whole - cost_split
            
            if score > best_score
                best_score = score
                best_point = split_idx
            end
        end
        
        return best_point, best_score
    end
    
    # Recursively apply binary segmentation
    segments = [(1, n)]
    
    for _ in 1:n_max_cp
        best_split_score = -Inf
        best_split_point = -1
        best_segment_idx = -1
        
        # Evaluate splits for each current segment
        for (i, (start_idx, end_idx)) in enumerate(segments)
            if end_idx - start_idx + 1 >= 2 * min_segment_length
                split_point, split_score = evaluate_split(start_idx, end_idx)
                
                if split_score > best_split_score
                    best_split_score = split_score
                    best_split_point = split_point
                    best_segment_idx = i
                end
            end
        end
        
        # If we can't find a meaningful split, stop
        if best_split_score <= 0 || best_split_point == -1
            break
        end
        
        # Add the change point
        push!(change_points, best_split_point)
        push!(change_point_scores, best_split_score)
        
        # Update segments
        start_idx, end_idx = segments[best_segment_idx]
        segments[best_segment_idx] = (start_idx, best_split_point)
        insert!(segments, best_segment_idx + 1, (best_split_point + 1, end_idx))
    end
    
    # Sort change points
    sorted_indices = sortperm(change_points)
    sorted_change_points = change_points[sorted_indices]
    
    # Get times and values at change points
    times_at_change_points = times[sorted_change_points]
    values_at_change_points = values[sorted_change_points]
    
    return sorted_change_points, times_at_change_points, values_at_change_points
end

"""
    cpd_local_detection(
    data::Matrix{Float64},
    n_max_cp::Int;
    type_of_detection="lsdd",
    type_of_curve="original",
    pt_derivative=0,
    size_win=2,
    method="peaks_prominence",
    number_of_bin=40,
    cost_func=:rss,
    min_segment_length=3
    )

Performs change point detection analysis using the specified algorithm on the provided data matrix.

# Arguments:

- `data::Matrix{Float64}`: A matrix containing the data for change point detection. Each row should correspond to a different time series or experimental condition.
- `n_max_cp::Int`: Maximum number of change points to detect.

# Key Arguments:

- `type_of_detection="lsdd"`: Algorithm for detecting change points. Options include:
  - `"lsdd"`: Least Squares Density Difference.
  - `"piecewise"`: Piecewise linear fitting on the specific growth rate (requires `pt_derivative` > 0).
  - `"fast_bs"`: Fast binary segmentation algorithm (much more efficient for large datasets).
- `type_of_curve="original"`: Specifies the curve on which to perform change point detection. Options include:
  - `"original"`: Use the original time series data.
  - `"deriv"`: Use the specific growth rate time series.
- `pt_derivative=0`: Number of points used to evaluate the derivative or specific growth rate:
  - `0`: Use numerical derivative.
  - `>0`: Use a sliding window approach with the specified window size.
- `size_win::Int`: Size of the window used for the change point detection algorithms.
- `method="peaks_prominence"`: Method for peak detection on the dissimilarity curve. Options include:
  - `"peaks_prominence"`: Orders peaks by prominence.
  - `"thr_scan"`: Uses a threshold to select peaks.
- `number_of_bin=40`: Number of bins used to generate the threshold if `method="thr_scan"`.
- `cost_func=:rss`: Cost function to use for fast_bs algorithm. Options include:
  - `:rss`: Residual sum of squares (mean-based cost function).
  - `:l1`: L1 norm (median-based cost function, more robust to outliers).
- `min_segment_length=3`: Minimum segment length for fast_bs algorithm.

# Output:

- An array containing the list of detected change points.

"""
function cpd_local_detection(
    data::Matrix{Float64},
    n_max_cp::Int;
    type_of_detection="lsdd",
    type_of_curve="original",
    pt_derivative=0,
    size_win=2,
    method="peaks_prominence",
    number_of_bin=40,
    cost_func=:rss,
    min_segment_length=3
)
    if type_of_detection == "fast_bs"
        # Prepare data for binary segmentation based on curve type
        if type_of_curve == "deriv"
            deriv = specific_gr_evaluation(data, pt_derivative)
            specific_gr_times = [
                (data[1, r] + data[1, (r+pt_derivative)]) / 2 for
                r = 1:1:(eachindex(data[2, :])[end].-pt_derivative)
            ]
            input_data = Matrix(transpose(hcat(specific_gr_times, deriv)))
        else
            input_data = data
        end
        
        # Run the fast binary segmentation algorithm
        list_of_cpds = fast_binary_segmentation(
            input_data,
            n_max_cp;
            cost_func=cost_func,
            min_segment_length=min_segment_length
        )
    elseif type_of_detection == "lsdd" && type_of_curve == "deriv"
        list_of_cpds = cpd_lsdd_profile(
            data,
            n_max_cp;
            type_of_curve="deriv",
            pt_deriv=pt_derivative,
            window_size=size_win,
            method=method,
            number_of_bin=number_of_bin,
        )
    elseif type_of_detection == "lsdd" && type_of_curve != "deriv"
        list_of_cpds = cpd_lsdd_profile(
            data,
            n_max_cp;
            pt_deriv=pt_derivative,
            window_size=size_win,
            type_of_curve="original",
            method=method,
            number_of_bin=number_of_bin,
        )
    else
        list_of_cpds = detect_list_change_points(
            data,
            n_max_cp;
            pt_deriv=pt_derivative,
            type_of_curve=type_of_curve,
            win_size=size_win,
            method=method,
            number_of_bin=number_of_bin,
        )
    end

    return list_of_cpds
end

function cpd_lsdd_profile(
    data::Matrix{Float64},
    n_max::Int;
    window_size=2,
    type_of_curve="original",
    pt_deriv=0,
    method="peaks_prominence",
    number_of_bin=40,
)

    selected_change_point_index = Any
    # evaluating the profile of lsdd on the data or on the derivative of the data
    if type_of_curve == "deriv"
        deriv = specific_gr_evaluation(data, pt_deriv)
        profile = ChangePointDetection.lsdd_profile(deriv; window=window_size)
    else
        profile = ChangePointDetection.lsdd_profile(data[2, :]; window=window_size)
    end

    # adding time to profile
    profile = convert.(Float64, profile)
    data_dissim = Matrix(transpose(hcat(data[1, 1:length(profile)], profile)))
    selected_change_point_index =
        peaks_detection(data_dissim, n_max; method=method, number_of_bin=number_of_bin)

    return selected_change_point_index
end


function detect_list_change_points(
    data::Matrix{Float64},
    n_max::Int;
    win_size=2,
    method="peaks_prominence",
    number_of_bin=40,
    type_of_curve="original",
    pt_deriv=7,
)

    if type_of_curve == "deriv"
        data_gr = specific_gr_evaluation(data, pt_deriv)
        specific_gr_times = [
            (data[1, r] + data[1, (r+pt_deriv)]) / 2 for
            r = 1:1:(eachindex(data[2, :])[end].-pt_deriv)
        ]
        data = Matrix(transpose(hcat(specific_gr_times,data_gr)))
    end

    curve_dissimilitary_deriv = curve_dissimilitary_lin_fitting(
        data, # dataset first row times second row OD
        1, # index of start
        win_size, # size sliding window
    )
    data_dissim = Matrix(
        transpose(
            hcat(
                data[1, convert.(Int, curve_dissimilitary_deriv[1, :])],
                curve_dissimilitary_deriv[2, :],
            ),
        ),
    )
    selected_change_point_index =
        peaks_detection(data_dissim, n_max; method=method, number_of_bin=number_of_bin)

    return selected_change_point_index
end


function peaks_detection(
    data::Matrix{Float64},
    n_max::Int;
    method="peaks_prominence",
    number_of_bin=40,
)
    if method == "peaks_prominence"
        index_of_peaks = findmaxima(data[2, :]; strict=true)
        array_prominence = peakproms(index_of_peaks[1], data[2, :])[2]
        index_prominence = peakproms(index_of_peaks[1], data[2, :])[1]

        if length(array_prominence) < n_max
            @warn "The max number of peaks is too much"
            top_prominence = sort(array_prominence)
        else
            top_prominence = sort(array_prominence)[((end-n_max)+1):end]
        end

        index_top_peaks = [findall(array_prominence .== i)[1] for i in top_prominence]
        selected_change_point_index = index_prominence[index_top_peaks]
        times_top_peaks = data[1, selected_change_point_index]
        values_top_peaks = data[2, selected_change_point_index]
    end

    if method == "thr_scan"
        if n_max == 1
            selected_change_point_index =
                argmax(data[2, :])
        else
            selected_change_point_list =
                getpoints_mod(data[2, :], number_of_bin=number_of_bin)
            lenght_cpd_list = length.(selected_change_point_list)

            if n_max > maximum(lenght_cpd_list)
               
                    @warn    "This number of peaks is to much selecting the max number detected"
                
                selected_change_point_index = selected_change_point_list[end]
            else
                selected_change_point_index =
                    selected_change_point_list[maximum(findlast(lenght_cpd_list .<= n_max))]

                if length(selected_change_point_index) != n_max
                    
                    @warn    "This number of peaks is not detected changing to nearest one smaller"
                    
                end
            end
        end
        times_top_peaks = data[1, selected_change_point_index]
        values_top_peaks = data[2, selected_change_point_index]
    end

    return selected_change_point_index, times_top_peaks, values_top_peaks
end

function curve_dissimilitary_lin_fitting(
    data::Matrix{Float64}, # dataset first row times second row OD
    start_time_Index::Int,
    size_wind::Int, # size sliding window
)
    discrepancy_measure_curve = [start_time_Index, 0.0]
    ending = convert(Int, (length(data[2, :]) - floor(size_wind / 2) * 2))

    for index_t = start_time_Index:ending
        # defining the window
        middle_index = convert(Int, (index_t + floor(size_wind / 2)))
        end_index = convert(Int, (index_t + floor(size_wind / 2) * 2))
        win_1_data = data[2, index_t:middle_index]
        win_2_data = data[2, middle_index:end_index]
        win_tot_data = data[2, index_t:end_index]
        win_1_times = data[1, index_t:middle_index]
        win_2_times = data[1, middle_index:end_index]
        win_tot_times = data[1, index_t:end_index]

        #fitting total data
        data_total = Matrix(transpose(hcat(win_tot_times, win_tot_data)))
        data_1 = Matrix(transpose(hcat(win_1_times, win_1_data)))
        data_2 = Matrix(transpose(hcat(win_2_times, win_2_data)))

        X_total = data_total[1, :]
        Y_total = data_total[2, :]
        X_1 = data_1[1, :]
        Y_1 = data_2[2, :]
        X_2 = data_1[1, :]
        Y_2 = data_2[2, :]
        N_1 = length(data_1[1, :])
        N_2 = length(data_2[1, :])
        N_tot = length(data_total[1, :])
        M_1 = [ones(N_1) X_1]
        M_2 = [ones(N_2) X_2]
        M_tot = [ones(N_tot) X_total]
        fit_1 = M_1 \ Y_1
        fit_2 = M_2 \ Y_2
        fit_total = M_tot \ Y_total




        # residual calculation
        res_total = sum([
            abs((
                data_total[2, ll] - fit_total[2] * data_total[1, ll] -
                fit_total[1]
            )) for ll = 1:length(data_total[1, :])
        ])


        # residual calculation
        res_win_1 = sum([
            abs((data_1[2, ll] - fit_1[2] * data_1[1, ll] - fit_1[1])) for
            ll = 1:length(data_1[1, :])
        ])

        #fitting win 2

        # residual calculation
        res_win_2 = sum([
            abs((data_2[2, ll] - fit_2[2] * data_2[1, ll] - fit_2[1])) for
            ll = 1:length(data_2[1, :])
        ])

        #evaluation of the cost
        cost = -res_total + res_win_1 + res_win_2
        discrepancy_measure_curve =
            hcat(discrepancy_measure_curve, [index_t + floor(size_wind / 2), cost])
        # stop when first change point is fitted
    end

    return discrepancy_measure_curve
end

function generation_of_combination_of_cpds(
    cp_list;
    n_fix=0
)
    if n_fix > length(cp_list)
        n_fix = length(cp_list) 
    end

    if n_fix == 0 #return al possible combinations 

        combinations_tot = collect(combinations(cp_list))

    else #return   combinations  with fixed lenght
        combinations_tot = collect(combinations(cp_list, n_fix))
    end



    return combinations_tot
end

export getpoints_mod
export cpd_local_detection
export cpd_lsdd_profile
export detect_list_change_points
export peaks_detection
export curve_dissimilitary_lin_fitting
export generation_of_combination_of_cpds
export fast_binary_segmentation
