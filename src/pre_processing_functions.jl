using Lowess
using CSV
using DataFrames

"""
    correction_OD_multiple_scattering(
    data::Matrix{Float64},
    calibration_curve::String;
    method="interpolation"
    )

Multiple scattering correction of a given time series.
# Arguments:

- `data`: Matrix of size 2xn, where n is the number of time points (single curve).
- `calibration_OD_curve="NA"`: String. Path for the calibration data (.csv file). It is used only if `multiple_scattering_correction=true`.
- `method`: String. Method of choice to perform the multiple scattering curve inference. Options: '"interpolation"' or '"exp_fit"' (adapted from Meyers, A., Furtmann, C., & Jose, J., *Enzyme and microbial technology*, 118, 1-5., 2018). 

# Output: `Matrix{Float64}`. Array with the corrected data. 

"""
function correction_OD_multiple_scattering(
    data::Matrix{Float64},
    calibration_curve::String;
    method="interpolation"
)

    od_calib = CSV.File(calibration_curve)
    names_of_cols = propertynames(od_calib)

    if method == "interpolation"
        Od_real = od_calib[names_of_cols[1]]
        od_calib_array = Matrix(transpose(hcat(Od_real, od_calib[names_of_cols[2]])))
        soterd_calib = sort!(od_calib_array, rev=false, dims=2)
        itp =
            interpolate(soterd_calib[1, :], soterd_calib[2, :], SteffenMonotonicInterpolation())
        extrap_spline = extrapolate(itp, 0)
        corrected_data = [extrap_spline(k) for k in data[2, :]]
        data_fin = Matrix(transpose(hcat(data[1, :], corrected_data)))
        # exponential fit from "Direct optical density determination of 
        #bacterial cultures in microplates for high-throughput screening applications"
    elseif method == "exp_fit"
        NL_model_exp_ms(x, p) = p[1] .* (1 .- exp.(.-x .* p[2]))
        p0 = [0.5, 0.5]
        xdata = od_calib[names_of_cols[1]]
        ydata = od_calib[names_of_cols[2]]
        fit = LsqFit.curve_fit(NL_model_exp_ms, xdata, ydata, p0)
        corrected_data = NL_model_exp_ms(data[2, :], fit.param)
        data_fin = Matrix(transpose(hcat(data[1, :], corrected_data)))

    end


    return data_fin
end



"""
    smoothing_data(
        data::Matrix{Float64};
        method="rolling_avg",
        pt_avg=7,
        thr_lowess=0.05
    )

# Arguments:

- `data`: Matrix of size 2xn, where n is the number of time points (single curve).
- `method="rolling_avg"`: String. Method of choice to smooth the data. Options: "NO", "rolling_avg" (rolling average of the data), and "lowess".
- `pt_avg=7`: Number of points to generate the initial condition or to do the rolling avgerage smoothing.
- `thr_lowess=0.05`: Float64. Argument of the lowess smoothing.

# Output: `Matrix{Float64}`. Array of smoothed data. 

"""
function smoothing_data(
    data::Matrix{Float64};
    method="rolling_avg",
    pt_avg=7,
    thr_lowess=0.05
)

    
    if method == "rolling_avg" && pt_avg < 3
        println("WARNING: the number of points to do rolling average is too low")
        println("changing the method of smoothing to lowess")
        method = "lowess"
    end

    if method == "rolling_avg"
        times = [
            sum(@view data[1, i:(i+pt_avg-1)]) / pt_avg for
            i = 1:(eachindex(data[1, :])[end]-(pt_avg-1))
        ]
        values = [
            sum(@view data[2, i:(i+pt_avg-1)]) / pt_avg for
            i = 1:(eachindex(data[2, :])[end]-(pt_avg-1))
        ]
        smoothed_data = Matrix(transpose(hcat(times, values)))
    elseif method == "lowess"

        model_fit = lowess_model(data[1, :], data[2, :], thr_lowess)
        smoothed_data = Matrix(transpose(hcat(data[1, :], model_fit)))

    else

        #  println(" Warning wrong smoothing input, this part is skipped")

        smoothed_data = copy(data)

    end
    return smoothed_data
end

function blank_subtraction(
    dfs_data::Any, # dataset first row times second row OD ,
    list_of_blank::Any;
    method="avg_blank"
)

    if method == "avg_blank"

        blank_values = mean([mean(dfs_data[k]) for k in list_of_blank])

    elseif method == "time_blank"

        blank_values =
            [mean([dfs_data[k][j] for k in list_of_blank]) for j in eachindex(times_data)]

    else

        blank_values = zeros(length(dfs_data[list_of_blank[1]]))

    end

    return blank_values
end

function average_replicate(dfs_data, times_data, properties_of_annotation, names_of_annotated_df)

    new_data = times_data

    list_replicate = unique(properties_of_annotation)
    list_replicate = filter!(e -> e != "b", list_replicate)
    list_replicate = filter!(e -> e != "X", list_replicate)


    for replicate_temp in list_replicate

        names_of_replicate_temp =
            Symbol.(
                names_of_annotated_df[findall(
                    x -> x == replicate_temp,
                    properties_of_annotation,
                )]
            )
        replicate_mean = [
            mean([dfs_data[k][j] for k in names_of_replicate_temp]) for
            j in eachindex(times_data)
        ]

        new_data = hcat(new_data, replicate_mean)

    end

    new_data = DataFrame(new_data, :auto)
    rename!(new_data, vcat(:Time, reduce(vcat, Symbol.(list_replicate))))
    names_of_cols = propertynames(new_data)
    dfs_data = new_data

    return dfs_data, names_of_cols


end

function remove_negative_value(data::Any)

    # index bigger than zero
    index_not_zero = findall(>(0), data)
    data = data[index_not_zero]

    return data, index_not_zero
end

function negative_value_correction(data::Any,
    list_of_blank::Any;
    method="remove",
    thr_negative=0.01,)


    if method == "blank_correction"

        times = data[1, :]
        values = data[2, :]
        index_neg = findall(x -> x < 0.0, values)
        number_of_neg = length(index_neg)
        replacement_values = abs.(sample(list_of_blank, number_of_neg) .- mean(list_of_blank))
        values[index_neg] .= replacement_values
        data_corrected = transpose(hcat(times, values))

    elseif method == "thr_correction"
        times = data[1, :]
        values = data[2, :]
        index_neg = findall(x -> x < thr_negative, values)
        values[index_neg] .= thr_negative
        data_corrected = transpose(hcat(times, values))
    else
        times = data[1, :]
        values = data[2, :]
        values_corrected, index_not_zero = remove_negative_value(values)
        data_corrected = transpose(hcat(times[index_not_zero], values_corrected))

    end

    return Matrix(data_corrected)
end

export correction_OD_multiple_scattering
export smoothing_data
export blank_subtraction
export average_replicate
export remove_negative_value
export negative_value_correction
