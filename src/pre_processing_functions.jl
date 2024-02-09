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



function smoothing_data(
    data::Matrix{Float64};
    method="rolling_avg",
    pt_avg=7,
    thr_lowess=0.05
)
    """
    rolling average smoothing of the data
    data  = matrix of data
    pt_avg = size of the windows of the rolling average
    """
    if  method == "rolling_avg" && pt_avg < 3
       println("WARNING: the number of points to do rolling average is too low")
       println("changing the method of smoothinf to lowess")
       method = "lowess"
    end

    if method == "rolling_avg"
        times = [
            sum(@view data[1, i:(i+pt_avg-1)]) / pt_avg for
            i = 1:(length(data[1, :])-(pt_avg-1))
        ]
        values = [
            sum(@view data[2, i:(i+pt_avg-1)]) / pt_avg for
            i = 1:(length(data[2, :])-(pt_avg-1))
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




function thr_negative_correction(
    data::Matrix{Float64}, # dataset first row times second row OD ,
    thr_negative::Float64, # the value at which data are setted if negative
)

    times = data[1, :]
    values = data[2, :]
    index_neg = findall(x -> x < thr_negative, values)
    values[index_neg] .= thr_negative
    data_corrected = transpose(hcat(times, values))

    return data_corrected
end

function blank_distrib_negative_correction(
    data::Matrix{Float64}, # dataset first row times second row OD ,
    blank_array::Vector{Float64},
)

    times = data[1, :]
    values = data[2, :]
    #println(values)
    index_neg = findall(x -> x < 0.0, values)

    number_of_neg = length(index_neg)
    replacement_values = abs.(sample(blank_array, number_of_neg) .- mean(blank_array))
    values[index_neg] .= replacement_values
    data_corrected = transpose(hcat(times, values))

    return data_corrected
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
            [mean([dfs_data[k][j] for k in list_of_blank]) for j = 1:length(times_data)]

    else

        blank_values = zeros(length(dfs_data[list_of_blank[1]]))

    end

    return blank_values
end


function average_replicate(dfs_data, times_data, properties_of_annotation, names_of_annotated_df)

    new_data = times_data

    list_replicate = unique(properties_of_annotation)
    list_replicate = filter!(e -> e != "b", list_replicate)


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
            j = 1:length(times_data)
        ]

        new_data = hcat(new_data, replicate_mean)

    end

    new_data = DataFrame(new_data, :auto)
    rename!(new_data, vcat(:Time, reduce(vcat, Symbol.(list_replicate))))
    names_of_cols = propertynames(new_data)
    dfs_data = new_data

    return dfs_data, names_of_cols


end

