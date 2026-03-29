using Kinbiont
using Random
using Distributions
using BenchmarkTools
using Plots

"""
This script benchmarks the performance of different change point detection algorithms
in KinBiont.jl, focusing on comparing the original LSDD algorithm with the new
fast binary segmentation algorithm.

The benchmark creates synthetic time series data with known change points and
measures the execution time and accuracy of each algorithm.
"""

# Set a random seed for reproducibility
Random.seed!(123)

# Function to generate synthetic data with change points
function generate_synthetic_data(n_points, change_points; noise_level=0.1)
    # Initialize time and values
    time = collect(1:n_points)
    values = zeros(n_points)
    
    # Initial segment
    slope = rand(Uniform(0.01, 0.05))
    intercept = rand(Uniform(0.1, 0.5))
    
    # Generate piecewise linear data with change points
    current_start = 1
    for cp in sort(change_points)
        # Fill segment with linear trend
        for i in current_start:cp
            values[i] = intercept + slope * (i - current_start + 1) + noise_level * randn()
        end
        
        # Change parameters for next segment
        intercept = values[cp]
        slope = rand(Uniform(-0.05, 0.05))
        current_start = cp + 1
    end
    
    # Fill the last segment
    for i in current_start:n_points
        values[i] = intercept + slope * (i - current_start + 1) + noise_level * randn()
    end
    
    # Return as matrix
    return Matrix(transpose(hcat(time, values))), sort(change_points)
end

# Benchmark settings
n_datasets = 5
dataset_sizes = [100, 1000, 5000, 10000]  # Different data sizes to test scalability
n_change_points = 3

# Results containers
lsdd_times = zeros(length(dataset_sizes))
fast_bs_times = zeros(length(dataset_sizes))
lsdd_accuracies = zeros(length(dataset_sizes))
fast_bs_accuracies = zeros(length(dataset_sizes))

# Function to calculate accuracy (measured by mean distance to true change points)
function calc_accuracy(detected_points, true_points)
    if isempty(detected_points) || isempty(true_points)
        return 0.0
    end
    
    # Pair each detected point with the closest true point
    total_dist = 0.0
    for dp in detected_points
        min_dist = minimum(abs.(dp .- true_points))
        total_dist += min_dist
    end
    
    # Return inverse of mean distance (higher is better)
    return 1.0 / (1.0 + total_dist / length(detected_points))
end

# Run benchmarks for each dataset size
for (i, size) in enumerate(dataset_sizes)
    println("\nTesting with dataset size: $size")
    
    lsdd_time_sum = 0.0
    fast_bs_time_sum = 0.0
    lsdd_accuracy_sum = 0.0
    fast_bs_accuracy_sum = 0.0
    
    for j in 1:n_datasets
        # Generate synthetic data with known change points
        true_change_points = sort(sample(2:size-1, n_change_points, replace=false))
        data, true_cp = generate_synthetic_data(size, true_change_points)
        
        # LSDD algorithm
        lsdd_result = @timed cpd_local_detection(
            data, 
            n_change_points, 
            type_of_detection="lsdd",
            type_of_curve="original"
        )
        lsdd_cp = lsdd_result[1][1]  # Extract change point indices
        lsdd_time = lsdd_result[2]   # Extract execution time
        
        # Fast binary segmentation
        fast_bs_result = @timed cpd_local_detection(
            data, 
            n_change_points, 
            type_of_detection="fast_bs",
            type_of_curve="original"
        )
        fast_bs_cp = fast_bs_result[1][1]  # Extract change point indices
        fast_bs_time = fast_bs_result[2]   # Extract execution time
        
        # Calculate accuracy
        lsdd_acc = calc_accuracy(lsdd_cp, true_cp)
        fast_bs_acc = calc_accuracy(fast_bs_cp, true_cp)
        
        # Accumulate results
        lsdd_time_sum += lsdd_time
        fast_bs_time_sum += fast_bs_time
        lsdd_accuracy_sum += lsdd_acc
        fast_bs_accuracy_sum += fast_bs_acc
        
        println("  Dataset $j: LSDD time: $(round(lsdd_time, digits=4))s, Fast BS time: $(round(fast_bs_time, digits=4))s")
        println("  Dataset $j: LSDD accuracy: $(round(lsdd_acc, digits=4)), Fast BS accuracy: $(round(fast_bs_acc, digits=4))")
    end
    
    # Average results
    lsdd_times[i] = lsdd_time_sum / n_datasets
    fast_bs_times[i] = fast_bs_time_sum / n_datasets
    lsdd_accuracies[i] = lsdd_accuracy_sum / n_datasets
    fast_bs_accuracies[i] = fast_bs_accuracy_sum / n_datasets
    
    println("Average for size $size:")
    println("  LSDD time: $(round(lsdd_times[i], digits=4))s")
    println("  Fast BS time: $(round(fast_bs_times[i], digits=4))s")
    println("  Speed improvement: $(round(lsdd_times[i] / fast_bs_times[i], digits=2))x")
    println("  LSDD accuracy: $(round(lsdd_accuracies[i], digits=4))")
    println("  Fast BS accuracy: $(round(fast_bs_accuracies[i], digits=4))")
end

# Plot time comparison
p1 = plot(
    dataset_sizes, 
    [lsdd_times fast_bs_times], 
    label=["LSDD" "Fast Binary Segmentation"],
    linewidth=2,
    xlabel="Dataset Size",
    ylabel="Execution Time (s)",
    title="Algorithm Performance Comparison",
    markers=[:circle :square],
    yscale=:log10
)

# Plot accuracy comparison
p2 = plot(
    dataset_sizes, 
    [lsdd_accuracies fast_bs_accuracies], 
    label=["LSDD" "Fast Binary Segmentation"],
    linewidth=2,
    xlabel="Dataset Size",
    ylabel="Accuracy",
    title="Algorithm Accuracy Comparison",
    markers=[:circle :square],
    ylims=(0, 1)
)

# Plot speedup
speedup = lsdd_times ./ fast_bs_times
p3 = plot(
    dataset_sizes, 
    speedup, 
    linewidth=2,
    xlabel="Dataset Size",
    ylabel="Speedup Factor",
    title="Fast BS Speedup vs LSDD",
    label="Speedup",
    marker=:diamond
)

# Combine plots
combined_plot = plot(p1, p2, p3, layout=(3,1), size=(800, 800))
display(combined_plot)
savefig(combined_plot, "change_point_algorithm_comparison.png")

println("\nBenchmark Summary:")
println("Dataset Sizes: $dataset_sizes")
println("LSDD Times: $(round.(lsdd_times, digits=4))")
println("Fast BS Times: $(round.(fast_bs_times, digits=4))")
println("Speedup: $(round.(speedup, digits=2))")
println("LSDD Accuracies: $(round.(lsdd_accuracies, digits=4))")
println("Fast BS Accuracies: $(round.(fast_bs_accuracies, digits=4))")