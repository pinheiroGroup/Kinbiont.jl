using Kinbiont
using Plots
using CSV
using DataFrames

"""
Example script demonstrating the fast binary segmentation algorithm for change point detection.
This algorithm is significantly more efficient for large datasets than the existing methods.

The performance improvement comes from:
1. Using precomputed cumulative sums for efficient cost calculation
2. Binary segmentation approach with dynamic programming principles
3. Optional early stopping when no significant change points are found
4. Support for both L1 (median-based) and L2 (mean-based) cost functions
"""

# Load example data
path_to_data = joinpath(dirname(dirname(Base.active_project())), "data_examples", "plate_data.csv")
data = CSV.read(path_to_data, DataFrame)

# Select a specific well for change point detection
well_id = "A1"  # Change this to analyze a different well
well_data = data[data.well_id .== well_id, :]

# Extract time and OD values and convert to matrix format
time = collect(well_data.time)
od = collect(well_data.od_corrected)
data_matrix = Matrix(transpose(hcat(time, od)))

# Define parameters for change point detection
n_max_cp = 3          # Maximum number of change points to detect
min_segment_length = 3  # Minimum length of a segment

# Run change point detection with standard LSDD algorithm
lsdd_start = time_ns()
cp_lsdd = cpd_local_detection(
    data_matrix,
    n_max_cp;
    type_of_detection="lsdd",
    type_of_curve="original",
    method="peaks_prominence"
)
lsdd_time = (time_ns() - lsdd_start) / 1e9
println("LSDD algorithm time: $(lsdd_time) seconds")

# Run change point detection with fast binary segmentation using RSS cost
fast_bs_start = time_ns()
cp_fast_bs_rss = cpd_local_detection(
    data_matrix,
    n_max_cp;
    type_of_detection="fast_bs",
    type_of_curve="original",
    cost_func=:rss,
    min_segment_length=min_segment_length
)
fast_bs_rss_time = (time_ns() - fast_bs_start) / 1e9
println("Fast binary segmentation (RSS) time: $(fast_bs_rss_time) seconds")
println("Speedup: $(lsdd_time / fast_bs_rss_time)x")

# Run change point detection with fast binary segmentation using L1 cost (more robust to outliers)
fast_bs_l1_start = time_ns()
cp_fast_bs_l1 = cpd_local_detection(
    data_matrix,
    n_max_cp;
    type_of_detection="fast_bs",
    type_of_curve="original",
    cost_func=:l1,
    min_segment_length=min_segment_length
)
fast_bs_l1_time = (time_ns() - fast_bs_l1_start) / 1e9
println("Fast binary segmentation (L1) time: $(fast_bs_l1_time) seconds")

# Plot the results
p = plot(time, od, label="Data", linewidth=2, title="Change Point Detection Comparison (Well $well_id)")

# Add vertical lines for LSDD change points
for (i, t) in enumerate(cp_lsdd[2])
    plot!([t, t], [minimum(od), maximum(od)], label=i==1 ? "LSDD" : "", linestyle=:dash, color=:red)
end

# Add vertical lines for fast_bs (RSS) change points
for (i, t) in enumerate(cp_fast_bs_rss[2])
    plot!([t, t], [minimum(od), maximum(od)], label=i==1 ? "Fast BS (RSS)" : "", linestyle=:dot, color=:blue)
end

# Add vertical lines for fast_bs (L1) change points
for (i, t) in enumerate(cp_fast_bs_l1[2])
    plot!([t, t], [minimum(od), maximum(od)], label=i==1 ? "Fast BS (L1)" : "", linestyle=:dashdot, color=:green)
end

# Show the comparison of different methods
display(p)

# Optional: Save the plot
savefig(p, "change_point_detection_comparison.png")

# Print the detected change points for each method
println("\nChange points (LSDD): ", cp_lsdd[2])
println("Change points (Fast BS RSS): ", cp_fast_bs_rss[2])
println("Change points (Fast BS L1): ", cp_fast_bs_l1[2])