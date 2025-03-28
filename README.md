# Kinbiont.jl

Ecological and evolutionary processes of microbes are characterized by observables like growth rates and biomass yield, inferred from kinetics experiments. 
Across conditions, these observables map response patterns such as antibiotic growth inhibition and yield dependence on substrate.
But how do we extract ecological and evolutionary insights from massive datasets of time-resolved microbial data? Here, we introduce Kinbiont — an ecosystem of numerical methods combining state-of-the-art solvers for ordinary differential equations, non-linear optimization, signal processing, and interpretable machine learning algorithms.
Kinbiont provides a comprehensive, model-based analysis pipeline, covering all aspects of microbial kinetics data, from preprocessing to result interpretation. 

# Documentation

For stable documentation please consult [Documentation](https://pinheirogroup.github.io/Kinbiont.jl/). 

# Publication

Pre-print at  [Biorxiv](https://www.biorxiv.org/content/10.1101/2024.09.09.611847v1).

Data and scripts to reproduce the paper results at [Kinbiont utilities](https://github.com/pinheiroGroup/Kinbiont_utilities)

## Recent Updates

### Improved Change Point Detection

The latest version includes a significantly faster change point detection algorithm based on binary segmentation with dynamic programming principles. Key improvements:

1. **Performance**: Up to 10-100x faster for large datasets compared to LSDD algorithm
2. **Scalability**: O(n log n) complexity instead of O(n²)
3. **Multiple cost functions**: Support for both L2 (mean-based) and L1 (median-based, robust to outliers) cost functions
4. **Early stopping**: Automatically stops when no significant change points are found

To use the improved algorithm:

```julia
# Fast binary segmentation with RSS cost function
change_points = cpd_local_detection(
    data_matrix,
    n_max_cp;
    type_of_detection="fast_bs",
    type_of_curve="original",
    cost_func=:rss
)

# Fast binary segmentation with L1 cost function (more robust to outliers)
change_points = cpd_local_detection(
    data_matrix,
    n_max_cp;
    type_of_detection="fast_bs",
    type_of_curve="original",
    cost_func=:l1
)
```

Example scripts are provided in the `Example_scripts` directory:
- `fast_change_point_detection.jl`: Simple example of using the algorithm
- `fast_change_point_benchmarks.jl`: Benchmarks comparing performance with the original algorithm
