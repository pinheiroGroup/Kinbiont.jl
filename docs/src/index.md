# Kinbiont 

```@raw html
<div style="text-align: center; margin-bottom: 20px; margin: auto; max-width: 320px;">
    <img alt="Kinbiont logo" src="./assets/logo.png">
</div>
```

Kinbiont is a Julia package designed as an end-to-end pipeline for biological discovery, enabling data-driven generation of hypotheses that can be tested in targeted experiments. Leveraging Julia's speed, flexibility, and growing popularity in biological, Kinbiont integrates advanced solvers for ordinary differential equations (ODEs), non-linear optimization methods, signal processing, and interpretable ML algorithms. Kinbiont can fit various models of microbial dynamics, including solutions in closed-form, ODEs, and user-defined models. Unlike existing tools, Kinbiont extends model-based parameter estimation to fits with segmentation, allowing for more detailed analysis of the microbial dynamics. The inferred parameters can then be analyzed using symbolic regression to discover equations that capture response patterns or decision tree to identify informative response patterns from large datasets.

New to Kinbiont? Start with the [Quick Start](@ref quickstart) — you can go from a CSV file to fitted growth parameters in under 10 lines of Julia.

```@contents
Pages = [
    "01_install/index.md",
    "02_quickstart/index.md",
    "03_data/index.md",
    "04_preprocessing/index.md",
    "05_fitting/index.md",
    "06_clustering/index.md",
    "07_ml/index.md",
    "08_simulations/index.md",
    "09_math/index.md",
    "api.md"
]
```
