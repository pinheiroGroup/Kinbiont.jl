# Kinbiont 

```@raw html
<div style="text-align: center; margin-bottom: 20px; margin: auto; max-width: 320px;">
    <img alt="Kinbiont logo" src="./assets/logo.png">
</div>
```

Kinbiont is a Julia package designed as an end-to-end pipeline for biological discovery, enabling data-driven generation of hypotheses that can be tested in targeted experiments. Leveraging Julia's speed, flexibility, and growing popularity in biological sciences, Kinbiont integrates advanced solvers for ordinary differential equations (ODEs), non-linear optimization methods, signal processing, and interpretable ML algorithms. Kinbiont can fit various models of microbial dynamics, including solutions in closed-form, ODEs, and user-defined models. Unlike existing tools, Kinbiont extends model-based parameter estimation to fits with segmentation, allowing for more detailed analysis of the microbial dynamics. The inferred parameters can then be analyzed using symbolic regression to discover equations that capture response patterns or decision tree to identify informative response patterns from large datasets.

```@contents
Pages = [
    "01_install/index.md",
    "02_data/index.md",
    "04_math/index.md",
    "05_examples_simulations/index.md",
    "06_examples_fit/index.md",
    "07_examples_ML/index.md",
    "api.md"
]
```
