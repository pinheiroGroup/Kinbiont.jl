# Kimchi - Julia Model-based Analyses of microbial Kinetics

```@raw html
<div style="text-align: center; margin-bottom: 20px;">
    <img alt="Kimchi logo" src="./assets/logo.png">
</div>
```

Kimchi (a Julia package for Model-based Analyses of microbial Kinetics) is a versatile software tool that utilizes Ordinary Differential Equations (ODEs) to fit bacterial growth data from plate reader experiments. 
With Kimchi it is possible to simulate, fit, perform model selection, and conduct sensitivity analysis for multi-well plate reader experiments.
The parameter fitting in Kimchi is defined as a constrained optimization problem, which is solved using a differential evolution box-constrained non-linear optimizer.

To address complex cases,  Kimchi preprocesses data by detecting change points in the differential operator (or in the solution of the ODE) within the time series. 
It then automatically assembles and fits a segmented ODE model, resulting in a fully interpretable representation of the time series.

```@contents
Pages = [
    "01_install/index.md",
    "02_data/index.md",
    "04_math/index.md",
    "05_examples/index.md"
]
```
