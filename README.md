<p align="center">
  <img src="https://github.com/pinheiroGroup/Kinbiont.jl/blob/main/static/kinbiont_logo.png">
</p>

Kinbiont (a Julia package for Model-based Analyses of microbial Kinetics) is a versatile software tool that utilizes Ordinary Differential Equations (ODEs) to fit bacterial growth data from plate reader experiments. 
With Kinbiont it is possible to simulate, fit, perform model selection, and conduct sensitivity analysis for multi-well plate reader experiments.
The parameter fitting in Kinbiont is defined as a constrained optimization problem, which is solved using a differential evolution box-constrained non-linear optimizer.

To address complex cases, Kinbiont preprocesses data by detecting change points in the differential operator (or in the solution of the ODE) within the time series. 
It then automatically assembles and fits a segmented ODE model, resulting in a fully interpretable representation of the time series.
