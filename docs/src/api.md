# API 

## Pre-processing functions
### Smoothing data
[`KinBiont.smoothing_data`](@ref)

```@docs
KinBiont.smoothing_data
```

### Correction for multiple scattering
[`KinBiont.correction_OD_multiple_scattering`](@ref)

```@docs
KinBiont.correction_OD_multiple_scattering
```

## Fitting one kinetics

### Log-Lin fitting
[`KinBiont.fitting_one_well_Log_Lin`](@ref)

```@docs
KinBiont.fitting_one_well_Log_Lin
```
### Analysis of segments
[`KinBiont.segment_gr_analysis`](@ref)

```@docs
KinBiont.segment_gr_analysis
```
### ODE fitting
#### Fitting a harcoded model

[`KinBiont.fitting_one_well_ODE_constrained`](@ref)

```@docs
KinBiont.fitting_one_well_ODE_constrained
```

#### Fitting a custom model

[`KinBiont.fitting_one_well_custom_ODE`](@ref)

```@docs
KinBiont.fitting_one_well_custom_ODE
```
#### ODE model selection

[`KinBiont.ODE_Model_selection`](@ref)

```@docs
KinBiont.ODE_Model_selection
```
#### ODE Morris sensitivity

[`KinBiont.one_well_morris_sensitivity`](@ref)

```@docs
KinBiont.one_well_morris_sensitivity
```

### NL Fitting


[`KinBiont.fit_NL_model`](@ref)

```@docs
KinBiont.fit_NL_model
```
[`KinBiont.fit_NL_model_with_sensitivity`](@ref)

```@docs
KinBiont.fit_NL_model_with_sensitivity
```

[`KinBiont.fit_NL_model_bootstrap`](@ref)

```@docs
KinBiont.fit_NL_model_bootstrap
```

#### NL Model selection
[`KinBiont.NL_model_selection`](@ref)

```@docs
KinBiont.NL_model_selection
```



### Segmented fitting 

#### ODE segmentation
[`KinBiont.selection_ODE_fixed_intervals`](@ref)

```@docs
KinBiont.selection_ODE_fixed_intervals
```
[`KinBiont.segmentation_ODE`](@ref)

```@docs
KinBiont.segmentation_ODE
```
#### NL segmentation
[`KinBiont.selection_NL_fixed_interval`](@ref)

```@docs
KinBiont.selection_NL_fixed_interval
```
[`KinBiont.segmentation_NL`](@ref)

```@docs
KinBiont.segmentation_NL
```


## Fitting one a .csv file

### Log-Lin fitting
[`KinBiont.fit_one_file_Log_Lin`](@ref)

```@docs
KinBiont.fit_one_file_Log_Lin
```

### Analysis of segments
[`KinBiont.segment_gr_analysis_file`](@ref)

```@docs
KinBiont.segment_gr_analysis_file
```


### ODE fitting

[`KinBiont.fit_file_ODE`](@ref)

```@docs
KinBiont.fit_file_ODE
```

[`KinBiont.fit_file_custom_ODE`](@ref)

```@docs
KinBiont.fit_file_custom_ODE
```

[`KinBiont.ODE_model_selection_file`](@ref)

```@docs
KinBiont.ODE_model_selection_file
```
### NL fitting
[`KinBiont.fit_NL_model_file`](@ref)

```@docs
KinBiont.fit_NL_model_file
```
[`KinBiont.fit_NL_model_selection_file`](@ref)

```@docs
KinBiont.fit_NL_model_selection_file
```
### Segmented fitting 
[`KinBiont.segmentation_ODE_file`](@ref)

```@docs
KinBiont.fit_NL_segmentation_file
```
[`KinBiont.fit_NL_segmentation_file`](@ref)

```@docs
KinBiont.segmentation_ODE_file
```

## ML downstrema analysis

[`KinBiont.downstream_decision_tree_regression`](@ref)

```@docs
KinBiont.downstream_decision_tree_regression
```
[`KinBiont.downstream_symbolic_regression`](@ref)

```@docs
KinBiont.downstream_symbolic_regression
```

## Simulations 
### ODE simulations
[`KinBiont.ODE_sim`](@ref)

```@docs
KinBiont.ODE_sim
```
### Stochastic simulations
[`KinBiont.stochastic_sim`](@ref)

```@docs
KinBiont.stochastic_sim
```
## Various
 [`KinBiont.specific_gr_evaluation`](@ref)

```@docs
KinBiont.specific_gr_evaluation
```



