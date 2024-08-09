# API 

## Pre-processing functions
### Smoothing data
[`Kinbiont.smoothing_data`](@ref)

```@docs
Kinbiont.smoothing_data
```

### Correction for multiple scattering
[`Kinbiont.correction_OD_multiple_scattering`](@ref)

```@docs
Kinbiont.correction_OD_multiple_scattering
```

## Fitting one kinetics

### Log-Lin fitting
[`Kinbiont.fitting_one_well_Log_Lin`](@ref)

```@docs
Kinbiont.fitting_one_well_Log_Lin
```
### Analysis of segments
[`Kinbiont.segment_gr_analysis`](@ref)

```@docs
Kinbiont.segment_gr_analysis
```
### ODE fitting
#### Fitting a harcoded model

[`Kinbiont.fitting_one_well_ODE_constrained`](@ref)

```@docs
Kinbiont.fitting_one_well_ODE_constrained
```

#### Fitting a custom model

[`Kinbiont.fitting_one_well_custom_ODE`](@ref)

```@docs
Kinbiont.fitting_one_well_custom_ODE
```
#### ODE model selection

[`Kinbiont.ODE_Model_selection`](@ref)

```@docs
Kinbiont.ODE_Model_selection
```
#### ODE Morris sensitivity

[`Kinbiont.one_well_morris_sensitivity`](@ref)

```@docs
Kinbiont.one_well_morris_sensitivity
```

### NL Fitting


[`Kinbiont.fit_NL_model`](@ref)

```@docs
Kinbiont.fit_NL_model
```
[`Kinbiont.fit_NL_model_with_sensitivity`](@ref)

```@docs
Kinbiont.fit_NL_model_with_sensitivity
```

[`Kinbiont.fit_NL_model_bootstrap`](@ref)

```@docs
Kinbiont.fit_NL_model_bootstrap
```

#### NL Model selection
[`Kinbiont.NL_model_selection`](@ref)

```@docs
Kinbiont.NL_model_selection
```



### Segmented fitting 

#### ODE segmentation
[`Kinbiont.selection_ODE_fixed_intervals`](@ref)

```@docs
Kinbiont.selection_ODE_fixed_intervals
```
[`Kinbiont.segmentation_ODE`](@ref)

```@docs
Kinbiont.segmentation_ODE
```
#### NL segmentation
[`Kinbiont.selection_NL_fixed_interval`](@ref)

```@docs
Kinbiont.selection_NL_fixed_interval
```
[`Kinbiont.segmentation_NL`](@ref)

```@docs
Kinbiont.segmentation_NL
```


## Fitting one a .csv file

### Log-Lin fitting
[`Kinbiont.fit_one_file_Log_Lin`](@ref)

```@docs
Kinbiont.fit_one_file_Log_Lin
```

### Analysis of segments
[`Kinbiont.segment_gr_analysis_file`](@ref)

```@docs
Kinbiont.segment_gr_analysis_file
```


### ODE fitting

[`Kinbiont.fit_file_ODE`](@ref)

```@docs
Kinbiont.fit_file_ODE
```

[`Kinbiont.fit_file_custom_ODE`](@ref)

```@docs
Kinbiont.fit_file_custom_ODE
```

[`Kinbiont.ODE_model_selection_file`](@ref)

```@docs
Kinbiont.ODE_model_selection_file
```
### NL fitting
[`Kinbiont.fit_NL_model_file`](@ref)

```@docs
Kinbiont.fit_NL_model_file
```
[`Kinbiont.fit_NL_model_selection_file`](@ref)

```@docs
Kinbiont.fit_NL_model_selection_file
```
### Segmented fitting 
[`Kinbiont.segmentation_ODE_file`](@ref)

```@docs
Kinbiont.fit_NL_segmentation_file
```
[`Kinbiont.fit_NL_segmentation_file`](@ref)

```@docs
Kinbiont.segmentation_ODE_file
```

## ML downstrema analysis

[`Kinbiont.downstream_decision_tree_regression`](@ref)

```@docs
Kinbiont.downstream_decision_tree_regression
```
[`Kinbiont.downstream_symbolic_regression`](@ref)

```@docs
Kinbiont.downstream_symbolic_regression
```

## Simulations 
### ODE simulations
[`Kinbiont.ODE_sim`](@ref)

```@docs
Kinbiont.ODE_sim
```
### Stochastic simulations
[`Kinbiont.stochastic_sim`](@ref)

```@docs
Kinbiont.stochastic_sim
```
## Various
 [`Kinbiont.specific_gr_evaluation`](@ref)

```@docs
Kinbiont.specific_gr_evaluation
```



