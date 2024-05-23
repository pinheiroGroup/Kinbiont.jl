# API 

## Pre-processing functions
### Smoothing data
[`Kimchi.smoothing_data`](@ref)

```@docs
Kimchi.smoothing_data
```

### Correction for multiple scattering
[`Kimchi.correction_OD_multiple_scattering`](@ref)

```@docs
Kimchi.correction_OD_multiple_scattering
```

## Fitting one kinetics

### Log-Lin fitting
[`Kimchi.fitting_one_well_Log_Lin`](@ref)

```@docs
Kimchi.fitting_one_well_Log_Lin
```

### ODE fitting
#### Fitting a harcoded model

[`Kimchi.fitting_one_well_ODE_constrained`](@ref)

```@docs
Kimchi.fitting_one_well_ODE_constrained
```

#### Fitting a custom model

[`Kimchi.fitting_one_well_custom_ODE`](@ref)

```@docs
Kimchi.fitting_one_well_custom_ODE
```
#### ODE model selection

[`Kimchi.ODE_Model_selection`](@ref)

```@docs
Kimchi.ODE_Model_selection
```
#### ODE Morris sensitivity

[`Kimchi.one_well_morris_sensitivity`](@ref)

```@docs
Kimchi.one_well_morris_sensitivity
```

### NL Fitting

### NL Fitting

[`Kimchi.fit_NL_model`](@ref)

```@docs
Kimchi.fit_NL_model
```
[`Kimchi.fit_NL_model_with_sensitivity`](@ref)

```@docs
Kimchi.fit_NL_model_with_sensitivity
```
[`Kimchi.fit_NL_model_MCMC_intialization`](@ref)

```@docs
Kimchi.fit_NL_model_MCMC_intialization
```
[`Kimchi.fit_NL_model_bootstrap`](@ref)

```@docs
Kimchi.fit_NL_model_bootstrap
```
[`Kimchi.NL_error_blanks`](@ref)

```@docs
Kimchi.NL_error_blanks
```
#### NL Model selection
[`Kimchi.NL_model_selection`](@ref)

```@docs
Kimchi.NL_model_selection
```



### Segmented fitting 

#### ODE segmentation
[`Kimchi.selection_ODE_fixed_intervals`](@ref)

```@docs
Kimchi.selection_ODE_fixed_intervals
```
[`Kimchi.segmentation_ODE`](@ref)

```@docs
Kimchi.segmentation_ODE
```
#### NL segmentation
[`Kimchi.selection_NL_fixed_interval`](@ref)

```@docs
Kimchi.selection_NL_fixed_interval
```
[`Kimchi.selection_NL_max_change_points`](@ref)

```@docs
Kimchi.selection_NL_max_change_points
```

## Plot a file
[`Kimchi.plot_data`](@ref)

```@docs
Kimchi.plot_data
```

## Fitting one a .csv file

### Log-Lin fitting
[`Kimchi.fit_one_file_Log_Lin`](@ref)

```@docs
Kimchi.fit_one_file_Log_Lin
```
### ODE fitting

[`Kimchi.fit_file_ODE`](@ref)

```@docs
Kimchi.fit_file_ODE
```

[`Kimchi.fit_file_custom_ODE`](@ref)

```@docs
Kimchi.fit_file_custom_ODE
```

[`Kimchi.ODE_model_selection_file`](@ref)

```@docs
Kimchi.ODE_model_selection_file
```
### NL fitting
[`Kimchi.fit_NL_model_file`](@ref)

```@docs
Kimchi.fit_NL_model_file
```
[`Kimchi.fit_NL_model_selection_file`](@ref)

```@docs
Kimchi.fit_NL_model_selection_file
```
### Segmented fitting 
[`Kimchi.segmentation_ODE_file`](@ref)

```@docs
Kimchi.fit_NL_segmentation_file
```
[`Kimchi.fit_NL_segmentation_file`](@ref)

```@docs
Kimchi.segmentation_ODE_file
```
## Simulations 
### ODE simulations
[`Kimchi.ODE_sim`](@ref)

```@docs
Kimchi.ODE_sim
```
### Stochastic simulations
[`Kimchi.stochastic_sim`](@ref)

```@docs
Kimchi.stochastic_sim
```
## Various
 [`Kimchi.specific_gr_evaluation`](@ref)

```@docs
Kimchi.specific_gr_evaluation
```



