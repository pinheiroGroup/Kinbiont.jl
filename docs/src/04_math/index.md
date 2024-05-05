# [The mathematical models](@id models)

1. [ODEs for bacterial growth](#ODE_list)
2. [ Stochastic models for bacterial growth](#stoch_model)
3. [Type of loss](#loss)

## ODEs for bacterial growth

TO DO

The models and their parameters are sumarrized in the following table
| Model                                  |  Parameters                                       |
| --------------------------------------- | -------------------------------------------------------- |
| Diauxic_piecewise_damped_logistic      | gr_1, N_max, shape_1, lag, linear_const, t_shift, gr_2, N_max_2, shape_2, end_second_lag, lag_2_gr |
| Diauxic_replicator_1                    | gr, N_max, lag, arbitrary_const, linear_const |
| Diauxic_replicator_2                    | gr, N_max, lag, arbitrary_const, linear_const, growth_stationary |
| HPM                                    | gr, exit_lag_rate, N_max |
| dHPM                                    | gr, exit_lag_rate, N_max,shape |
| HPM_3_death                            | gr, exit_lag_rate, inactivation_rate, death_rate |
| HPM_3_death_resistance                 | gr, exit_lag_rate, inactivation_rate, death_rate, n_res, n_max |
| HPM_3_inhibition                       | gr, exit_lag_rate, inactivation_rate |
| HPM_inhibition                         | gr, inhibition_rate, gr_inhibition, N_max |
| HPM_exp                                | gr, exit_lag_rate |
| ODEs_HPM_SR                            | gr, gr_phage, scale, death_rate, resistance_rate |
| baranyi_exp                            | gr, lag_time, shape |
| baranyi_richards                       | gr, N_max, lag_time, shape |
| baranyi_roberts                        | gr, N_max, lag_time, shape_1, shape_2 |
| bertalanffy_richards                   | gr, N_max, shape |
| exponential                             | gr |
| four_piecewise                          | gr, gr_2, gr_3, gr_4, lag, t_decay_gr, t_stationary |
| gbsm_piecewise                         | gr, a_1, b_1, c, a_2, b_2 |
| gompertz                               | gr, N_max |
| hyper_gompertz                         | gr, N_max, shape |
| hyper_logistic                         | doubling_time, gr, N_max, shape |
| huang                                  | gr, N_max, lag |
| logistic                               | gr, N_max |
| logistic                               | gr, N_max, shape |
| ode_von_bertalanffy                    | alpha, beta, a, b |
| piecewise_damped_logistic              | gr, N_max, lag, shape, linear_const |
| triple_piecewise                       | gr, gr_2, gr_3, lag, t_stationary |
| triple_piecewise_bertalanffy_richards  | gr, gr_lag, t_lag, t_stationary, gr_stat, shape, N_max |
| triple_piecewise_damped_logistic       | gr, gr_2, gr_3, lag, t_stationary, N_max |
| triple_piecewise_sublinear             | gr, gr_2, gr_3, lag, t_stationary, N_max |

## Stochastic models for bacterial growth


In the stochastic version of the growth models, the growth rate of each population component (denoted as $\mu_i$) is evaluated based on the concentration of the limiting nutrient. The user is required to specify the starting amount of nutrients and the volume of the solution. Various kinetic growth models are considered.

### Monod Model

The Monod model is described by the following equation:


$\mu([\text{Nut.}]; k_1, \mu_\text{max}) = \mu_\text{max} \frac{[\text{Nut.}]}{k_1 + [\text{Nut.}]}.$

### Haldane Model

The Haldane model is expressed as:


$\mu([\text{Nut.}]; k_1, k_2, \mu_\text{max}) = \mu_\text{max} \frac{[\text{Nut.}]}{k_1 + [\text{Nut.}] + \frac{k_2}{[\text{Nut.}]^2}}.$

### Blackman Model

The Blackman model is given by:


$\mu([\text{Nut.}]; k_1, \mu_\text{max}) = \mu_\text{max} \frac{[\text{Nut.}]}{k_1 + [\text{Nut.}]}.$

### Tessier Model

The Tessier model is represented as:


$\mu([\text{Nut.}]; k_1, \mu_\text{max}) = \mu_\text{max} (1 - e^{k_1[\text{Nut.}] }).$

### Moser Model

The Moser model is defined by:


$\mu([\text{Nut.}]; k_1, k_2, \mu_\text{max}) = \mu_\text{max} \frac{[\text{Nut.}]^{k_2}}{k_1 + [\text{Nut.}]^{k_2}}.$

### Aiba-Edwards Model

The Aiba-Edwards model is given by:

$\mu([\text{Nut.}]; k_1, k_2, \mu_\text{max}) = \mu_\text{max} \frac{[\text{Nut.}]}{k_1 + [\text{Nut.}]} e^{-\frac{[\text{Nut.}]}{k_2}}.$

### Verhulst Model

The Verhulst model is defined as:


$\mu(N; N_\text{max}, \mu_\text{max}) = \mu_\text{max} \left(1 - \frac{N}{N_\text{max}}\right).$

Where $[\text{Nut.}]$ is the limiting nutrient concentration, $\mu_\text{max}$ is the maximum possible growth rate, $k_1$ and $k_2$ are numerical constants with meanings depending


## Type of loss functions

`type_of_loss = "L2"`: Minimize the L2 norm of the difference between the numerical solution of an ODE and the given data.


 $\mathcal{D}(D(t_i), \bar{N}(t_i, \{P\})) = \left(D(t_i) - \bar{N}(t_i, \{P\})\right)^2$
 
`type_of_loss = "L2_derivative"`: Minimize the L2 norm of the difference between the derivatives of the numerical solution of an ODE and the corresponding derivatives of the data.

`type_of_loss ="RE" `: Minimize the relative error between the solution and data 

$\mathcal{D}(D(t_i), \bar{N}(t_i, \{P\})) = 0.5 \cdot \left(1 - \frac{D(t_i)}{\bar{N}(t_i, \{P\})}\right)^2$

`type_of_loss = "blank_weighted_L2"` : Minimize a weighted version of the L2 norm, where the difference between the solution and data is weighted based on a distribution obtained from empirical blank data. 

$\mathcal{D}(D(t_i), \bar{N}(t_i, \{P\})) = (1 - P(D(t_i) - \bar{N}(t_i, \{P\})|\text{noise})) \cdot \left(D(t_i) - \bar{N}(t_i, \{P\})\right)^2$


