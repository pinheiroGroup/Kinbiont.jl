# [The mathematical models](@id models)

1. [NL models for bacterial growth](#NL_list)
2. [ODEs for bacterial growth](#ODE_list)
3. [ Stochastic models for bacterial growth](#stoch_model)
4. [Error functions](#loss)
## NL models for bacterial growth



KinBiont employs the following NL model to fit growth curves


- **Exponential**

$$N(t) = p_1 \cdot e^{p_2 \cdot t}$$

- **Gompertz**

$$N(t) = p_1 \cdot e^{-e^{-p_2 \cdot (t - p_3)}}$$

- **Logistic**

$$N(t) = \frac{p_3}{1 + \left( \frac{p_3}{p_2} - 1 \right) \exp\left( - p_4 \cdot (t - p_1) \right)}$$

- **Richards model**

$$N(t) = \frac{p_1}{[1 + p_2 \cdot e^{-p_3 \cdot (t - p_4)}]^{\frac{1}{p_2}}}$$

- **Weibull**

$$N(t) = p_1 - (p_1 - p_2) \cdot e^{-(p_3 \cdot t)^{p_4}}$$

- **Morgan**

$$N(t) = \frac{p_1 \cdot p_2^{p_3} + p_4 \cdot t^{p_3}}{p_2^{p_3} + t^{p_3}}$$

- **Bertalanffy**

$$N(t) = p_1 + (p_2 - p_1) \cdot (1 - e^{-p_3 \cdot t})^{\frac{1}{p_4}}$$

- **Piece-wise linear-logistic**

$$\begin{cases}
  N(t) = p_2, & t < t_\text{lag} \\
  N(t) = \frac{p_3}{1 + \left( \frac{p_3}{p_2} - 1 \right) \exp\left( - p_4 \cdot (t - p_1) \right)}, & t_\text{lag} \leq t
\end{cases}$$

- **Piece-wise exponential-logistic**

$$\begin{cases}
  N(t) = \text{const}_1 \exp{(p_0 t)}, & t < t_\text{lag} \\
  N(t) = \frac{p_3}{1 + \left( \frac{p_3}{\text{const}_1 \exp^{(p_0 t_\text{lag})}} - 1 \right) \exp\left( - p_4 \cdot (t - p_1) \right)}, & t_\text{lag} \leq t
\end{cases}$$




To call these models use the string present in this table, the parameters will be returned in the same order of this table.


| **Model Name**                   | **Parameters List**                       | **String to call**                     |
| ------------------------------- | ------------------------------------- | ---------------------------------- |
| Exponential                     | $p_1, p_2$                          | `"NL_exponential"`                 |
| Gompertz                        | $p_1, p_2, p_3$                     | `"NL_Gompertz"`                    |
| Logistic                        | $p_1, p_2, p_3, p_4$                | `"NL_logistic"`                    |
| Richards model                  | $p_1, p_2, p_3, p_4$                | `"NL_Richards"`                    |
| Weibull                         | $p_1, p_2, p_3, p_4$                | `"NL_Weibull"`                     |
| Morgan                          | $p_1, p_2, p_3, p_4$                | `"NL_Morgan"`                      |
| Bertalanffy                     | $p_1, p_2, p_3, p_4$                | `"NL_Bertalanffy"`                 |
| piece-wise linear-logistic      | $p_1, p_2, p_3, p_4, t_\text{lag}$  | `"NL_piecewise_lin_logistic"`      |
| piece-wise exponential-logistic | $p_0, p_1, p_3, p_4, t_\text{lag}, \text{const}_1$ | `"NL_piecewise_exp_logistic"` |

## ODEs for bacterial growth

The models implemented in KinBiont are the following:



- **Exponential**:

$$\frac{d N(t)}{dt} =\mu N(t)$$


- **Hyper Gompertz**:

$$\frac{d N(t)}{dt} = \mu_{\text{max}} \cdot \log \left( \frac{N_{\text{max}}}{N(t)} \right)^{(1-n)}$$

- **Hyper Logistic**:

$$\frac{d N(t)}{dt} = \frac{\mu_{\text{max}}}{N_{\text{max}}} \cdot N(t)^{(1-n)} (N(t) - N_{\text{max}})^{(1+n)}$$

- **Von Bertalanffy ODE**:

$$\frac{d N(t)}{dt} = N(t) \cdot (p_1 \cdot p_2 \cdot t^{p_2-1}) - p_3 \cdot p_4 \cdot t^{p_4-1}$$

- **Bertalanffy-Richards**:

$$\frac{d N(t)}{dt} = \frac{t^n}{t^n + \lambda^n} \cdot \mu_{\text{max}} \left( 1 - \frac{N(t)}{N_{\text{max}}} \right) \cdot N(t)$$

- **Logistic**:

$$\frac{d N(t)}{dt} = \mu_{\text{max}} \left( 1 - \frac{N(t)}{N_{\text{max}}} \right) \cdot N(t)$$

- **Adjusted Logistic**:

$$\frac{d N(t)}{dt} = \mu_{\text{max}} \left( 1 - \left(\frac{N(t)}{N_{\text{max}}}\right) ^n \right) \cdot N(t)$$

- **Gompertz**:

$$\frac{d N(t)}{dt} = \mu_{\text{max}} \cdot N(t) \cdot \log \left( \frac{N_{\text{max}}}{N(t)} \right)$$

- **Baranyi Richards**:

$$\frac{d N(t)}{dt} = \frac{t^n}{t^n + \lambda^n} \cdot \mu_{\text{max}} \left( 1 - \frac{N(t)}{N_{\text{max}}} \right) \cdot N(t)$$

- **Baranyi Roberts**:

$$\frac{d N(t)}{dt} = \frac{t^n}{t^n + \lambda^n} \cdot \mu_{\text{max}} \left( 1 - \left( \frac{N(t)}{N_{\text{max}}} \right)^m \right) \cdot N(t)$$

- **Piece-wise Adjusted Logistic**:

$$\frac{d N(t)}{dt} = 
  \begin{cases} 
  \text{const.} \cdot N(t) & t < t_{\text{lag}} \\ 
  \mu_{\text{max}} \left( 1 - \left( \frac{N(t)}{N_{\text{max}}} \right)^m \right) \cdot N(t) & t \geq t_{\text{lag}}
\end{cases}$$

- **Triple Piece-wise Adjusted Logistic**:

$$\frac{d N(t)}{dt} = 
  \begin{cases} 
  \text{const}_1 \cdot N(t) & t < t_{\text{lag}} \\ 
  \mu_{\text{max}} \left( 1 - \left( \frac{N(t)}{N_{\text{max}}} \right)^m \right) \cdot N(t) & t_{\text{lag}} \leq t < t_{\text{stat}} \\ 
  \text{const}_2 \cdot N(t) & t \geq t_{\text{stat}}
\end{cases}$$

- **Triple Piece-wise**:

$$\frac{d N(t)}{dt} = 
  \begin{cases} 
  \text{const}_1 \cdot N(t) & t < t_{\text{lag}} \\ 
  \mu_{\text{max}} \cdot N(t) & t_{\text{lag}} \leq t < t_{\text{stat}} \\ 
  \text{const}_2 \cdot (1 - \log \left( \frac{N(t)}{N_{\text{max}}} \right)) & t \geq t_{\text{stat}}
\end{cases}$$

- **Triple Piece-wise Exponential**:

$$\frac{d N(t)}{dt} = 
  \begin{cases} 
  \text{const}_1 \cdot N(t) & t < t_{\text{lag}} \\ 
  \mu_{\text{max}} \cdot N(t) & t_{\text{lag}} \leq t < t_{\text{stat}} \\ 
  \text{const}_2 \cdot N(t) & t \geq t_{\text{stat}}
\end{cases}$$

- **Four Piece-wise Exponential**:

$$\frac{d N(t)}{dt} = 
  \begin{cases} 
  \text{const}_1 \cdot N(t) & t < t_1 \\ 
  \mu_{\text{max}} \cdot N(t) & t_1 \leq t < t_2 \\ 
  \text{const}_2 \cdot N(t) & t_2 \leq t < t_3 \\ 
  \text{const}_3 \cdot N(t) & t \geq t_3
\end{cases}$$

- **Heterogeneous Population Model (HPM McKellar)**:

$$\begin{cases}
  N(t) = N_1(t) + N_2(t) \\
  \frac{d N_1(t)}{dt} = - \text{rate}_{\text{lag}} \cdot N_1(t) \\
  \frac{d N_2(t)}{dt} = \text{rate}_{\text{lag}} \cdot N_1(t) + \mu_{\text{max}} \cdot N_2(t) \cdot \left(1 - \frac{N_1(t) + N_2(t)}{N_{\text{max}}}\right)
\end{cases}$$

- **Exponential Heterogeneous Population Model**:

$$\begin{cases}
  N(t) = N_1(t) + N_2(t) \\
  \frac{d N_1(t)}{dt} = - \text{rate}_{\text{lag}} \cdot N_1(t) \\
  \frac{d N_2(t)}{dt} = \text{rate}_{\text{lag}} \cdot N_1(t) + \mu_{\text{max}} \cdot N_2(t) 
\end{cases}$$

- **Adjusted Heterogeneous Population Model**:

$$\begin{cases}
  N(t) = N_1(t) + N_2(t) \\
  \frac{d N_1(t)}{dt} = - \text{rate}_{\text{lag}} \cdot N_1(t) \\
  \frac{d N_2(t)}{dt} = \text{rate}_{\text{lag}} \cdot N_1(t) + \mu_{\text{max}} \cdot N_2(t) \cdot \left(1 - \frac{N_1(t) + N_2(t)}{N_{\text{max}}}\right)
\end{cases}$$

- **Heterogeneous Population Model with Inhibition**:

$$\begin{cases}
  N(t) = N_1(t) + N_2(t) + N_3(t), \\
  \frac{d N_1(t)}{dt} = - \text{rate}_{\text{lag}} \cdot N_1(t), \\
  \frac{d N_2(t)}{dt} = \text{rate}_{\text{lag}} \cdot N_1(t) + \mu_{\text{max}} \cdot N_2(t) - \text{rate}_{\text{inhibition}} \cdot N_2(t), \\
  \frac{d N_3(t)}{dt} = \text{rate}_{\text{inhibition}} \cdot N_2(t),
\end{cases}$$


- **Adjusted Heterogeneous Population Model with Inhibition**:

$$\begin{cases}
  N(t) = N_1(t) + N_2(t) + N_3(t), \\
  \frac{d N_1(t)}{dt} = - \text{rate}_{\text{lag}} \cdot N_1(t), \\
  \frac{d N_2(t)}{dt} = \text{rate}_{\text{lag}} \cdot N_1(t) + \mu_{\text{max}} \cdot N_2(t) - \text{rate}_{\text{inhibition}} \cdot N_2(t), \\
  \frac{d N_3(t)}{dt} = \text{rate}_{\text{inhibition}} \cdot N_2(t),
\end{cases}$$


- **Heterogeneous Population Model with Inhibition and Death**:

$$\begin{cases}
  N(t) = N_1(t) + N_2(t) + N_3(t), \\
  \frac{d N_1(t)}{dt} = - \text{rate}_{\text{lag}} \cdot N_1(t), \\
  \frac{d N_2(t)}{dt} = \text{rate}_{\text{lag}} \cdot N_1(t) + \mu_{\text{max}} \cdot N_2(t) - \text{rate}_{\text{inhibition}} \cdot N_2(t), \\
  \frac{d N_3(t)}{dt} = - \text{rate}_{\text{death}} \cdot N_3(t) + \text{rate}_{\text{inhibition}} \cdot N_2(t),
\end{cases}$$

- **Heterogeneous Population Model for Inhibition + Death + Resistance**:

$$\begin{cases}
  N(t) = N_1(t) + N_2(t) + N_3(t), \\
  \frac{d N_1(t)}{dt} = - \text{rate}_{\text{lag}} \cdot N_1(t), \\
  \frac{d N_2(t)}{dt} = \text{rate}_{\text{lag}} \cdot N_1(t) + \mu_{\text{max}} \cdot N_2(t) - \text{rate}_{\text{inhibition}} \cdot N_2(t), \\
  \frac{d N_3(t)}{dt} = - \text{rate}_{\text{death}} \cdot N_3(t) \left(1 - \frac{N_3(t)}{N_{\text{res}}}\right) + \text{rate}_{\text{inhibition}} \cdot N_2(t),
\end{cases}$$

To call these models use the string present in this table, the parameters will be returned in the same order of this table.



| **Model Name**                               | **Parameters**                                                                                       | **String to call**                             |
|----------------------------------------------|------------------------------------------------------------------------------------------------------|------------------------------------------|
| Exponential ODE                              | `label_exp`, `well`, `model`, `gr`, `th_max_gr`, `emp_max_gr`, `loss`           | `"exponential"`                         |
| Hyper Gompertz                              | `label_exp`, `well`, `model`, `gr`, `N_max`, `shape`, `th_max_gr`, `emp_max_gr`, `loss`           | `"hyper_gompertz"`                         |
| Hyper Logistic                              | `label_exp`, `well`, `model`, `doubling_time`, `gr`, `N_max`, `shape`, `th_max_gr`, `emp_max_gr`, `loss` | `"hyper_logistic"`                         |
| Von Bertalanffy ODE                         | `label_exp`, `well`, `model`, `alpha`, `beta`, `a`, `b`, `th_max_gr`, `emp_max_gr`, `loss`         | `"ode_von_bertalanffy"`                    |
| Bertalanffy-Richards                        | `label_exp`, `well`, `model`, `gr`, `N_max`, `shape`, `th_max_gr`, `emp_max_gr`, `loss`           | `"bertalanffy_richards"`                  |
| Logistic                                    | `label_exp`, `well`, `model`, `gr`, `N_max`, `th_max_gr`, `emp_max_gr`, `loss`                    | `"logistic"`                               |
| Adjusted Logistic                           | `label_exp`, `well`, `model`, `gr`, `N_max`, `shape`, `th_max_gr`, `emp_max_gr`, `loss`           | `"alogistic"`                              |
| Gompertz                                    | `label_exp`, `well`, `model`, `gr`, `N_max`, `th_max_gr`, `emp_max_gr`, `loss`                    | `"gompertz"`                              |
| Baranyi Richards                            | `label_exp`, `well`, `model`, `gr`, `N_max`, `lag_time`, `shape`, `th_max_gr`, `emp_max_gr`, `loss`| `"baranyi_richards"`                       |
| Baranyi Roberts                             | `label_exp`, `well`, `model`, `gr`, `N_max`, `lag_time`, `shape_1`, `shape_2`, `th_max_gr`, `emp_max_gr`, `loss` | `"baranyi_roberts"`                        |
| Piece-wise Adjusted Logistic                | `label_exp`, `well`, `model`, `gr`, `N_max`, `lag`, `shape`, `linear_const`, `th_max_gr`, `emp_max_gr`, `loss` | `"piecewise_adjusted_logistic"`        |
| Triple Piece-wise Adjusted Logistic         | `label_exp`, `well`, `model`, `gr`, `N_max`, `lag`, `shape`, `linear_const`, `t_stationary`, `linear_lag`, `th_max_gr`, `emp_max_gr`, `loss` | `"triple_piecewise_adjusted_logistic"` |
| Triple Piece-wise                           | `label_exp`, `well`, `model`, `gr`, `gr_2`, `gr_3`, `lag`, `t_stationary`, `th_max_gr`, `emp_max_gr`, `loss` | `"ODE_triple_piecewise"`                   |
| Triple Piece-wise Exponential               | `label_exp`, `well`, `model`, `gr`, `gr_2`, `gr_3`, `lag`, `t_stationary`, `th_max_gr`, `emp_max_gr`, `loss` | `"ODE_triple_piecewise_exponential"`       |
| Four Piece-wise Exponential                 | `label_exp`, `well`, `model`, `gr`, `gr_2`, `gr_3`, `gr_4`, `lag`, `t_decay_gr`, `t_stationary`, `th_max_gr`, `emp_max_gr`, `loss` | `"ODE_four_piecewise"`                     |
| Diauxic Piecewise Adjusted Logistic         | `label_exp`, `well`, `model`, `gr_1`, `N_max`, `shape_1`, `lag`, `linear_const`, `t_shift`, `gr_2`, `N_max_2`, `shape_2`, `end_second_lag`, `lag_2_gr`, `th_max_gr`, `emp_max_gr`, `loss` | `"Diauxic_piecewise_adjusted_logistic"`                     |
| Heterogeneous Population Model (HPM McKellar) | `label_exp`, `well`, `model`, `gr`, `exit_lag_rate`, `N_max`, `th_max_gr`, `emp_max_gr`, `loss`   | `"HPM"`                          |
| Exponential Heterogeneous Population Model (HPM McKellar) | `label_exp`, `well`, `model`, `gr`, `exit_lag_rate`, `th_max_gr`, `emp_max_gr`, `loss`   | `"HPM_exp"`                          |
| Adjusted Heterogeneous Population Model     | `label_exp`, `well`, `model`, `gr`, `exit_lag_rate`, `N_max`, `shape`, `th_max_gr`, `emp_max_gr`, `loss` | `"aHPM"`                 |
| Heterogeneous Population Model with Inhibition | `label_exp`, `well`, `model`, `gr`, `exit_lag_rate`, `inactivation_rate`, `th_max_gr`, `emp_max_gr`, `loss`| `"HPM_3_inhibition"`                    |
| Heterogeneous Population Model with Inhibition and Death | `label_exp`, `well`, `model`, `gr`, `exit_lag_rate`, `inactivation_rate`, `death_rate`, `th_max_gr`, `emp_max_gr`, `loss` | `"HPM_3_death"`                         |
| Heterogeneous Population Model for Inhibition + Death + Resistance | `label_exp`, `well`, `model`, `gr`, `exit_lag_rate`, `inactivation_rate`, `death_rate`, `n_res`, `shape`, `th_max_gr`, `emp_max_gr`, `loss` | `"aHPM_3_death_resistance"`             |

## Stochastic models for bacterial growth


In the stochastic version of the growth models, the growth rate of each population component (denoted as $\mu_i$) is evaluated based on the concentration of the limiting nutrient and then the number of birth event is evaluated with the Poisson approximation. The user is required to specify the starting amount of nutrients and the volume of the solution. Various kinetic growth models are considered. Note that these models can be used only during simulations
In the following, $[\text{Nut.}]$ is the limiting nutrient concentration, $\mu_\text{max}$ is the maximum possible growth rate, $k_1$ and $k_2$ are numerical constants with meanings depending on the selected model.

- Monod Model




$\mu([\text{Nut.}]; k_1, \mu_\text{max}) = \mu_\text{max} \frac{[\text{Nut.}]}{k_1 + [\text{Nut.}]}.$

- Haldane Model




$\mu([\text{Nut.}]; k_1, k_2, \mu_\text{max}) = \mu_\text{max} \frac{[\text{Nut.}]}{k_1 + [\text{Nut.}] + \frac{k_2}{[\text{Nut.}]^2}}.$

- Blackman Model




$\mu([\text{Nut.}]; k_1, \mu_\text{max}) = \mu_\text{max} \frac{[\text{Nut.}]}{k_1 + [\text{Nut.}]}.$

- Tessier Model




$\mu([\text{Nut.}]; k_1, \mu_\text{max}) = \mu_\text{max} (1 - e^{k_1[\text{Nut.}] }).$

- Moser Model




$\mu([\text{Nut.}]; k_1, k_2, \mu_\text{max}) = \mu_\text{max} \frac{[\text{Nut.}]^{k_2}}{k_1 + [\text{Nut.}]^{k_2}}.$

- Aiba-Edwards Model



$\mu([\text{Nut.}]; k_1, k_2, \mu_\text{max}) = \mu_\text{max} \frac{[\text{Nut.}]}{k_1 + [\text{Nut.}]} e^{-\frac{[\text{Nut.}]}{k_2}}.$

- Verhulst Model




$\mu(N; N_\text{max}, \mu_\text{max}) = \mu_\text{max} \left(1 - \frac{N}{N_\text{max}}\right).$




## Error functions

`type_of_loss = "L2"`: Minimize the L2 norm of the difference between the numerical solution of the desired model and the given data.

$$\mathcal{D}(D(t_i), \bar{N}(t_i, \{P\})) = \frac{1}{n} \sum_{i=1}^n \left(D(t_i) - \bar{N}(t_i, \{P\})\right)^2$$

where $n$ is the number of data points.

`type_of_loss = "RE"`: Minimize the relative error between the solution and data.

$$\mathcal{D}(D(t_i), \bar{N}(t_i, \{P\})) = \frac{1}{n} \sum_{i=1}^n 0.5 \cdot \left(1 - \frac{D(t_i)}{\bar{N}(t_i, \{P\})}\right)^2$$

where $n$ is the number of data points.

`type_of_loss = "L2_derivative"`: Minimize the L2 norm of the difference between the specific growth rate of the numerical solution of the desired model and the corresponding derivatives of the data.

$$\mathcal{D}\left(\frac{dD(t_i)}{dt}, \frac{d\bar{N}(t_i, \{P\})}{dt}\right) = \frac{1}{n} \sum_{i=1}^n \left(\frac{dD(t_i)}{dt} - \frac{d\bar{N}(t_i, \{P\})}{dt}\right)^2$$

where $n$ is the number of data points.

`type_of_loss = "blank_weighted_L2"`: Minimize a weighted version of the L2 norm, where the difference between the solution and data is weighted based on a distribution obtained from empirical blank data.

$$\mathcal{D}(D(t_i), \bar{N}(t_i, \{P\})) = \frac{1}{n} \sum_{i=1}^n \left(1 - P(D(t_i) - \bar{N}(t_i, \{P\})|\text{noise})\right) \cdot \left(D(t_i) - \bar{N}(t_i, \{P\})\right)^2$$

where $P(D(t_i) - \bar{N}(t_i, \{P\})|\text{noise})$ is the probability distribution of the empirical blank data, and $n$ is the number of data points.

`type_of_loss = "L2_log"`: Minimize the logarithm of the L2 norm of the difference between the numerical solution of the desired model and the given data.

$$\mathcal{D}(D(t_i), \bar{N}(t_i, \{P\})) = \log\left(\frac{1}{n} \sum_{i=1}^n \left(D(t_i) - \bar{N}(t_i, \{P\})\right)^2\right)$$

where $n$ is the number of data points.

`type_of_loss = "RE_log"`: Minimize the logarithm of the relative error between the solution and data.

$$\mathcal{D}(D(t_i), \bar{N}(t_i, \{P\})) = \log\left(\frac{1}{n} \sum_{i=1}^n 0.5 \cdot \left(1 - \frac{D(t_i)}{\bar{N}(t_i, \{P\})}\right)^2\right)$$

where $n$ is the number of data points.

`type_of_loss = "L2_std_blank"`: Minimize the L2 norm of the difference between the numerical solution of the desired model and the data, normalized by the standard deviation of empirical blank data.

$$\mathcal{D}(D(t_i), \bar{N}(t_i, \{P\})) = \frac{1}{n} \sum_{i=1}^n \left(\frac{D(t_i) - \bar{N}(t_i, \{P\})}{\text{std\_blank}}\right)^2$$

where $\text{std\_blank}$ is the standard deviation of the empirical blank data, and $n$ is the number of data points.

---