# [The mathematical models](@id models)

In Kinbiont it is possible to simulate and fit the bacterial growth with any ordinary differential equations system.  We can broadly divide the classes of possible mathematical models into the following:

```@contents
Pages = ["index.md"]
Depth = 2
```

In this section we show some of the harcoded models of Kinbiont.jl but note that you can input any custom model both has analytic functio both as ODE.



## NL models for bacterial growth

In general, NL fitting is preferred to ODE fitting in the following cases:
1.  Since it is faster when you have to analyze a large dataset
2. When you do not trust initial conditions (e.g., initial inoculum under the detection limit of the instrument). ODE fit needs to fix the initial condition of data on the first (or an average of the first) time point of data, and this could lead to errors.

In this case, we are supposed to know the analytic formula of microbial growth; in particular, we have implemented some models from "Statistical evaluation of mathematical models for microbial growth" and added some piecewise models. They are:

- **Exponential**  

  $$N(t) = N_0 \cdot e^{\mu \cdot t}$$

  where $\mu$ is the growth rate, and $N_0$ is the starting condition.

- **Gompertz**  

  $$N(t) = N_{\text{max}} \cdot e^{-e^{-\mu \cdot (t - t_{\text{L}})}}$$

  where $\mu$ is the growth rate, $N_{\text{max}}$ is the total growth, and $t_{\text{L}}$ is the lag time.

- **Logistic**  

  $$N(t) = \frac{N_{\text{max}}}{1 + \left( \frac{N_{\text{max}}}{N_0} - 1 \right) \exp\left( - \mu \cdot t \right)}$$

  where $\mu$ is the growth rate, $N_0$ is the starting condition, and $N_{\text{max}}$ is the total growth.

- **Richards model**  

  $$N(t) = \frac{N_{\text{max}}}{[1 + \nu \cdot e^{-\mu \cdot (t - t^{\text{L}})}]^{\frac{1}{\nu}}}$$

  where $\mu$ is the growth rate, $N_{\text{max}}$ is the total growth, $t_{\text{L}}$ is the lag time, and $\nu$ is a shape constant.

- **Weibull**  

  $$N(t) = N_{\text{max}} - (N_{\text{max}} - N_0) \cdot e^{-(\mu \cdot t)^{\nu}}$$

  where $\mu$ is the growth rate, $N_0$ is the starting condition, $N_{\text{max}}$ is the total growth, and $\nu$ is a shape constant.

- **Morgan**  

  $$N(t) = \frac{N_0 \cdot K^{\nu} + N_{\text{max}} \cdot t^{\nu}}{K^{\nu} + t^{\nu}}$$

  where $N_0$ is the starting condition, $N_{\text{max}}$ is the total growth, and $K$ and $\nu$ are shape constants.

- **Bertalanffy**  

  $$N(t) = N_0 + (N_{\text{max}} - N_0) \cdot (1 - e^{-\mu \cdot t})^{\frac{1}{\nu}}$$

  where $N_0$ is the starting condition, $N_{\text{max}}$ is the total growth, $\mu$ is the growth rate, and $\nu$ is a shape constant.

- **Piece-Wise Linear-Logistic**  

  $$N(t) = 
  \begin{cases} 
  N_0, & t < t_{\text{L}} \\ 
  N(t) = \frac{N_{\text{max}}}{1 + \left( \frac{N_{\text{max}}}{N_0} - 1 \right) \exp\left( - \mu \cdot (t - t_{\text{L}}) \right)}, & t_{\text{L}} \leq t 
  \end{cases}$$

  where $N_0$ is the starting condition, $N_{\text{max}}$ is the total growth, $\mu$ is the growth rate, and $t_{\text{L}}$ is the lag time.

- **Piece-wise Exponential-Logistic** 

  $$N(t) = 
  \begin{cases} 
  N_0 \exp{(\mu_0 \cdot t)}, & t < t_{\text{L}} \\ 
  \frac{N_{\text{max}}}{1 + \left( \frac{N_{\text{max}}}{N_0 \exp{(\mu_0 \cdot t_{\text{L}})}} - 1 \right) \exp\left( - \mu \cdot (t - t_{\text{L}}) \right)}, & t_{\text{L}} \leq t 
  \end{cases}$$

  where $N_0$ is the starting condition, $N_{\text{max}}$ is the total growth, $\mu$ is the growth rate, $t_{\text{L}}$ is the lag time, and $\mu_0$ is the growth rate during the lag phase.




To call these models use the string present in this table, the parameters will be returned in the same order of this table.


| **Model Name**                   | **Parameters List**                       | **String to call**                     |
| ------------------------------- | ------------------------------------- | ---------------------------------- |
| Exponential                     | $N_0, \mu$                          | `"NL_exponential"`                 |
| Gompertz                        | $N_{\text{max}}, \mu, t_{\text{L}}$                     | `"NL_Gompertz"`                    |
| Logistic                        | $N_{\text{max}}, N_0,\mu$                | `"NL_logistic"`                    |
| Richards model                  | $N_{\text{max}}, \mu,\nu,t_{\text{L}}$                | `"NL_Richards"`                    |
| Weibull                         | $N_{\text{max}}, N_0,\mu,\nu$                | `"NL_Weibull"`                     |
| Morgan                          | $N_{\text{max}}, N_0,K,\nu$                | `"NL_Morgan"`                      |
| Bertalanffy                     | $N_{\text{max}}, N_0,\mu,\nu$                | `"NL_Bertalanffy"`                 |
| piece-wise linear-logistic      | $N_0, N_{\text{max}},\mu, t_\text{L}$  | `"NL_piecewise_lin_logistic"`      |
| piece-wise exponential-logistic | $N_0, N_{\text{max}},\mu, t_\text{L}, t_\text{L},\mu_0$ | `"NL_piecewise_exp_logistic"` |





For a general idea of the properties of models, consult the following table:


| **Model Name**                   | **Has Lag?** | **Is Piecewise?** | **Has Stationary Phase?** |
| --------------------------------- | ---------- | -------------- | ------------------- |
| Exponential                       | No         | No             | No                    |
| Gompertz                          | Yes        | No             | Yes                   |
| Logistic                          | No         | No             | Yes                   |
| Richards model                    | Yes        | No             | Yes                   |
| Weibull                           | No         | No             | Yes                   |
| Morgan                            | No         | No             | Yes                   |
| Bertalanffy                       | No         | No             | Yes                   |
| Piece-wise linear-logistic        | Yes        | Yes            | Yes                   |
| Piece-wise exponential-logistic   | Yes        | Yes            | Yes                   |

If undecided between different models, please use the model selection function.


## ODEs for bacterial growth

The models implemented in Kinbiont are the following:



- **Exponential**:

$$\frac{d N(t)}{dt} =\mu N(t)$$
where $\mu$ is the growth rate.


- **Hyper Gompertz**:

$$\frac{d N(t)}{dt} = \mu \, \log \left( \frac{N_{\text{max}}}{N(t)} \right)^{(1-n)}$$
where $\mu$ is the growth rate, $N_{\text{max}}$ the total growth and $n$ a shape constant.
- **Hyper Logistic**:

$$\frac{d N(t)}{dt} = \frac{\mu}{N_{\text{max}}} \, N(t)^{(1-n)} (N(t) - N_{\text{max}})^{(1+n)}$$
where $\mu$ is the growth rate, $N_{\text{max}}$ the total growth and $n$ a shape constant.

- **Bertalanffy-Richards**:
$$\frac{d N(t)}{dt} =  \frac{\mu}{N_\text{max}^n } \cdot  \left ( N_\text{max}^n - N^n(t) \right)$$

where $\mu$ is the growth rate, $N_{\text{max}}$ the total growth, and $n$ a shape constant.

- **Logistic**:

$$\frac{d N(t)}{dt} = \mu \left( 1 - \frac{N(t)}{N_{\text{max}}} \right) \, N(t)$$
where $\mu$ is the growth rate, and $N_{\text{max}}$ the total growth.

- **Adjusted Logistic**:

$$\frac{d N(t)}{dt} = \mu \left( 1 - \left(\frac{N(t)}{N_{\text{max}}}\right) ^n \right) \, N(t)$$
where $\mu$ is the growth rate, $N_{\text{max}}$ the total growth and $n$ a shape constant.
- **Gompertz**:

$$\frac{d N(t)}{dt} = \mu \, N(t) \, \log \left( \frac{N_{\text{max}}}{N(t)} \right)$$
where $\mu$ is the growth rate, and $N_{\text{max}}$ the total growth.

- **Baranyi Richards**:

$$\frac{d N(t)}{dt} = \frac{t^n}{t^n + \lambda^n} \, \mu \left( 1 - \frac{N(t)}{N_{\text{max}}} \right) \, N(t)$$
where $\mu$ is the growth rate, $N_{\text{max}}$ the total growth, $\lambda$ is the lag time and $n$ a shape constant.

- **Baranyi Roberts**:

$$\frac{d N(t)}{dt} = \frac{t^n}{t^n + \lambda^n} \, \mu \left( 1 - \left( \frac{N(t)}{N_{\text{max}}} \right)^m \right) \, N(t)$$
where $\mu$ is the growth rate, $N_{\text{max}}$ the total growth, $\lambda$ is the lag time,  $n$ and $m$  are shape constants.

- **Piece-wise Adjusted Logistic**:

$$\frac{d N(t)}{dt} = 
  \begin{cases} 
  \text{const.} \, N(t) & t < t_{\text{L}} \\ 
  \mu \left( 1 - \left( \frac{N(t)}{N_{\text{max}}} \right)^m \right) \, N(t) & t \geq t_{\text{L}}
\end{cases}$$
where $\mu$ is the growth rate, $N_{\text{max}}$ the total growth, $t_\text{L}$ is the lag time,    $m$  is shape constant, and $c$ the growth  rate during the lag phase (can be 0).
- **Triple Piece-wise Adjusted Logistic**:

$$\frac{d N(t)}{dt} = 
  \begin{cases} 
  c_1 \cdot N(t) & \text{for } t < t_{\text{L}}, \\ 
  \mu \left( 1 - \left( \frac{N(t)}{N_{\text{max}}} \right)^m \right) \cdot N(t) & \text{for } t_{\text{L}} \leq t < t_{\text{stat}}, \\ 
  c_2 \cdot N(t) & \text{for } t \geq t_{\text{stat}},
\end{cases}$$

 where $\mu$ is the growth rate, $N_{\text{max}}$ the total growth, $t_\text{L}$ is the lag time,    $m$  is a shape constant,  $c_1$ the growth rate during the lag phase (can be 0), $t_{\text{stat}} $ the time when stationary phase starts, and $c_2$ the growth rate during the stationary phase.
- **Triple Piece-wise**:

$$\frac{d N(t)}{dt} = 
  \begin{cases} 
  c_1 \cdot N(t) & \text{for } t < t_{\text{L}}, \\ 
  \mu \cdot N(t) & \text{for } t_{\text{L}} \leq t < t_{\text{stat}},\\ 
  c_2 \cdot \left(1 - \log \left( \frac{N(t)}{N_{\text{max}}} \right)\right) & \text{for } t \geq t_{\text{stat}},
\end{cases}$$

where $\mu$ is the growth rate, $N_{\text{max}}$ the total growth, $t_\text{L}$ is the lag time,       $c_1$ the growth rate during the lag phase (can be 0), $t_{\text{stat}} $ the time when stationary phase starts, and $c_2$ the growth rate during the stationary phase.

- **Triple Piece-wise Exponential**:

$$\frac{d N(t)}{dt} = 
  \begin{cases} 
  c_1 \cdot N(t) & \text{for } t < t_{\text{L}}, \\ 
  \mu \cdot N(t) & \text{for } t_{\text{L}} \leq t < t_{\text{stat}}, \\ 
  c_2 \cdot N(t) & \text{for } t \geq t_{\text{stat}},
\end{cases}$$

where $\mu$ is the growth rate, $N_{\text{max}}$ the total growth, $t_\text{L}$ is the lag time,    $c_1$ the growth  rate during the lag phase (can be 0), $t_{\text{stat}} $ the time when stationary phase starts, and $c_2$ the growth rate during the stationary phase.
- **Four Piece-wise Exponential**:

$$\frac{d N(t)}{dt} = 
  \begin{cases} 
  c_1 \cdot N(t) & \text{for } t < t_1, \\ 
  \mu \cdot N(t) & \text{for } t_1 \leq t < t_2, \\ 
  c_2 \cdot N(t) & \text{for } t_2 \leq t < t_3,\\ 
  c_3 \cdot N(t) & \text{for } t \geq t_3,
\end{cases}$$

where $\mu$ is the growth rate, $N_{\text{max}}$ the total growth, $t_1$ is the lag time,    $c_1$ the growth rate during the lag phase (can be 0), $t_2 $ the time when a growth after exponential growths,  $c_2$ the growth rate during this phase, $t_3$ the start of stationary phase and, $c_3$ the growth  rate during the stationary phase.

- **Heterogeneous Population Model (HPM)**:
$$\begin{cases}
  N(t) = N_1(t) + N_2(t), \\
  \frac{d N_1(t)}{dt} = - r_{\text{L}} \cdot N_1(t), \\
  \frac{d N_2(t)}{dt} = r_{\text{L}} \cdot N_1(t) + \mu \cdot N_2(t) \cdot \left(1 - \frac{N_1(t) + N_2(t)}{N_{\text{max}}}\right),
\end{cases}$$

where $\mu$ is the growth rate, $N_{\text{max}}$ the total growth, and $r_\text{L}$ is the lag rate (i.e. the rate of transition between $N_1(t)$ and $N_2(t)$).     
Note that these models assume that the cells are in two states: $N_1(t)$ dormant cells (the cells are not able to reproduce because they are in the lag phase) and $N_2(t)$ active cells, which are able to duplicate. At the start, all the cells are assumed in the dormant state (i.e., $N_{1}(start) = OD(start)$, and $N_{2}(start) = 0.0$) .

- **Exponential Heterogeneous Population Model**:

$$\begin{cases}
  N(t) = N_1(t) + N_2(t) \\
  \frac{d N_1(t)}{dt} = - \text{r}_{\text{L}} \, N_1(t) \\
  \frac{d N_2(t)}{dt} = \text{r}_{\text{L}} \, N_1(t) + \mu \, N_2(t) 
\end{cases}$$

where similarly to the HPM model, $N_1$ and $N_2$ refer to the populations of dormant and active cells, respectively. $\mu$ is the growth rate, and the lag rate $r_\text{L}$ denotes the transition between the $N_1$ and $N_2$ populations.     Here, we also assume that all cells are in the dormant state at the start (i.e., $N_{1}(t = 0) = \text{OD}(t = 0)$, and $N_{2}(t = 0) = 0$).


- **Adjusted Heterogeneous Population Model**:

$$\begin{cases}
  N(t) = N_1(t) + N_2(t), \\
  \frac{d N_1(t)}{dt} = - r_{\text{L}} \cdot N_1(t), \\
  \frac{d N_2(t)}{dt} = r_{\text{L}} \cdot N_1(t) + \mu \cdot N_2(t) \cdot \left(1 - \left(\frac{N_1(t) + N_2(t)}{N_{\text{max}}}\right)^m \right),
\end{cases}$$

where $\mu$ is the growth rate, and $N_{\text{max}}$ the total growth.
Note that these models assume that the cells are in two states: $N_1(t)$ dormant cells (the cells are not able to reproduce because they are in the lag phase) and $N_2(t)$ active cells, which are able to duplicate.At the start, all the cells are assumed in the dormant state (i.e., $N_{1}(start) = OD(start)$, and $N_{2}(start) = 0.0$) .

- **Heterogeneous Population Model with Inhibition**:

$$\begin{cases}
  N(t) = N_1(t) + N_2(t)+ N_3(t)  \\
  \frac{d N_1(t)}{dt} = - r_{\text{L}} \cdot N_1(t) \\
  \frac{d N_2(t)}{dt} = r_{\text{L}} \cdot N_1(t) - r_{\text{I}} \cdot  N_2(t)\\
  \frac{d N_3(t)}{dt} =  r_{\text{I}} \cdot  N_2(t)
\end{cases}$$


where $\mu$ is the growth rate, $N_{\text{max}}$ the total growth,  $r_\text{L}$ is the lag rate (i.e. the rate of transition between $N_1(t)$ and $N_2(t)$) and $r_{\text{I}}$ a shape constant.     
Note that these models assume that the cells are in two states: $N_1(t)$ dormant cells (the cells are not able to reproduce because they are in the lag phase), $N_2(t)$ active cells, which are able to duplicate, and  inactive cells $N_3(t)$. At the start, all the cells are assumed in the dormant state (i.e., $N_{1}(\text{start}) = OD(\text{start})$, $N_{2}(\text{start}) = 0.0$, and $N_{3}(\text{start}) = 0.0$).

- **Heterogeneous Population Model with Inhibition and Death**:

$$\begin{cases} N(t) = N_1(t) + N_2(t) + N_3(t), \\
\frac{d N_1(t)}{dt} = - r_{\text{L}} \cdot N_1(t), \\
\frac{d N_2(t)}{dt} = r_{\text{L}} \cdot N_1(t) + \mu \cdot N_2(t) - r_{\text{I}} \cdot N_2(t), \\
\frac{d N_3(t)}{dt} = - r_{\text{D}} \cdot N_3(t) + r_{\text{I}} \cdot N_2(t), 
\end{cases}$$

where $\mu$ is the growth rate, $r_\text{L}$ is the lag rate (i.e. the rate of transition between $N_1(t)$ and $N_2(t)$) ,  $r_\text{inhibition}$ is the  rate of which cell are inhibited (i.e. the rate of transition between $N_2(t)$ and $N_3(t)$), and $r_{\text{D}}$ is the  rate of which cell are die.


Note that these models assume that the cells are in three states: $N_1(t)$ dormant cells (the cells are not able to reproduce because they are in the lag phase), $N_2(t)$ active cells, which are able to duplicate, and inactive cells $N_3(t)$ that die at a rate $r_{\text{D}}$. At the start, all the cells are assumed in the dormant state (i.e., $N_{1}(\text{start}) = OD(\text{start})$, $N_{2}(\text{start}) = 0.0$, and $N_{3}(\text{start}) = 0.0$).

- **Heterogeneous Population Model with Inhibition, Death and Resistance**:

$$\begin{cases}
N(t) = N_1(t) + N_2(t) + N_3(t), \\
\frac{d N_1(t)}{dt} = - r_{\text{L}} \cdot N_1(t), \\
\frac{d N_2(t)}{dt} = r_{\text{L}} \cdot N_1(t) + \mu \cdot N_2(t) - r_{\text{I}} \cdot N_2(t), \\
\frac{d N_3(t)}{dt} = - r_{\text{D}} \cdot N_3(t) \left(1 - \frac{N_3(t)}{N_{\text{res}}}\right) + r_{\text{I}} \cdot N_2(t), \quad \text{with} \quad N_3(t) \leq N_{\text{res}}
\end{cases}$$

where $\mu$ is the growth rate, $r_\text{L}$ is the lag rate (i.e. the rate of transition between $N_1(t)$ and $N_2(t)$) ,  $r_\text{inhibition}$ is the  rate of which cell are inhibited (i.e. the rate of transition between $N_2(t)$ and $N_3(t)$),  $r_{\text{D}}$ is the  rate of which cell are die, and $N_{\text{res}}$ it the number of cells that will be inactive but do not die.



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
| Heterogeneous Population Model with Inhibition, Death and Resistance | `label_exp`, `well`, `model`, `gr`, `exit_lag_rate`, `inactivation_rate`, `death_rate`, `n_res`, `shape`, `th_max_gr`, `emp_max_gr`, `loss` | `"aHPM_3_death_resistance"`             |


In the following table, you can find a general description of the properties of the hardcoded ODE models of Kinbiont:


| **Model Name**                               | **Has Lag?** | **Is Piecewise?** | **Has Stationary Phase?** | **Has Inhibition?** | **Is Monotonic?** | **Supposes Multiple States?** |
|----------------------------------------------|-------------|------------------|--------------------------|-------------------|----------------|----------------------------|
| Exponential ODE                              | No          | No               | No                       | No                | Yes            | No                         |
| Hyper Gompertz                               | No          | No               | Yes                      | No                | Yes            | No                         |
| Hyper Logistic                               | No          | No               | Yes                      | No                | Yes            | No                         |
| Von Bertalanffy ODE                          | No          | No               | Yes                      | No                | Yes            | No                         |
| Bertalanffy-Richards                         | No          | No               | Yes                      | No                | Yes            | No                         |
| Logistic                                     | No          | No               | Yes                      | No                | Yes            | No                         |
| Adjusted Logistic                            | No          | No               | Yes                      | No                | Yes            | No                         |
| Gompertz                                     | No          | No               | Yes                      | No                | Yes            | No                         |
| Baranyi Richards                             | Yes         | No               | Yes                      | No                | Yes            | No                         |
| Baranyi Roberts                              | Yes         | No               | Yes                      | No                | Yes            | No                         |
| Piece-wise Adjusted Logistic                 | Yes         | Yes              | Yes                      | No                | No             | No                         |
| Triple Piece-wise Adjusted Logistic          | Yes         | Yes              | Yes                      | No                | No             | No                         |
| Triple Piece-wise                            | Yes         | Yes              | Yes                      | No                | No             | No                         |
| Triple Piece-wise Exponential                | Yes         | Yes              | Yes                      | No                | No             | No                         |
| Four Piece-wise Exponential                  | Yes         | Yes              | Yes                      | No                | No             | No                         |
| Diauxic Piecewise Adjusted Logistic          | Yes         | Yes              | Yes                      | No                | No             | No                         |
| Heterogeneous Population Model (HPM McKellar) | Yes        | No               | Yes                      | No                | Yes            | Yes                        |
| Exponential Heterogeneous Population Model (HPM McKellar) | Yes | No | No | No | Yes | Yes |
| Adjusted Heterogeneous Population Model      | Yes         | No               | Yes                      | No                | Yes            | Yes                        |
| Heterogeneous Population Model with Inhibition | Yes       | No               | Yes                      | Yes               | Yes             | Yes                        |
| Heterogeneous Population Model with Inhibition and Death | Yes | No | Yes | Yes | No | Yes |
| Heterogeneous Population Model with Inhibition, Death and Resistance | Yes | No | Yes | Yes | No | Yes |

If undecided between different models, the model selection function should be used.

## Stochastic models for bacterial growth


In the stochastic version of the growth models, the growth rate of each population component (denoted as $\mu_i$) is evaluated based on the concentration of the limiting nutrient, and then the number of birth events is evaluated with the Poisson approximation. The user is required to specify the starting amount of nutrients and the volume of the solution. Various kinetic growth models are considered. **Note that these models can be used only during simulations.**
In the following, we use $\nu$ to represent the limiting nutrient concentration throughout, $\mu_\text{max}$ denotes the maximum possible growth rate, $k_1$ (for $i=1,2$) is a numerical constant whose specific meaning depends on the model, $N$ indicates the number of present cells, and $N_\text{max}$ is the carrying capacity in the Verhulst model.

- Monod Model:


$\mu(\nu;k_1,\mu_\text{max})  = \displaystyle{\mu_\text{max} \frac{\nu}{k_1+\nu}}.$

- Haldane Model:



$\mu(\nu;k_1,k_2,\mu_\text{max}) = \displaystyle{\mu_\text{max} \frac{\nu}{k_1+\nu +\frac{\nu^2}{k_2}}}.$

- Blackman Model:


$\mu(\nu;k_1,\mu_\text{max}) = \displaystyle{\mu_\text{max} \frac{\nu}{k_1}}.$

- Tessier Model:



$\mu(\nu;k_1,\mu_\text{max}) = \displaystyle{\mu_\text{max} (1 - e^{k_1\nu })}.$

- Moser Model:



$\mu(\nu;k_1,k_2,\mu_\text{max})  = \displaystyle{\mu_\text{max} \frac{\nu^{k_2}}{k_1+\nu^{k_2}}}.$

- Aiba-Edwards Model:

$\mu([\text{Nut.}]; k_1, k_2, \mu_\text{max}) = \mu_\text{max} \frac{[\text{Nut.}]}{k_1 + [\text{Nut.}]} e^{-\frac{[\text{Nut.}]}{k_2}}.$

- Verhulst Model:


$\mu(N;N_\text{max},\mu_\text{max})  = \displaystyle{\mu_\text{max} \left(1-\frac{N}{N_\text{max}}\right)}.$



## ODEs system for bacterial growth


In this section, we present some examples of ODEs multidimensional system hardcoded in Kinbiont. Note that these are just examples since you can define custom models:

- **SIR Model** (Susceptible-Infected-Recovered)  
$$\begin{cases}
    \frac{dS}{dt} = -\beta S I \\
    \frac{dI}{dt} = \beta S I - \gamma I \\
    \frac{dR}{dt} = \gamma I
  \end{cases}$$
  Parameters: Infection rate ($\beta$), Recovery rate ($\gamma$).

- **SIR with Birth and Death (SIR_BD)**  
$$\begin{cases}
    \frac{dS}{dt} = -\beta S I + b S - d S \\
    \frac{dI}{dt} = \beta S I - \gamma I - d I \\
    \frac{dR}{dt} = \gamma I - d R
  \end{cases}$$
  Parameters: Infection rate ($\beta$), Recovery rate ($\gamma$), Birth rate ($b$), Death rate ($d$).

- **SIS Model** (Susceptible-Infected-Susceptible)  
$$\begin{cases}
    \frac{dS}{dt} = -\beta S I + \gamma I \\
    \frac{dI}{dt} = \beta S I - \gamma I
  \end{cases}$$
  Parameters: Infection rate ($\beta$), Recovery rate ($\gamma$).

- **Lotka-Volterra Predator-Prey Model**  
$$\begin{cases}
    \frac{dP}{dt} = \alpha P - \beta P C \\
    \frac{dC}{dt} = -\delta C + \gamma C P
  \end{cases}$$
  Parameters: Prey birth rate ($\alpha$), Predation rate ($\beta$), Predator death rate ($\delta$), Predator reproduction rate ($\gamma$).

- **Lotka-Volterra with Substrate Limitation**  
$$\begin{cases}
    \frac{dP}{dt} = \alpha P \frac{S}{S + K} - \beta P C \\
    \frac{dC}{dt} = -\delta C + \gamma C P \\
    \frac{dS}{dt} = -\alpha P \frac{S}{S + K}
  \end{cases}$$
  Parameters: Growth rate ($\alpha$), Half-saturation ($K$), Predation rate ($\beta$), Predator mortality ($\delta$), Predator efficiency ($\gamma$).

- **Monod Chemostat Model** (Microbial Growth in a Chemostat)  
$$\begin{cases}
    \frac{dX}{dt} = \mu X - D X \\
    \frac{dS}{dt} = D (S_{\text{in}} - S) - \frac{\mu X}{Y} - m X
  \end{cases}$$
  where  
$$\mu = \mu_m \frac{S}{K_s + S}$$
  Parameters: Substrate affinity ($K_s$), Maintenance coefficient ($m$), Yield coefficient ($Y$), Max growth rate ($\mu_m$), Dilution rate ($D$), Substrate inflow ($S_{\text{in}}$).

- **Droop Model** (Nutrient Quota Model)  
$$\begin{cases}
    \frac{dX}{dt} = \mu X - D X \\
    \frac{dS}{dt} = \rho X - D S + D S_{\text{in}} \\
    \frac{dQ}{dt} = \rho - \mu Q
  \end{cases}$$
  where  
$$\mu = \mu_m \left(1 - \frac{Q_0}{Q}\right)$$
  and  
$$\rho = \rho_m \frac{S}{K_s + S}$$
  Parameters: Growth rate ($\mu_m$), Nutrient uptake rate ($\rho_m$), Half-saturation ($K_s$), Dilution rate ($D$), Minimum quota ($Q_0$), Substrate inflow ($S_{\text{in}}$).

- **Synthetic Chemostat Model** (Including Biological Inertia)  
$$\begin{cases}
    \frac{dx}{dt} = Y q_s - a_0 r x - D x \\
    \frac{ds}{dt} = D (s_r - s) - q_s x \\
    \frac{dr}{dt} = (Y q_s - a_0 r) \left(\frac{s}{K_r + s} - r\right)
  \end{cases}$$
  where  
$$q_s = r \frac{Q_s K_s}{K_s + s} + (1 - r) \frac{Q_s' K_s'}{K_s' + s}$$
  Parameters: Yield ($Y$), Biological inertia ($a_0$), Dilution rate ($D$), Nutrient uptake coefficients ($Q_s, Q_s'$), Saturation constants ($K_s, K_s'$), Half-saturation constant for $r$($K_r$).

- **Monod-Ierusalimsky** 
This model describes microbial growth, substrate consumption, and product formation using Monod-Ierusalimsky kinetics.  


The specific growth rate $\mu$ follows the Monod-Ierusalimsky kinetics:  
$$\mu = \mu_m \cdot \frac{s}{K_s + s} \cdot \frac{K_p}{K_p + p}$$
where:  
- The first fraction represents substrate-limited growth (Monod equation).  
- The second fraction accounts for product inhibition (Ierusalimsky modification).  

The effective biomass yield is given by:  
$$Y = \frac{Y_{max} \cdot D}{D + m \cdot Y_{max}}$$

Finally, the  System Dynamics is specified by:  

$$\frac{dx}{dt} = \mu x - D x$$

where $\mu x$ represents microbial growth, and $-D x$ accounts for biomass washout due to dilution.  

$$\frac{ds}{dt} = D (s_r - s) - \frac{\mu x}{Y} - m x$$  

$$\frac{dp}{dt} = Y_p \mu x - D p$$

where $Y_p$ is the product yield coefficient, and $-D p$ accounts for product washout.  

The system  state variables are:  
- Biomass concentration: $x = u_1$  
- Substrate concentration: $s = u_2$  
- Product concentration: $p = u_3$  


The system dynamics are governed by the following parameters:  

| Parameter | Description |
|-----------|-------------|
| $K_s$  | Substrate affinity constant |
| $K_p$  | Product inhibition constant |
| $\mu_m$ | Maximum specific growth rate |
| $Y_{max}$ | Maximum yield coefficient |
| $Y_p$  | Product yield coefficient |
| $D$  | Dilution rate |
| $s_r$  | Inflow substrate concentration |
| $m$  | Maintenance coefficient |



| **Model Name**                   | **Parameters List**                                 | **String to Call**                   |
|-----------------------------------|----------------------------------------------------|--------------------------------------|
| SIR Model                        | $\beta, \gamma$                                | `SIR`                           |
| SIR with Birth/Death             | $\beta, \gamma, b, d$                          | `SIR_BD`                        |
| SIS Model                         | $\beta, \gamma$                                | `SIS`                           |
| Lotka-Volterra                   | $\alpha, \beta, \delta, \gamma$                | `Lotka_Volterra`                |
| Lotka-Volterra with Substrate     | $\alpha, \beta, \delta, \gamma, K$            | `Lotka_Volterra_with_substrate`                  |
| Monod Chemostat                  | $K_s, m, Y, \mu_m, D, S_{\text{in}}$          | `Monod_Chemostat`               |
| Droop Model                       | $\mu_m, \rho_m, K_s, D, S_{\text{in}}, Q_0$   | `Droop`                         |
| Synthetic Chemostat               | $Y, a_0, D, Q_s, Q_s', K_s, K_s', K_r$       | `Synthetic_Chemostat`           |
| Monod-Ierusalimsky            | $K_s, K_p, \mu, Y_{max}, Y_p", D, s_r, m$   | `Monod_Ierusalimsky`          |

## Cybernetic models for bacterial growth

In this case, Kinbiont does not take as input an ODE function but a data structure, and it proceeds to construct the ODE system.
Note that the same results in theory could be obtained by declaring a custom ODEs system that behaves in the same way.
The data struct is composed in the following way: 

```julia
model = Kinbiont_Cybernetic_Model(
    Bio_mass_conc = 1.0,  # Initial biomass concentration
    Substrate_concentrations = [3.0, 3.0],  # Concentrations of 2 substrates
    Protein_concentrations = [0.0, 0.0],  # Initial protein concentrations
    allocation_rule = threshold_switching_rule,  # Dynamic resource allocation rule
    reaction = nothing,  # No specific reaction function provided
    cost = nothing,  # No cost function
    protein_thresholds = 0.01,  # Protein activation threshold
    a = [0.1, 0.4],  # Synthesis rate constants for proteins
    b = [0.00001, 0.000001],  # Degradation constants for proteins
    V_S = [0.7, 0.1],  # Substrate utilization rates
    k_S = [0.1, 0.11],  # Saturation constants for substrates
    Y_S = [0.07, 0.11]  # Yield coefficients for biomass per substrate
)
```
Then the ODEs system that will be constructed is the following:
- **Substrate Consumption:**  
$$\frac{dS_i}{dt} = -\frac{V_{S_i} P_i S_i}{k_{S_i} + S_i} \cdot u_1, \quad \forall i \in \{1, \dots, n\}$$
- **Protein Synthesis:**  
$$\frac{dP_i}{dt} = a_i \cdot \text{alloc}_i \cdot k_{S_i} - b_i P_i u_1, \quad \forall i \in \{1, \dots, n\}$$

- **Biomass Growth:**  
$$\frac{dN}{dt} = -\sum_{i=1}^{n} Y_{S_i} \frac{dS_i}{dt} \cdot N$$

Where $n$ is the number of substrates, $V_{S_i}$ i-th substrate utilization rates, $k_{S_i}$ i-th saturation constants for substrates, $Y_{S_i}$ i-th yield coefficients for biomass per substrate, $a_{i}$  i-th synthesis rate for proteins and  $b_{i}$ i-th degradation constants for proteins. 


    
The `allocation_rule` in the cybernetic model defines how resources are allocated to different substrates, proteins, or other components of the system. You can create custom rules based on various factors such as substrate concentration, protein levels, or cost-benefit analysis. Below are some examples of how to define a custom allocation rule:
    
```julia
# Function for proportional allocation based on substrate concentrations
function proportional_allocation_rule(a, b, V_S, k_S, Y_S, P, S, cost, protein_thresholds)
# Normalize substrate concentrations to create allocation vector
alloc = S ./ sum(S)
return alloc
end
```
This rule allocates resources to substrates proportionally based on their concentrations.
For the moment all substrate follow a Monod-like concentration effect on growth rate, so please leave the option `reaction = nothing`. 

## Reaction networks

In the case the user want to input the system as a network of reaction, Kinbiont.jl relies on [Catalyst.jl](https://github.com/SciML/Catalyst.jl) to generate the ODE problem to be fitted and solved.
A network and its parameters should be declared as the following:
```julia
u0 = [:S => 301, :E => 100, :SE => 0, :P => 0]
ps = [:kB => 0.00166, :kD => 0.0001, :kP => 0.1]

model_Michaelis_Menten = @reaction_network begin
    kB, S + E --> SE
    kD, SE --> S + E
    kP, SE --> P + E
end 
```
For other examples on how declare a reaction network please consult:[Catalyst.jl documentation](https://docs.sciml.ai/Catalyst/stable/).
