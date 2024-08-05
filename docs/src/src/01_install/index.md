# [Installation & requirements](@id installation)

## Manual installation
Download the folder from Github. 
First, it is necessary to install the Julia from https://julialang.org/.   
Next, the user need to copy the project folder in the chosen working directory. 

1. Using REPL or the COMMAND LINE move to the working directory.  
2. If you use the COMMAND LINE, to start a Julia session run the command:

> julia

3. To enter in the Pkg REPL  type 

>]  

4. Type the command 
> activate .

5. To activate the KinBiont project, type
> instantiate

6. at the start of the code or notebook (where you are going to do analyses) you should write 

```
using KinBiont
```

## Package installation

With any package manager run
```
Pkg.add(KinBiont)
```
and then in your script use:

```
using KinBiont
```
## Requirements & Dependencies

- **Julia**: (>1.10)
- **AbstractTrees**: 0.4.5
- **BetaML**: 0.12.0
- **CSV**: 0.10.14
- **ChangePointDetection**: 1.2.0
- **Combinatorics**: 1.0.2
- **DataFrames**: 1.6.1
- **DecisionTree**: 0.12.4
- **DelimitedFiles**: 1.9.1
- **DifferentialEquations**: 7.13.0
- **Distributions**: 0.25.109
- **Interpolations**: 0.13.6
- **Lowess**: 0.1.0
- **LsqFit**: 0.15.0
- **MLJ**: 0.20.6
- **MLJDecisionTreeInterface**: 0.4.2
- **Missings**: 1.2.0
- **NaNMath**: 1.0.2
- **Optim**: 1.9.4
- **Optimization**: 3.26.3
- **OptimizationBBO**: 0.3.0
- **OptimizationMultistartOptimization**: 0.2.0
- **OrdinaryDiffEq**: 6.85.0
- **Peaks**: 0.5.2
- **Random**: 1.10.0
- **SciMLBase**: 2.42.0
- **SpecialFunctions**: 2.4.0
- **Statistics**: 1.10.0
- **StatsBase**: 0.34.3
- **SymbolicRegression**: 0.22.2
- **SymbolicUtils**: 1.6.0
- **Tables**: 1.11.1
- **Zygote**: 0.6.70



