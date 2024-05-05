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

5. To activate the JMAKi project, type
> instantiate

6. at the start of the code or notebook (where you are going to do analyses) you should write 

```
using DifferentialEquations
using Optimization
using Plots
using Random
using CSV
using DataFrames
using Statistics
using Optim
using OptimizationBBO
using NaNMath
using StatsBase
using Tables
using Distributions
using Interpolations
using Peaks
using ChangePointDetection
using Lowess  
using  LsqFit
using Combinatorics
using MLJ
using SymbolicRegression
using Zygote
using DelimitedFiles
using BetaML
using DecisionTree
using MLJDecisionTreeInterface
include("your_path_to_JMAKi_main_folder/src/functions.jl")


```
this last step is Temporary before the official realese
## Package installation
## Requirements
### Dependencies

- Julia (1.7,1.8,1.9)
- DifferentialEquations
- Optimization
- Plots
- Random
- CSV
- DataFrames
- Statistics
-  Optim
-   OptimizationBBO
- NaNMath
- StatsBase
- Tables
- Distributions
- Interpolations
- Peaks
- ChangePointDetection
- Lowess
- LsqFit
- Combinatorics


