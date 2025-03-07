# KinBiont.jl Development Guidelines

## Build & Test Commands
- Run Julia REPL in project: `julia --project=.`
- Build documentation: `julia --project=docs --color=yes docs/make.jl`
- Run examples: `julia --project=. Example_scripts/examples_single_fit.jl`

## Code Style Guidelines
- **Imports**: Use `using` for packages, `include()` for local files
- **Types**: Type function parameters when possible, avoid overuse of `Any`
- **Naming**: Use `snake_case` for functions and variables
- **Documentation**: Use triple-quote docstrings describing function, args, and returns
- **Formatting**: 4-space indentation, reasonable line length (<100 chars)
- **Error handling**: Prefer `@warn` for non-fatal issues, check input bounds
- **Module Structure**: Export public functions with `export` statements

## Example Pattern
```julia
function process_data(input_data::DataFrame, threshold::Float64=0.5)
    # Function implementation...
    return processed_data
end
```