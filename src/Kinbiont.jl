module Kinbiont
  include("functions.jl")

  # Unified matrix API
  include("api/types.jl")
  include("api/model_registry.jl")
  include("api/preprocessing.jl")
  include("api/fitting.jl")
  include("api/io.jl")
end
