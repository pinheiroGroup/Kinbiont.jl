module Kinbiont
  include("functions.jl")

  # Unified matrix API
  include("api/types.jl")
  include("api/irregular.jl")
  include("api/model_registry.jl")
  include("api/preprocessing.jl")
  include("api/clustering_helpers.jl")
  include("api/fitting.jl")
  include("api/io.jl")
  include("api/batch.jl")
end
