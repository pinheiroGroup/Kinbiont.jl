using Documenter, Kinbiont

push!(LOAD_PATH,"../src/")
#Documenter.Remotes.URL("../src/", nothing)

makedocs(sitename="Kinbiont Documentation",
        modules = [Kinbiont],
        source  = "src/",
        checkdocs=:none,
        format = Documenter.HTML(
          size_threshold = nothing,
          assets = ["assets/favicon.ico"]
        )
)
