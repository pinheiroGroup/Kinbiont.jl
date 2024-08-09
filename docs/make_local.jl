using Documenter, Kinbiont
push!(LOAD_PATH,"../src/")
Documenter.Remotes.URL("../src/", nothing)
makedocs(sitename="Kinbiont Documentation",
        modules = [Kinbiont],
        source  = "src/",
        checkdocs=:none,
        remotes = nothing,
        format = Documenter.HTML(
        assets = ["assets/favicon.ico"]
        )
)