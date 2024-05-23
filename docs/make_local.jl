using Documenter, Kimchi
push!(LOAD_PATH,"../src/")
Documenter.Remotes.URL("../src/", nothing)
makedocs(sitename="Kimchi Documentation",
        modules = [Kimchi],
        source  = "src/",
        checkdocs=:none,
        remotes = nothing,
        format = Documenter.HTML(
        assets = ["assets/favicon.ico"]
        )
)