using Documenter, KinBiont
push!(LOAD_PATH,"../src/")
Documenter.Remotes.URL("../src/", nothing)
makedocs(sitename="KinBiont Documentation",
        modules = [KinBiont],
        source  = "src/",
        checkdocs=:none,
        remotes = nothing,
        format = Documenter.HTML(
        assets = ["assets/favicon.ico"]
        )
)