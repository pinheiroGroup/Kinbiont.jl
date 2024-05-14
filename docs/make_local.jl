using Documenter, JMAKi
push!(LOAD_PATH,"../src/")
Documenter.Remotes.URL("../src/", nothing)
makedocs(sitename="JMAKi Documentation",
        modules = [JMAKi],
        source  = "src/",
        checkdocs=:none,
        remotes = nothing,
        format = Documenter.HTML(
        assets = ["assets/favicon.ico"]
        )
)