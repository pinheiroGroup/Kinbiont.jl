using Documenter, Kinbiont

makedocs(sitename="Kinbiont Documentation",
        modules = [Kinbiont],
        checkdocs=:none,
        format = Documenter.HTML(
            size_threshold = nothing,
            assets = ["assets/favicon.ico"]
        )
)
