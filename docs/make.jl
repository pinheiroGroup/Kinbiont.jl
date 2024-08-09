using Documenter, Kinbiont

makedocs(sitename="Kinbiont Documentation",
        modules = [Kinbiont],
        checkdocs=:none,
        format = Documenter.HTML(
            assets = ["assets/favicon.ico"]
        )
)
