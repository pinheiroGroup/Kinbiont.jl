using Documenter, Kimchi

makedocs(sitename="Kimchi Documentation",
        modules = [Kimchi],
        checkdocs=:none,
        format = Documenter.HTML(
            assets = ["assets/favicon.ico"]
        )
)
