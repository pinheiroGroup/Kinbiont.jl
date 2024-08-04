using Documenter, KinBiont

makedocs(sitename="KinBiont Documentation",
        modules = [KinBiont],
        checkdocs=:none,
        format = Documenter.HTML(
            assets = ["assets/favicon.ico"]
        )
)
