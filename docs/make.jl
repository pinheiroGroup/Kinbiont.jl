using Documenter, JMAKi

makedocs(sitename="JMAKi Documentation",
        modules = [JMAKi],
        checkdocs=:none,
        format = Documenter.HTML(
            assets = ["assets/favicon.ico"]
        )
)
