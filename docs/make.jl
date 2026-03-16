using Documenter, DocumenterCitations, Helios

bib = CitationBibliography(joinpath(@__DIR__, "src", "references.bib"))

makedocs(;
    sitename="Helios.jl",
    pages = [
        "Home" => "index.md",
        "Usage" => [
            "Getting started" => "usage/getting-started.md"
        ],
        "API" => [
            "api/solar_position.md",
            "api/atmosphere.md",
            "api/irradiance.md",
            "api/clearsky.md",
        ],
        "References" => "references.md"
    ],
    format = Documenter.HTML(
        edit_link="master",
        assets=["assets/favicon.ico"]
    ),
    repo=Remotes.GitHub("nicomignoni", "Helios.jl"),
    plugins=[bib]
)

deploydocs(
    repo = "github.com/nicomignoni/Helios.jl.git",
)
