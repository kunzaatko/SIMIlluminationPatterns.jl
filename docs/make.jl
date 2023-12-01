using SIMIlluminationPatterns
using Documenter, DocumenterCitations

DocMeta.setdocmeta!(SIMIlluminationPatterns, :DocTestSetup, :(using SIMIlluminationPatterns); recursive=true)

bib = CitationBibliography(
    joinpath(@__DIR__, "src", "refs.bib");
    # style=:authoryear
)

makedocs(;
    modules=[SIMIlluminationPatterns],
    authors="Martin Kunz <martinkunz@email.cz> and contributors",
    repo="https://github.com/kunzaatko/SIMIlluminationPatterns.jl/blob/{commit}{path}#{line}",
    sitename="SIMIlluminationPatterns.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://kunzaatko.github.io/SIMIlluminationPatterns.jl",
        edit_link="trunk",
        assets=String[]
    ),
    pages=[
        "Home" => "index.md",
        "Theory" => "pages/01_theory.md",
        "General Interface" => "pages/02_interface.md",
        "Illumination Patterns" => [
            "Harmonic" => "pages/03_patterns/01_harmonic.md"
            "Nonlinear SIM" => "pages/03_patterns/02_nonlinear_sim.md"
            "Blind SIM" => "pages/03_patterns/03_blindsim.md"
        ],
        "References" => [
            "API" => "pages/04_apireference.md",
            "Bibliography" => "pages/05_bibliography.md"
        ]
    ],
    plugins=[bib],
    # NOTE: doctesting is done in the `runtests.jl` so it is not necessary to do here
    doctest=false
)

deploydocs(;
    repo="github.com/kunzaatko/SIMIlluminationPatterns.jl",
    devbranch="trunk"
)
