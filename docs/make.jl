using SIMIlluminationPatterns
using Documenter, DocumenterCitations

# NOTE: When updating, must update also in `test/runtests.jl` <18-12-24> 
DocMeta.setdocmeta!(SIMIlluminationPatterns, :DocTestSetup, :(
        using SIMIlluminationPatterns;
        using SIMIlluminationPatterns.Synthetic;
        using Distributions;
        using TestImages
    ); recursive=true)

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
            "Harmonic" => "pages/03_patterns/01_harmonic.md",
            "Nonlinear SIM" => "pages/03_patterns/02_nonlinear_sim.md",
            "Blind SIM" => "pages/03_patterns/03_blindsim.md",
        ],
        "Synthetic Data" => "pages/04_synthetic_data.md",
        "References" => [
            "API" => "pages/05_apireference.md",
            "Bibliography" => "pages/06_bibliography.md",
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
