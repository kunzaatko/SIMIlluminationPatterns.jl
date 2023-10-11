using SIMIlluminationPatterns
using Documenter

DocMeta.setdocmeta!(SIMIlluminationPatterns, :DocTestSetup, :(using SIMIlluminationPatterns); recursive=true)

makedocs(;
    modules=[SIMIlluminationPatterns],
    authors="Martin Kunz <martinkunz@email.cz> and contributors",
    repo="https://github.com/kunzaatko/SIMIlluminationPatterns.jl/blob/{commit}{path}#{line}",
    sitename="SIMIlluminationPatterns.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://kunzaatko.github.io/SIMIlluminationPatterns.jl",
        edit_link="trunk",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/kunzaatko/SIMIlluminationPatterns.jl",
    devbranch="trunk",
)
