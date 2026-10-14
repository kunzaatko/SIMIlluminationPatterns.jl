using SIMIlluminationPatterns
using Documenter, DocumenterCitations, DocumenterInterLinks, MakieMaestro

MakieMaestro.Themes.width!(25u"cm")

makie_doc_blocks = MakieMaestro.MakieDocBlocks(;
    formats=[:png]
)

# NOTE: Links can be explored trough the REPL with `links(query)`
links = InterLinks(
    "Julia" => "https://docs.julialang.org/en/v1/",
    "Unitful" => "https://juliaphysics.github.io/Unitful.jl/dev/",
    "TransferFunctions" => "https://kunzaatko.github.io/TransferFunctions.jl/stable/",
    "ComponentArrays" => "https://docs.sciml.ai/ComponentArrays/stable/",
    "DomainSets" => "https://juliaapproximation.github.io/DomainSets.jl/dev/",
    "StaticArrays" => "https://juliaarrays.github.io/StaticArrays.jl/stable/"
)

DocMeta.setdocmeta!(SIMIlluminationPatterns, :DocTestSetup, :(
        include(joinpath(@__DIR__, "../test/doctestsetup.jl"))
    ); recursive=true)

bib = CitationBibliography(
    joinpath(@__DIR__, "src", "refs.bib");
    # style=:authoryear
)

makedocs(;
    modules=[SIMIlluminationPatterns],
    authors="Martin Kunz <martinkunz@email.cz> and contributors",
    repo=Remotes.GitHub("kunzaatko", "SIMIlluminationPatterns.jl"),
    sitename="SIMIlluminationPatterns.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://kunzaatko.github.io/SIMIlluminationPatterns.jl",
        edit_link="trunk",
        assets=String[]
        # assets=["assets/favicon.ico"],
    ),
    pages=[
        "Home" => "index.md",
        "Manual" => "manual.md",
        "Harmonic" => "harmonic.md",
        "References" => [
            "API" => "apireference.md",
            "Bibliography" => "bibliography.md",
        ]
    ],
    warnonly=[:missing_docs],
    plugins=[bib, links],
    doctest=false # tests run in `test/runtests.jl`
)

deploydocs(;
    repo="github.com/kunzaatko/SIMIlluminationPatterns.jl",
    devbranch="trunk"
)

if !haskey(ENV, "GITHUB_ACTIONS")
    build_path = joinpath(@__DIR__, "build")
    cached_path = joinpath(@__DIR__, "cached")
    @info "Making cached docs at $cached_path"
    cp(build_path, cached_path; force=true)
end
