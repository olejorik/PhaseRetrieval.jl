using PhaseRetrieval
using Documenter

makedocs(;
    modules=[PhaseRetrieval],
    authors="Oleg Soloviev",
    repo="https://github.com/olejorik/PhaseRetrieval.jl/blob/{commit}{path}#L{line}",
    sitename="PhaseRetrieval.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://olejorik.github.io/PhaseRetrieval.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/olejorik/PhaseRetrieval.jl",
)
