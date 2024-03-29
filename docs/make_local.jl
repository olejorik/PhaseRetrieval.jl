using Pkg
Pkg.add(url = "https://github.com/olejorik/AlternatingProjections.jl", rev= "develop"); 
Pkg.add(url = "https://github.com/olejorik/SampledDomains.jl");
Pkg.develop(PackageSpec(path=pwd()))
Pkg.instantiate()
# cd(@__DIR__)
# Pkg.activate(".")
# push!(LOAD_PATH,"../src/")
using PhaseRetrieval
@show mm # just to test that the package is loaded
using Documenter

DocMeta.setdocmeta!(PhaseRetrieval, :DocTestSetup, :(using PhaseRetrieval); recursive=true)


makedocs(;
    modules=[PhaseRetrieval],
    authors="Oleg Soloviev",
    repo="https://github.com/olejorik/PhaseRetrieval.jl/blob/{commit}{path}#L{line}",
    sitename="PhaseRetrieval.jl",
    # format=Documenter.HTML(;
        # prettyurls=get(ENV, "CI", "false") == "true",
        # canonical="https://olejorik.github.io/PhaseRetrieval.jl",
        # assets=String[],
    # ),
    format = Documenter.HTML(
        # Use clean URLs, unless built as a "local" build
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical = "https://olejorik.github.io/PhaseRetrieval.jl/stable/",
        assets = ["assets/favicon.ico"],
        highlights = ["yaml"],
    ),
    clean = false,
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/olejorik/PhaseRetrieval.jl.git",
    target = "build",
)

Pkg.rm("AlternatingProjections")
Pkg.rm("SampledDomains")
Pkg.rm("PhaseRetrieval")