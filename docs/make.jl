# using Pkg
# Pkg.develop(PackageSpec(path=pwd()))
# Pkg.instantiate()
# cd(@__DIR__)
# Pkg.activate(".")
# push!(LOAD_PATH,"../src/")
using PhaseRetrieval
@show mm # just to test that the package is loaded
using Documenter, Literate

DocMeta.setdocmeta!(PhaseRetrieval, :DocTestSetup, :(using PhaseRetrieval); recursive=true)

# TODO Generate tutorials from literate files
@info "current dir =$(@__DIR__)"
tutorials_folder = (@__DIR__) * "/../tutorials"
docs_tutorials_folder = (@__DIR__) * "/src/tutorials"
@info tutorials_folder
for f in readdir(tutorials_folder; join=true)
    Literate.markdown(f, docs_tutorials_folder)
end

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
    # draft=true,
    checkdocs=:exported,
    # doctest=:fix,
    warnonly=true,
    format=Documenter.HTML(;
        # Use clean URLs, unless built as a "local" build
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://olejorik.github.io/PhaseRetrieval.jl/stable/",
        assets=["assets/favicon.ico"],
        highlights=["yaml"],
    ),
    clean=false,
    pages=[
        "Home" => "index.md",
        "About" => "about.md",
        "Getting Started" => "tutorials/GettingStarted.md",
        "Manual" => [
            "Forward model" => "tutorials/Hardware.md",
            "Inverse Problem" => [
                "Introduction" => "Inverse.md",
                "Gonsalves's method" => "tutorials/Gonsalves.md",
                "AP-based methods" => "tutorials/AP.md",
            ],
        ],
    ],
)

deploydocs(; repo="github.com/olejorik/PhaseRetrieval.jl.git", target="build")
