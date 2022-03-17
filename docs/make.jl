using QUBOConstraints
using Documenter

DocMeta.setdocmeta!(QUBOConstraints, :DocTestSetup, :(using QUBOConstraints); recursive=true)

makedocs(;
    modules=[QUBOConstraints],
    authors="azzaare <jf@baffier.fr> and contributors",
    repo="https://github.com/JuliaConstraints/QUBOConstraints.jl/blob/{commit}{path}#{line}",
    sitename="QUBOConstraints.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://JuliaConstraints.github.io/QUBOConstraints.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/JuliaConstraints/QUBOConstraints.jl",
    devbranch="main",
)
