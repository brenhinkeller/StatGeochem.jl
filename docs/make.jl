using StatGeochemBase
using Documenter

DocMeta.setdocmeta!(StatGeochemBase, :DocTestSetup, :(using StatGeochemBase); recursive=true)

makedocs(;
    modules=[StatGeochemBase],
    authors="C. Brenhin Keller",
    repo="https://github.com/brenhinkeller/StatGeochemBase.jl/blob/{commit}{path}#{line}",
    sitename="StatGeochemBase.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://brenhinkeller.github.io/StatGeochemBase.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/brenhinkeller/StatGeochemBase.jl",
    devbranch = "main",
)
