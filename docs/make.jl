using BatchReactor
using Documenter

DocMeta.setdocmeta!(BatchReactor, :DocTestSetup, :(using BatchReactor); recursive=true)

makedocs(;
    modules=[BatchReactor],
    authors="Vinod Janardhanan",
    repo="https://github.com/vinodjanardhanan/BatchReactor.jl/blob/{commit}{path}#{line}",
    sitename="BatchReactor.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://vinodjanardhanan.github.io/BatchReactor.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/vinodjanardhanan/BatchReactor.jl",
    devbranch="main",
)
