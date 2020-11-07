using QWNNM
using Documenter

makedocs(;
    modules=[QWNNM],
    authors="ramsterdam91",
    repo="https://github.com/ramsterdam91/QWNNM.jl/blob/{commit}{path}#L{line}",
    sitename="QWNNM.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://ramsterdam91.github.io/QWNNM.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/ramsterdam91/QWNNM.jl",
)
