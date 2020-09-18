using Uchiyama
using Documenter

makedocs(;
    modules=[Uchiyama],
    authors="Pierre Navaro",
    repo="https://github.com/pnavaro/Uchiyama.jl/blob/{commit}{path}#L{line}",
    sitename="Uchiyama.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://pnavaro.github.io/Uchiyama.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/pnavaro/Uchiyama.jl",
)
