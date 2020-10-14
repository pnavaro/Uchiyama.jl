using Uchiyama
using Documenter

using Literate
using Plots

# generate examples
examples = ["event_driven.jl", "hard_spheres.jl", 
            "periodic_hard_spheres.jl", "periodic_uchiyama.jl"]
output = joinpath(@__DIR__, "src")
examples_dir = joinpath(@__DIR__, "..", "examples")
for example in examples
    jl_file = joinpath(examples_dir, example)
    Literate.markdown(jl_file, output, documenter=true)
    Literate.notebook(jl_file, output, execute=false)
    #Literate.script(EXAMPLE, OUTPUT)
end

makedocs(;
    modules=[Uchiyama],
    authors="Nathalie Ayi and Pierre Navaro",
    repo="https://github.com/pnavaro/Uchiyama.jl/blob/{commit}{path}#L{line}",
    sitename="Uchiyama.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://pnavaro.github.io/Uchiyama.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Uchiyama model in box" => "event_driven.md",
        "Uchiyama model periodic" => "periodic_uchiyama.md",
        "Hard spheres in box" => "hard_spheres.md",
        "Hard spheres periodic" => "periodic_hard_spheres.md",
    ],
)

deploydocs(;
    repo="github.com/pnavaro/Uchiyama.jl",
)
