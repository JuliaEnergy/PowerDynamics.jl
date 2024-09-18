using OpPoDyn
using Documenter

DocMeta.setdocmeta!(OpPoDyn, :DocTestSetup, :(using OpPoDyn); recursive=true)

makedocs(;
    modules=[OpPoDyn],
    authors="Hans Würfel <git@wuerfel.io> and contributors",
    sitename="OpPoDyn.jl",
    linkcheck=true,
    pagesonly=true,
    format=Documenter.HTML(;
        canonical="https://juliaenergy.github.io/OpPoDyn.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Modeling Concepts" => "ModelingConcepts.md",
    ],
)

deploydocs(;
    repo="github.com/JuliaEnergy/OpPoDyn.jl",
    devbranch="main",
    push_preview=true,
)
