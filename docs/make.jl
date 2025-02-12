using OpPoDyn
using Documenter
using Literate
using DocumenterInterLinks

links = InterLinks(
    "NetworkDynamics" => "https://juliadynamics.github.io/NetworkDynamics.jl/stable/",
)

DocMeta.setdocmeta!(OpPoDyn, :DocTestSetup, :(using OpPoDyn); recursive=true)

# generate examples
example_dir = joinpath(@__DIR__, "examples")
outdir = joinpath(@__DIR__, "src", "generated")
isdir(outdir) && rm(outdir, recursive=true)
mkpath(outdir)

for example in filter(contains(r".jl$"), readdir(example_dir, join=true))
    Literate.markdown(example, outdir)
    Literate.script(example, outdir; keep_comments=true)
end


makedocs(;
    modules=[OpPoDyn],
    authors="Hans Würfel <git@wuerfel.io> and contributors",
    sitename="OpPoDyn.jl",
    linkcheck=true,
    pagesonly=true,
    plugins=[links],
    format=Documenter.HTML(;
        canonical="https://juliaenergy.github.io/OpPoDyn.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Modeling Concepts" => "ModelingConcepts.md",
        "Examples" => [
            "generated/ieee9bus.md",]
    ],
    warnonly=[:cross_references, :missing_docs, :docs_block],
)

deploydocs(;
    repo="github.com/JuliaEnergy/OpPoDyn.jl",
    devbranch="main",
    push_preview=true,
)
