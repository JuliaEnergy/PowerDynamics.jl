using PowerDynamics
using PowerDynamics.Library
using Documenter
using Literate
using DocumenterInterLinks
using ModelingToolkit
using Graphs
using NetworkDynamics
using OrdinaryDiffEqRosenbrock, OrdinaryDiffEqNonlinearSolve
using CairoMakie

links = InterLinks(
    "NetworkDynamics" => "https://juliadynamics.github.io/NetworkDynamics.jl/stable/",
)

DocMeta.setdocmeta!(PowerDynamics, :DocTestSetup, :(using PowerDynamics); recursive=true)

# generate examples
example_dir = joinpath(@__DIR__, "examples")
tutorial_dir = joinpath(@__DIR__, "tutorials")
outdir = joinpath(@__DIR__, "src", "generated")
isdir(outdir) && rm(outdir, recursive=true)
mkpath(outdir)

for dir in (example_dir, tutorial_dir)
    for script in filter(contains(r".jl$"), readdir(dir, join=true))
        Literate.markdown(script, outdir)
        Literate.script(script, outdir; keep_comments=true)
    end
end

kwargs = (;
    modules=[PowerDynamics, PowerDynamics.Library],
    authors="Hans Würfel, Tim Kittel, Jan Liße, Sabine Auer, Anton Plietzsch and contributors",
    sitename="PowerDynamics.jl",
    linkcheck=true,
    pagesonly=true,
    plugins=[links],
    format=Documenter.HTML(;
        canonical="https://juliaenergy.github.io/PowerDynamics.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Modeling Concepts" => "ModelingConcepts.md",
        "Initialization" => "initialization.md",
        # "Component Library" => "Library.md",
        "Tutorials" => [
            # "generated/custom_bus.md",
        ],
        "Examples" => [
            # "generated/ieee9bus.md",
            "IEEE39 Part I: Modeling" => "generated/ieee39_part1.md",
            "IEEE39 Part II: Initialization" => "generated/ieee39_part2.md",
        ],
        "API" => "API.md",
        "🔗 NetworkDynamics.jl Docs" => "networkdynamics_forward.md",
    ],
    warnonly=[:missing_docs],
)
kwargs_warnonly = (; kwargs..., warnonly=true)

if haskey(ENV,"GITHUB_ACTIONS")
    success = true
    thrown_ex = nothing
    try
        makedocs(; kwargs...)
    catch e
        @info "Strict doc build failed, try again with warnonly=true"
        global success = false
        global thrown_ex = e
        makedocs(; kwargs_warnonly...)
    end

    deploydocs(; repo="github.com/JuliaEnergy/PowerDynamics.jl.git",
            devbranch="main", push_preview=true)

    success || throw(thrown_ex)
else # local build
    makedocs(; kwargs_warnonly...)
end
