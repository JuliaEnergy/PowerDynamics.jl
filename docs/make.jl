using PowerDynamics
using PowerDynamics.Library
using Documenter
using Literate
using DocumenterInterLinks
using ModelingToolkitBase
using Graphs
using NetworkDynamics
using OrdinaryDiffEqRosenbrock, OrdinaryDiffEqNonlinearSolve
using CairoMakie

links = InterLinks(
    "NetworkDynamics" => "https://juliadynamics.github.io/NetworkDynamics.jl/stable/",
    "DiffEq" => "https://docs.sciml.ai/DiffEqDocs/stable/",
    "ModelingToolkit" => "https://docs.sciml.ai/ModelingToolkit/stable/",
)
#=
# manually look for interlinks for debuggin purposes
links["NetworkDynamics"]("get_callbacks") # search for name in all
links["NetworkDynamics"]("Sparsity") # search for name in all
=#

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

doc = makedocs(;
    modules=[PowerDynamics, PowerDynamics.Library, PowerDynamics.Library.ComposableInverter],
    authors="Hans Würfel, Tim Kittel, Jan Liße, Sabine Auer, Anton Plietzsch and contributors",
    sitename="PowerDynamics.jl",
    build=haskey(ENV, "DOCUMENTER_DRAFT") ? "build_draft" : "build",
    linkcheck=true,
    pagesonly=true,
    plugins=[links],
    format=Documenter.HTML(;
        canonical="https://juliaenergy.github.io/PowerDynamics.jl",
        edit_link="main",
        assets=String["assets/custom.css"],
    ),
    pages=[
        "Home" => "index.md",
        "Modeling Concepts" => "ModelingConcepts.md",
        "Initialization" => "initialization.md",
        "Component Library" => "Library.md",
        "vs_and_cs_models.md",
        "Tutorials" => [
            "Julia Setup for New Users" => "julia_setup.md",
            "Getting Started" => "generated/getting_started.md",
            "Typical Simulation Workflow" => "generated/typical_simulation_workflow.md",
            "Custom Generator Bus" => "generated/custom_bus.md",
            "Custom Transmission Line" => "generated/custom_line.md",
        ],
        "Advanced Examples" => [
            "generated/ieee9bus.md",
            "IEEE39 Part I: Modeling" => "generated/ieee39_part1.md",
            "IEEE39 Part II: Initialization" => "generated/ieee39_part2.md",
            "IEEE39 Part III: Simulation" => "generated/ieee39_part3.md",
            "IEEE39 Part IV: Parameter Tuning" => "generated/ieee39_part4.md",
            "Linear Analysis of 4-Bus EMT System" => "generated/linear_analysis.md",
            "EMT Toymodel" => "generated/emt_toymodel.md",
            "Zero-Impedance Circuit Breaker" => "generated/zero_imp_breaker.md",
        ],
        "API" => "API.md",
        "🔗 NetworkDynamics.jl Docs" => "networkdynamics_forward.md",
    ],
    draft=haskey(ENV, "DOCUMENTER_DRAFT"),
    linkcheck_ignore = [
        r"^\.\./assets/OpenIPSL_valid/.*\.png$",  # Match ../assets/OpenIPSL_valid/*.png
        "https://marketplace.visualstudio.com/items?itemName=julialang.language-julia", # curl blocked?
        ],
    warnonly=true,
    debug=true, # return doc object
)

if haskey(ENV, "GITHUB_ACTIONS")
    deploydocs(; repo="github.com/JuliaEnergy/PowerDynamics.jl.git",
            devbranch="main", push_preview=true)
    errors = setdiff(doc.internal.errors, [:missing_docs])
    if !isempty(errors)
        error("makedocs encountered errors: $(errors)")
    end
end
