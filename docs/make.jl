using Documenter
using PowerDynamics
using BlockSystems
using Literate

DocMeta.setdocmeta!(PowerDynamics, :DocTestSetup, :(using PowerDynamics); recursive=true)

Literate.markdown(joinpath(@__DIR__, "src", "MARiE.jl"), joinpath(@__DIR__, "src", "generated"))

makedocs(
    modules = [PowerDynamics],
    authors = "Tim Kittel, Jan Liße, Sabine Auer and further contributors.",
    linkcheck=false,
    sitename = "PowerDynamics.jl",
    pages = [
        "General" => "index.md",
        "Architecture" => "ARCHITECTURE.md",
        "language_conventions.md",
        "powergrid_model.md",
        "node_types.md",
        "custom_node_types.md",
        "Modular Inverters" => "generated/MARiE.md",
        "line_types.md",
        "states_solutions.md",
        "simulations.md",
        "error_types.md",
        "import_export.md",
        "fullindex.md",
        "contact.md",
    ],
    format = Documenter.HTML(canonical = "https://juliaenergy.github.io/PowerDynamics.jl/stable/"),
    warnonly=[:missing_docs]
   )
deploydocs(
    repo = "github.com/JuliaEnergy/PowerDynamics.jl.git",
    devbranch = "main"
)
