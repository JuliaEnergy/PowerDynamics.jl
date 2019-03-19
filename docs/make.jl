
using Documenter
using PowerDynBase, PowerDynSolve

makedocs(
    modules = [PowerDynBase, PowerDynSolve],
    authors = "Tim Kittel and contributors.",
    linkcheck=true,
    sitename = "PowerDynamics.jl",
    pages = [
        "General" => "index.md",
        "language_conventions.md",
        "node_dynamics_types.md",
        "node_types.md",
        "custom_node_types.md",
        "states_solutions.md",
        "error_types.md",
        "fullindex.md",
        "contact.md",
    ],
   format = Documenter.HTML(canonical = "https://juliaenergy.github.io/PowerDynamics.jl/stable/"),
   )
deploydocs(
    repo = "github.com/JuliaEnergy/PowerDynamics.jl.git",
    target = "build"
)
