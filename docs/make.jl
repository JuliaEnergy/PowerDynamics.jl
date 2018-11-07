
using Documenter
using PowerDynBase, PowerDynSolve

makedocs(
    # options
    modules = [PowerDynBase, PowerDynSolve],
    # html options
    format = :html,
    sitename = "PowerDynamics.jl",
    pages = [
        "General" => "index.md",
        "language_conventions.md",
        "node_dynamics_types.md",
        "node_types.md",
        "custom_node_types.md",
        "states_solutions.md",
        "error_types.md",
        # "Documentation of Internals" => [
        #     "internalsBase.md",
        #     "internalsSolve.md",
        # ],
        "fullindex.md",
        "contact.md",
    ],
   html_canonical = "https://juliaenergy.github.io/PowerDynamics.jl/latest/",
   )
deploydocs(
    repo = "github.com/JuliaEnergy/PowerDynamics.jl.git",
    target = "build",
    julia = "1.0",
    deps = nothing,
    make = nothing
)
