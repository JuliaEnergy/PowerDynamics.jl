
using Documenter
using PowerDynBase, PowerDynSolve

makedocs(
    # options
    modules = [PowerDynBase, PowerDynSolve],
    # html options
    format = :html,
    sitename = "PowerDynamics.jl",
    pages = [
        "Dynamic Power System Analysis" => [
            "PowerDynamics.jl" => "index.md",
            "language_conventions.md",
            "node_dynamics_types.md",
            "node_types.md",
            "custom_node_types.md",
            "error_types.md",
        ],
        "Documentation of Internals" => [
            "internalsBase.md",
            "internalsSolve.md",
        ],
        "fullindex.md",
        "contact.md"
    ])
