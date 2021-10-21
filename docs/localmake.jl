#=
julia localmake.jl
=#

using Pkg
Pkg.activate(@__DIR__)
Pkg.develop(PackageSpec(path=dirname(@__DIR__))) # adds the package this script is called from
push!(LOAD_PATH, Base.Filesystem.abspath("../"))
using Documenter
using PowerDynamics

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
        "line_types.md",
        "states_solutions.md",
        "simulations.md",
        "error_types.md",
        "import_export.md",
        "fullindex.md",
        "contact.md",
    ]
   )