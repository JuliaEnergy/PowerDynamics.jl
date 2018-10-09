using Documenter
using PowerDynBase

makedocs(
    # options
    modules = [PowerDynBase],
    # html options
    format = :html,
    sitename = "PowerDynBase.jl",
    pages = ["index.md"])
