using PowerDynamics
using PowerDynamics.Library
using Documenter
using Literate
using DocumenterInterLinks

links = InterLinks(
    "NetworkDynamics" => "https://juliadynamics.github.io/NetworkDynamics.jl/stable/",
)

DocMeta.setdocmeta!(PowerDynamics, :DocTestSetup, :(using PowerDynamics); recursive=true)

# generate examples
example_dir = joinpath(@__DIR__, "examples")
outdir = joinpath(@__DIR__, "src", "generated")
isdir(outdir) && rm(outdir, recursive=true)
mkpath(outdir)

for example in filter(contains(r".jl$"), readdir(example_dir, join=true))
    Literate.markdown(example, outdir)
    Literate.script(example, outdir; keep_comments=true)
end

kwargs = (;
    modules=[PowerDynamics, PowerDynamics.Library],
    authors="Hans Würfel <git@wuerfel.io> and contributors",
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
        "API.md",
        "Examples" => [
            "generated/ieee9bus.md",]
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

    deploydocs(; repo="github.com/JuliaEnergy/PowerDynamics.jl",
            devbranch="main", push_preview=true)

    success || throw(thrown_ex)
else # local build
    makedocs(; kwargs_warnonly...)
end
