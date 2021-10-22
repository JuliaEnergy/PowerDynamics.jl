#=
julia localmake.jl
=#

using Pkg
Pkg.activate(@__DIR__)
Pkg.develop(PackageSpec(path=dirname(@__DIR__))) # adds the package this script is called from
push!(LOAD_PATH, Base.Filesystem.abspath("../"))
using Documenter
using PowerDynamics

include("make.jl")

using LiveServer; serve(dir=joinpath(@__DIR__, "build"))