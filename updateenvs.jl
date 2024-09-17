#!julia --startup-file=no

using Pkg
cd(@__DIR__)
Pkg.activate(@__DIR__)
@assert Pkg.TOML.parsefile(Base.active_project())["name"] == "OpPoDyn"

subp = ["OpPoDyntesting", "test", "docs"]
for p in subp
    @info "Activate and update $p environment"
    Pkg.activate(p)
    Pkg.update()
    Pkg.resolve()
end
Pkg.activate(@__DIR__)
