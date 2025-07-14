using Pkg
Pkg.activate(@__DIR__)    # Ensures Project.toml/Manifest.toml are in example/env/
Pkg.develop(path=joinpath(@__DIR__, "..", ".."))  # Point to root of MarinePowerDynamics.jl
Pkg.add(PackageSpec(name="TimerOutputs", version="0.5.24")) # known dependency issue
Pkg.instantiate()
Pkg.precompile()