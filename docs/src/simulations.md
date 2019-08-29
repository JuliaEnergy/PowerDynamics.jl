# Simulations

After having defined a [`PowerGrid`](@ref) a simulation can be run.

Currently PowerDynamics supports the following simulations:
- [`Perturbation`](@ref)
- [`LineFault`](@ref)
- [`NodeShortCircuit`](@ref)
- [`PowerPerturbation`](@ref)

```@autodocs
Modules = [PowerDynamics]
Pages   = ["simulations.jl","NodeShortCircuit.jl","PowerPerturbation.jl"]
Order   = [:function, :type]
```
