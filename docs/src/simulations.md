# Simulations

After having defined a [`PowerGrid`](@ref) a simulation can be run.

Currently PowerDynamics supports the following simulations:
- [`Perturbation`](@ref)
- [`LineFault`](@ref)
- [`NodeShortCircuit`](@ref)
- PowerPerturbation (power scale down/up of a generator node)

```@autodocs
Modules = [PowerDynamics]
Pages   = ["simulations.jl","NodeShortCircuit.jl"]
Order   = [:function, :type]
```
