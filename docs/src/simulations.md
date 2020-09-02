# Simulations or Fault Scenarios

After having defined a [`PowerGrid`](@ref) a simulation can be run.

Currently PowerDynamics supports the following simulations:

- [`ChangeInitialConditions`](@ref)
```@eval
using InteractiveUtils, PowerDynamics, Markdown
perturbations = subtypes(PowerDynamics.AbstractPerturbation)
join(["* [`$n`](@ref PowerDynamics.$n)" for n in perturbations], "\n") |> Markdown.parse
```


## Detailed Simulation/ Fault Type Documentation

```@autodocs
Modules = [PowerDynamics]
Pages   = ["ChangeInitialConditions.jl","PowerPerturbation.jl","NodeParameterChange","LineFailure.jl","NodeShortCircuit.jl","AbstractPerturbation.jl"]
%Filter = t -> typeof(t) === DataType && (t <: PowerDynamics.AbstractPerturbation || %t<:PowerDynamics.ChangeInitialConditions)
Order   = [:function, :type]
```