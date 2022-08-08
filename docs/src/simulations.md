# Simulations or Fault Scenarios

After having defined a [`PowerGrid`](@ref) a simulation can be run.

Currently PowerDynamics supports the following fault/simulation types:

```@eval
using InteractiveUtils, PowerDynamics, Markdown
perturbations = [ChangeInitialConditions; subtypes(PowerDynamics.AbstractPerturbation)]
join(["* [`$n`](@ref PowerDynamics.$n)" for n in perturbations], "\n") |> Markdown.parse
```

A specific scenario can be simulated for a given powergrid using the `simulate` function
```@docs
simulate
```


## Detailed Simulation/ Fault Type Documentation
### Change Initial Conditions
```@docs
ChangeInitialConditions
```

### `AbstractPertubation` types
```@autodocs
Modules = [PowerDynamics]
Filter = t -> typeof(t) === DataType && t<:PowerDynamics.AbstractPerturbation
```

### LineFault
Deprecated! Use [`LineFailure`](@ref) instead!
```@docs
LineFault
```
