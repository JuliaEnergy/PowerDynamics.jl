# Node Types

The currently implementes node types are

```@eval
using InteractiveUtils, PowerDynBase, Markdown
nodetypes = subtypes(PowerDynBase.AbstractNodeParameters)
join(["* [`$n`](@ref PowerDynBase.$n)" for n in nodetypes], "\n") |> Markdown.parse
```

They are all subtypes of [`PowerDynBase.AbstractNodeParameters`](@ref).

## Detailed Node Type Documentation

```@docs
PowerDynBase.AbstractNodeParameters
```

```@autodocs
Modules = [PowerDynBase]
Filter = t -> typeof(t) === DataType && t !== PowerDynBase.AbstractNodeParameters && t <: PowerDynBase.AbstractNodeParameters
```
