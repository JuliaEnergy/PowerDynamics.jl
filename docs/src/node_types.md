# Node Types

The currently implementes node types are

```@eval
using InteractiveUtils, PowerDynamics, Markdown
nodetypes = subtypes(PowerDynamics.AbstractNodeParameters)
join(["* [`$n`](@ref PowerDynamics.$n)" for n in nodetypes], "\n") |> Markdown.parse
```

## Detailed Node Type Documentation

```@docs
PowerDynamics.AbstractNodeParameters
```

```@autodocs
Modules = [PowerDynamics]
Filter = t -> typeof(t) === DataType
```
