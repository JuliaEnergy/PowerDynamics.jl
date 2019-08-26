# Node Types

The currently implemented node types are

```@eval
using InteractiveUtils, PowerDynamics, Markdown
nodetypes = subtypes(PowerDynamics.AbstractNode)
join(["* [`$n`](@ref PowerDynamics.$n)" for n in nodetypes], "\n") |> Markdown.parse
```

## Detailed Node Type Documentation

```@autodocs
Modules = [PowerDynamics]
Filter = t -> typeof(t) === DataType && t <: PowerDynamics.AbstractNode
```
