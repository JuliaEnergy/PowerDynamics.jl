# Line Types

The currently implemented line types are

```@eval
using InteractiveUtils, PowerDynamics, Markdown
linetypes = subtypes(PowerDynamics.AbstractLine)
join(["* [`$n`](@ref PowerDynamics.$n)" for n in linetypes], "\n") |> Markdown.parse
```

## Detailed Line Type Documentation

```@autodocs
Modules = [PowerDynamics]
Filter = t -> typeof(t) === DataType && t <: PowerDynamics.AbstractLine
```
