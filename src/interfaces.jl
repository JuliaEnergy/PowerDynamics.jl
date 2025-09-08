function isterminal(sys::System)
    vars = Set(getname.(unknowns(sys))) == Set([:u_r, :u_i, :i_r, :i_i])
    inputs = isempty(ModelingToolkit.unbound_inputs(sys))
    outputs = isempty(ModelingToolkit.unbound_outputs(sys))
    vars && inputs && outputs
end
"""
    isinjectormodel(sys::System)

Check if a `System` satisfies the injector model interface.

An injector model must contain a [`Terminal`](@ref) named `:terminal`. Injector
models represent components like generators, loads, and other devices that
connect to a single bus. They can have arbitrary internal complexity as long as they
have exactly one terminal.

```
   (t)    ┌──────────┐
    o─────┤ Injector │
:terminal └──────────┘
```

See also: [`Terminal`](@ref)
"""
function isinjectormodel(sys::System)
    # HACK: use try catch instead of hasproperty https://github.com/SciML/ModelingToolkit.jl/issues/3016
    try
        return isterminal(sys.terminal)
    catch e
        if e isa ArgumentError && contains(e.msg, "variable terminal does not exist")
            return false
        else
            rethrow(e)
        end
    end
end


function isbusbar(sys::System)
    # TODO: consider overloading constructors in libary to add metadata
    # @set! busbar.metadata = (; modelname=:busbar)
    vars = getname.(unknowns(sys)) ⊇ Set([:u_r, :u_i, :i_r, :i_i])
    inputs = Set(getname.(ModelingToolkit.unbound_inputs(sys))) == Set([:i_r, :i_i])
    outputs = Set(getname.(ModelingToolkit.unbound_outputs(sys))) == Set([:u_r, :u_i])
    vars && inputs && outputs
end
"""
    isbusmodel(sys::System)

Check if a `System` satisfies the bus model interface.

A bus model must contain a component named `:busbar` that satisfies the busbar
interface. Bus models represent the complete dynamics of a power system bus and
can be transformed into a [`VertexModel`](@extref NetworkDynamics.VertexModel-Tuple{})
using [`compile_bus`](@ref).

```
┌───────────────────────────┐
│BusModel     ┌────────────┐│
│           ┌─┤ Injector 1 ││
│┌────────┐ │ └────────────┘│
││ BusBar ├─o               │
│└────────┘ │               │
│ :busbar   └ ...           │
│                           │
└───────────────────────────┘
```
Note: The BusModel musst contain exaclty one `BusBar`, the rest of the structure is free.
For example, you could also put a Brach between an injector and a Busbar or have multiple
injectors and controllers connected.

See also: [`compile_bus`](@ref), [`BusBar`](@ref), [`MTKBus`](@ref)
"""
function isbusmodel(sys::System)
    # HACK: use try catch instead of hasproperty https://github.com/SciML/ModelingToolkit.jl/issues/3016
    try
        return isbusbar(sys.busbar)
    catch e
        if e isa ArgumentError && contains(e.msg, "variable busbar does not exist")
            return false
        else
            rethrow(e)
        end
    end
end


function islineend(sys::System)
    vars = getname.(unknowns(sys)) ⊇ Set([:u_r, :u_i, :i_r, :i_i])
    inputs = Set(getname.(ModelingToolkit.unbound_inputs(sys))) == Set([:u_r, :u_i])
    outputs = Set(getname.(ModelingToolkit.unbound_outputs(sys))) == Set([:i_r, :i_i])
    vars && inputs && outputs
end
"""
    islinemodel(sys::System)

Check if a `System` satisfies the line model interface.

A line model must contain two components named `:src` and `:dst` that both
satisfy the line end interface. Line models represent transmission lines and can
be transformed into an [`EdgeModel`](@extref NetworkDynamics.EdgeModel-Tuple{}) using [`compile_line`](@ref).

```
┌──────────────────────────────────────┐
│LineModel     ┌────────┐              │
│            ┌─┤ Branch ├─┐            │
│┌─────────┐ │ └────────┘ │ ┌─────────┐│
││ LineEnd ├─o            o─┤ LineEnd ││
│└─────────┘ │            │ └─────────┘│
│   :src     └    ....    ┘    :dst    │
│                                      │
└──────────────────────────────────────┘
```
Note: Between the `LineEnd`s there can be arbeitrary structures, for example branches in
series or parallel.

See also: [`compile_line`](@ref), [`LineEnd`](@ref), [`MTKBus`](@ref)
"""
function islinemodel(sys::System)
    # HACK: use try catch instead of hasproperty https://github.com/SciML/ModelingToolkit.jl/issues/3016
    try
        return islineend(sys.src)
    catch e
        if e isa ArgumentError && contains(e.msg, "variable src does not exist")
            return false
        else
            rethrow(e)
        end
    end
    try
        return islineend(sys.dst)
    catch e
        if e isa ArgumentError && contains(e.msg, "variable dst does not exist")
            return false
        else
            rethrow(e)
        end
    end
end

"""
    isbranchmodel(sys::System)

Check if a `System` satisfies the branch model interface.

A branch model must contain two [`Terminal`](@ref) components named `:src` and
`:dst`. Branch models represent two-port network elements like transmission
lines, transformers, and other connecting devices.

```
 (t) ┌────────┐ (t)
  o──┤ Branch ├──o
:src └────────┘ :dst
```

See also: [`Terminal`](@ref)
"""
function isbranchmodel(sys::System)
    # HACK: use try catch instead of hasproperty https://github.com/SciML/ModelingToolkit.jl/issues/3016
    try
        return isterminal(sys.src)
    catch e
        if e isa ArgumentError && contains(e.msg, "variable src does not exist")
            return false
        else
            rethrow(e)
        end
    end
    try
        return isterminal(sys.dst)
    catch e
        if e isa ArgumentError && contains(e.msg, "variable dst does not exist")
            return false
        else
            rethrow(e)
        end
    end
end
