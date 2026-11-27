function isterminal(sys::System)
    vars = Set(getname.(unknowns(sys))) == Set([:u_r, :u_i, :i_r, :i_i])
    inputs = isempty(ModelingToolkitBase.unbound_inputs(sys))
    outputs = isempty(ModelingToolkitBase.unbound_outputs(sys))
    vars && inputs && outputs
end
"""
    isinjectormodel(sys::System)

Check if a `System` satisfies the injector model interface.

An injector model must contain a [`Terminal`](@ref) which connects it to the bus, by default
named `:terminal` (see [`get_injector_terminal`](@ref)). Injector models represent
components like generators, loads, and other devices that connect to a single bus. They can
have arbitrary internal complexity, including further terminals inside.

```asciiart
   (t)    ┌──────────┐
    o─────┤ Injector │
:terminal └──────────┘
```

A [`BusBar`](@ref) has a terminal too, but it is not an injector.

See also: [`Terminal`](@ref), [`injector_terminal`](@ref)
"""
function isinjectormodel(sys::System)
    isbusbar(sys) && return false
    name = get_injector_terminal(sys)
    subsystems = ModelingToolkitBase.get_systems(sys)
    idx = findfirst(s -> getname(s) == name, subsystems)
    isnothing(idx) ? false : isterminal(subsystems[idx])
end

# metadata key for the name of the injector terminal
struct InjectorTerminal end

"""
    get_injector_terminal(sys::System)

Name of the terminal which connects an injector to the bus: `:terminal` unless set
otherwise with [`set_injector_terminal`](@ref). Does not check whether `sys` is an injector.
"""
get_injector_terminal(sys::System) = ModelingToolkitBase.getmetadata(sys, InjectorTerminal, :terminal)

"""
    injector_terminal(sys::System)

The terminal which connects the injector `sys` to the bus, as namespaced subsystem, for
example for `connect(busbar.terminal, injector_terminal(inj))`. Throws if `sys` is no
injector model.

See also: [`get_injector_terminal`](@ref), [`isinjectormodel`](@ref)
"""
function injector_terminal(sys::System)
    isinjectormodel(sys) || throw(ArgumentError("$(getname(sys)) is not an injector model"))
    getproperty(sys, get_injector_terminal(sys))
end

"""
    set_injector_terminal(sys::System, name::Symbol)

Return a copy of `sys` whose terminal `name` connects it to the bus, for an injector with
several terminals. The copy has to satisfy the injector interface.

```julia
inj = set_injector_terminal(inj, :port)
injector_terminal(inj) # inj.port
```

See also: [`get_injector_terminal`](@ref), [`isinjectormodel`](@ref)
"""
function set_injector_terminal(sys::System, name::Symbol)
    newsys = ModelingToolkitBase.setmetadata(sys, InjectorTerminal, name)
    if !isinjectormodel(newsys)
        throw(ArgumentError("$(getname(sys)) with injector terminal $name does not satisfy the \
                             injector interface, it needs a Terminal subsystem $name."))
    end
    newsys
end


function isbusbar(sys::System)
    # TODO: consider overloading constructors in libary to add metadata
    # @set! busbar.metadata = (; modelname=:busbar)
    vars = getname.(unknowns(sys)) ⊇ Set([:u_r, :u_i, :i_r, :i_i])
    inputs = Set(getname.(ModelingToolkitBase.unbound_inputs(sys))) == Set([:i_r, :i_i])
    outputs = Set(getname.(ModelingToolkitBase.unbound_outputs(sys))) == Set([:u_r, :u_i])
    vars && inputs && outputs
end
"""
    isbusmodel(sys::System)

Check if a `System` satisfies the bus model interface.

A bus model must contain a component named `:busbar` that satisfies the busbar
interface. Bus models represent the complete dynamics of a power system bus and
can be transformed into a [`VertexModel`](@extref NetworkDynamics.VertexModel-Tuple{})
using [`compile_bus`](@ref).

```asciiart
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
    inputs = Set(getname.(ModelingToolkitBase.unbound_inputs(sys))) == Set([:u_r, :u_i])
    outputs = Set(getname.(ModelingToolkitBase.unbound_outputs(sys))) == Set([:i_r, :i_i])
    vars && inputs && outputs
end
"""
    islinemodel(sys::System)

Check if a `System` satisfies the line model interface.

A line model must contain two components named `:src` and `:dst` that both
satisfy the line end interface. Line models represent transmission lines and can
be transformed into an [`EdgeModel`](@extref NetworkDynamics.EdgeModel-Tuple{}) using [`compile_line`](@ref).

```asciiart
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

```asciiart
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
