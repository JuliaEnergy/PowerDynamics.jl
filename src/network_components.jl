####
#### Network level Bus representation
####

"""
    compile_bus(sys::System; verbose=false, name=nothing, assume_io_coupling=false, check=true, kwargs...)

Create a VertexModel from an `System` that satisfies the bus model interface.

As a convenience, a system which satisfies the *injector* model interface (see
[`isinjectormodel`](@ref)) instead is wrapped into a single-injector bus automatically:
`compile_bus(machine)` is fully equivalent to `compile_bus(MTKBus(machine))`, including the
resulting name (`:bus`, the [`MTKBus`](@ref) default) — pass `name` to override it.

# Arguments
- `sys::System`: The system must satisfy the bus model interface (see [`isbusmodel`](@ref))
  or the injector model interface, in which case it is wrapped in an [`MTKBus`](@ref)
- `verbose::Bool=false`: Enable verbose output during creation
- `name`: Name for the bus (defaults to the name of the bus system, i.e. `:bus` when a bare
  injector was wrapped)
- `assume_io_coupling::Bool=false`: If true, assume output depends on inputs (see NetworkDynamics.jl docs)
- `check::Bool=true`: If false, skip component validation checks
- `current_source::Bool=false`: If true, compile for use as an injector node via a [`LoopbackConnection`](@extref NetworkDynamics.LoopbackConnection)
  (i.e. switches from current in voltage out to voltage in current out)
- `kwargs...`: Additional keyword arguments passed to the Bus constructor

# Returns
- A [`VertexModel`](@extref NetworkDynamics.VertexModel-Tuple{}) representing the bus

```asciiart
                                                  ╔═════════════════════════╗
                                                  ║ VertexModel (compiled)  ║
            ┌────────────────────┐      Network   ║  ┌────────────────────┐ ║
            │MTKBus   ┌─────────┐│     interface  ║  │MTKBus   ┌─────────┐│ ║
            │        ┌┤Generator││                ║  │        ┌┤Generator││ ║
            │┌──────┐│└─────────┘│      current ────→│┌──────┐│└─────────┘│ ║
compile_bus(││BusBar├o           │) =>            ║  ││BusBar├o           │ ║
            │└──────┘│┌────┐     │      voltage ←────│└──────┘│┌────┐     │ ║
            │        └┤Load│     │                ║  │        └┤Load│     │ ║
            │         └────┘     │                ║  │         └────┘     │ ║
            └────────────────────┘                ║  └────────────────────┘ ║
                                                  ╚═════════════════════════╝
```

See also: [`MTKBus`](@ref)
"""
function compile_bus(
    sys::System;
    verbose=false,
    name=nothing,
    assume_io_coupling=false,
    check=true,
    current_source=false,
    ff_to_constraint=!current_source,
    extin=nothing,
    mtkcompile=nothing,
    Vbase=nothing,
    kwargs...)
    if !isbusmodel(sys)
        if isinjectormodel(sys)
            sys = MTKBus(sys)
        else
            throw(ArgumentError("The system must satisfy the bus model interface!"))
        end
    end
    # resolved *after* wrapping, so the sugar is equivalent to the explicit spelling
    name = isnothing(name) ? getname(sys) : name
    # `Vbase` convenience kwarg: sugar for the `busbar₊Vbase => …` default override
    !isnothing(Vbase) && haskey(kwargs, Symbol("busbar₊Vbase")) &&
        throw(ArgumentError("Pass either `Vbase` or `busbar₊Vbase`, not both."))
    vbase_override = isnothing(Vbase) ? (;) : NamedTuple{(Symbol("busbar₊Vbase"),)}((Vbase,))
    sys = add_systembase(sys) # no-opt if present
    io = _busio(sys, :busbar, current_source)
    vertexf = VertexModel(sys, io.in, io.out; verbose, name, assume_io_coupling, check, ff_to_constraint, extin, mtkcompile)
    _resolve_busbar_vbase!(vertexf, current_source)
    compile_bus(vertexf; copy=false, check, kwargs..., vbase_override...)
end

# The busbar's `Vbase` fallback is declared as a *weak init formula* (see `BusBase`) so that
# its provenance stays distinguishable from a value the user set, which always arrives as a
# real `default`. `compile_bus` is the one place that knows whether a bus is a *satellite*
# (an injector node attached to a hub bus via a `LoopbackConnection`), so it is where that
# fallback is resolved into its final storage class:
#
#   satellite  → drop the fallback and inherit from the hub instead. The swap is
#                unconditional: a user-set `Vbase` is a `default` and beats the inherited
#                weak formula on its own, so no `has_default` test is needed here.
#   plain bus  → *materialize* the fallback into a real `default`. This is not just tidiness:
#                `default_from` (used by line ends and by satellites) copies the source's
#                `default` in a pre-pass that runs before any init formula, so it cannot see
#                a weak formula. A bus that kept its fallback as a formula would silently
#                fail to feed its neighbours — `default_from` would skip and the line end
#                would be left with no value at all.
#
# Either way the fallback formula is consumed here, so it never reaches init.
function _resolve_busbar_vbase!(vertexf, satellite)
    sym = :busbar₊Vbase
    # get rid of weak init formula set by BusBar(), capture its value
    fallback = _strip_weak_initf!(vertexf, sym)
    if satellite
        set_default_from!(vertexf, sym, (:hub, sym))
    elseif !has_default(vertexf, sym)
        if !isnothing(fallback)
            set_default!(vertexf, sym, fallback)
        elseif sym ∈ psym(vertexf)
            # `Vbase` exists but has neither a default nor the fallback formula
            throw(ArgumentError("Could not resolve `$sym` of bus `$(vertexf.name)`: it has no \
                default and no fallback to materialize. Pass `compile_bus(…; Vbase=…)` to set \
                it explicitly."))
        end
        # `sym` absent entirely: a hand-rolled, purely per-unit bus without bases. skip.
    end
    vertexf
end
function _strip_weak_initf!(c, sym)
    has_initformula(c) || return nothing
    formulas = collect(get_initformulas(c))
    idx = findfirst(f -> f.weak && f.outsym == [sym] && isempty(f.sym), formulas)
    isnothing(idx) && return nothing
    f = popat!(formulas, idx)
    isempty(formulas) ? delete_initformulas!(c) : set_initformula!(c, Tuple(formulas))
    out = SymbolicView(zeros(1), f.outsym)
    f(out, Float64[])
    out[sym]
end

"""
    compile_bus(template::VertexModel; copy=true, pf=nothing, name=template.name, check=true, pairs...)

Similar to the `Bus` constructor, but takes a pre-compiled `VertexModel`. It copies the VertexModel
and applies the keyword arguments. This is useful when you want to create new bus models based on a template.
"""
function compile_bus(template::VertexModel; copy=true, vidx=nothing, pf=nothing, name=template.name, check=true, pairs...)
    vertexf = copy ? Base.copy(template) : template
    if name != template.name
        vertexf = VertexModel(vertexf; name, allow_output_sym_clash=true, check)
    end

    # is done in ND constructor too, but needs special handling because compile_line calls this
    if !isnothing(vidx)
        set_graphelement!(vertexf, vidx)
    end
    if !isnothing(pf)
        if ispfmodel(pf)
            set_pfmodel!(vertexf, pf)
        else
            throw(ArgumentError("The provided pf model is invalid!"))
        end
    end
    for (v, d) in pairs
        set_default!(vertexf, v, d)
    end

    vertexf
end

"""
    simplify_mtkbus(sys::System; busbar=:busbar)

Structurally simplify a bus model `System` by eliminating equations.

Closely matches what `VertexModel` does, but returns the `System` after
the simplifications rather than compiling it into a `VertexModel`.
"""
function simplify_mtkbus(sys::System; busbar=:busbar, current_source=false)
    @argcheck isbusmodel(sys) "The system must satisfy the bus model interface!"
    io = _busio(sys, busbar, current_source)
    mtkcompile(sys; inputs=io.in, outputs=io.out, simplify=false)
end

function _busio(sys::System, busbar, current_source)
    i = [getproperty(sys, busbar; namespace=false).i_r,
         getproperty(sys, busbar; namespace=false).i_i]
    u = [getproperty(sys, busbar; namespace=false).u_r,
         getproperty(sys, busbar; namespace=false).u_i]
    if current_source
        (; in=u, out=i)
    else
        (; in=i, out=u)
    end
end


####
#### Network level Line representation
####
"""
    compile_line(sys::System; verbose=false, name=nothing, assume_io_coupling=false, check=true, kwargs...)

Create an EdgeModel from a `System` that satisfies the line model interface.

As a convenience, a system which satisfies the *branch* model interface (see
[`isbranchmodel`](@ref)) instead is wrapped into a single-branch line automatically:
`compile_line(piline)` is fully equivalent to `compile_line(MTKLine(piline))`, including the
resulting name (`:line`, the [`MTKLine`](@ref) default) — pass `name` to override it.

# Arguments
- `sys::System`: The system must satisfy the line model interface (see [`islinemodel`](@ref))
  or the branch model interface, in which case it is wrapped in an [`MTKLine`](@ref)
- `verbose::Bool=false`: Enable verbose output during creation
- `name`: Name for the line (defaults to the name of the line system, i.e. `:line` when a
  bare branch was wrapped)
- `assume_io_coupling::Bool=false`: If true, assume output depends on inputs (see NetworkDynamics.jl docs)
- `check::Bool=true`: If false, skip component validation checks
- `kwargs...`: Additional keyword arguments passed to the Line constructor

# Returns
- An [`EdgeModel`](@extref NetworkDynamics.EdgeModel-Tuple{}) representing the line


```asciiart

                                                     ╔══════════════════════════════╗
                                                     ║ EdgeModel (compiled)         ║
             ┌─────────────────────────────┐     src ║ ┌──────────────────────────┐ ║ dst
             │MTKLine   ┌───────┐          │  vertex ║ │MTKLine   ┌────┐          │ ║ vertex
             │         ┌┤BranchA├┐         │         ║ │         ┌┤    ├┐         │ ║
             │┌───────┐│└───────┘│┌───────┐│     u ───→│┌───────┐│└────┘│┌───────┐│←─── u
compile_line(││LineEnd├o         o┤LineEnd││) =>     ║ ││LineEnd├o      o┤LineEnd││ ║
             │└───────┘│┌───────┐│└───────┘│     i ←───│└───────┘│┌────┐│└───────┘│───→ i
             │  :src   └┤BranchB├┘  :dst   │         ║ │         └┤    ├┘         │ ║
             │          └───────┘          │         ║ │          └────┘          │ ║
             └─────────────────────────────┘         ║ └──────────────────────────┘ ║
                                                     ╚══════════════════════════════╝
```


See also: [`MTKLine`](@ref)
"""
function compile_line(
    sys::System;
    verbose=false,
    name=nothing,
    assume_io_coupling=false,
    check=true,
    mtkcompile=nothing,
    kwargs...
)
    if !islinemodel(sys)
        if isbranchmodel(sys)
            sys = MTKLine(sys)
        else
            throw(ArgumentError("The system must satisfy the line model interface!"))
        end
    end
    # resolved *after* wrapping, so the sugar is equivalent to the explicit spelling
    name = isnothing(name) ? getname(sys) : name
    sys = add_systembase(sys) # no-opt if present
    io = _lineio(sys, :src, :dst)
    edgef = EdgeModel(sys, io.srcin, io.dstin, io.srcout, io.dstout; verbose, name, assume_io_coupling, check, mtkcompile)
    compile_line(edgef; copy=false, check, kwargs...)
end
function compile_line(edgef::EdgeModel; copy=true, src=nothing, dst=nothing, name=edgef.name, check=true, pairs...)
    if copy
        edgef = Base.copy(edgef)
    end
    if name != edgef.name
        edgef = EdgeModel(edgef; name, allow_output_sym_clash=true, check)
    end

    # is done in ND constructor too, but needs special handling because compile_line calls this
    if !isnothing(src) && !isnothing(dst)
        set_graphelement!(edgef, (;src, dst))
    end
    for (v, d) in pairs
        set_default!(edgef, v, d)
    end

    edgef
end

"""
    simplify_mtkline(sys::System; src=:src, dst=:dst)

Structurally simplify a line model `System` by eliminating equations.

Closely matches what `EdgeModel` does, but returns the `System` after
the simplifications rather than compiling it into an `EdgeModel`.
"""
function simplify_mtkline(sys::System; src=:src, dst=:dst)
    @argcheck islinemodel(sys) "The system must satisfy the line model interface!"
    io = _lineio(sys, src, dst)
    in = vcat(io.srcin, io.dstin)
    out = vcat(io.srcout, io.dstout)
    mtkcompile(sys; inputs=in, outputs=out, simplify=false)
end
function _lineio(sys::System, src, dst)
    (;srcin=[getproperty(sys, src; namespace=false).u_r,
             getproperty(sys, src; namespace=false).u_i],
     dstin=[getproperty(sys, dst; namespace=false).u_r,
            getproperty(sys, dst; namespace=false).u_i],
     srcout=[getproperty(sys, src; namespace=false).i_r,
             getproperty(sys, src; namespace=false).i_i],
     dstout=[getproperty(sys, dst; namespace=false).i_r,
             getproperty(sys, dst; namespace=false).i_i,])
end
