####
#### Network level Bus representation
####

"""
    compile_bus(sys::System; verbose=false, name=getname(sys), assume_io_coupling=false, check=true, kwargs...)

Create a VertexModel from an `System` that satisfies the bus model interface.

# Arguments
- `sys::System`: The system must satisfy the bus model interface (see [`isbusmodel`](@ref))
- `verbose::Bool=false`: Enable verbose output during creation
- `name`: Name for the bus (defaults to system name)
- `assume_io_coupling::Bool=false`: If true, assume output depends on inputs (see NetworkDynamics.jl docs)
- `check::Bool=true`: If false, skip component validation checks
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
function compile_bus(sys::System; verbose=false, name=getname(sys), assume_io_coupling=false, check=true, kwargs...)
    if !isbusmodel(sys)
        msg = "The system must satisfy the bus model interface!"
        if isinjectormodel(sys)
            msg *= " $(get_name(sys)) satisfies the component interface, did you mean to use `MTKBus`?"
        end
        throw(ArgumentError(msg))
    end
    io = _busio(sys, :busbar)
    vertexf = VertexModel(sys, io.in, io.out; verbose, name, assume_io_coupling, check)
    compile_bus(vertexf; copy=false, check, kwargs...)
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
    Bus(args...; kwargs...) => compile_bus(args...; kwargs...)

DEPRECATED! Use [`compile_bus`](@ref) instead.
"""
function Bus(args...; kwargs...)
    @warn "`Bus(...)` is deprecated in favor of `compile_bus(...)`" maxlog=1
    compile_bus(args...; kwargs...)
end

"""
    simplify_mtkbus(sys::System; busbar=:busbar)

Structurally simplify a bus model `System` by eliminating equations.

Closely matches what `VertexModel` does, but returns the `System` after
the simplifications rather than compiling it into a `VertexModel`.
"""
function simplify_mtkbus(sys::System; busbar=:busbar)
    @argcheck isbusmodel(sys) "The system must satisfy the bus model interface!"
    io = _busio(sys, busbar)
    mtkcompile(sys; inputs=io.in, outputs=io.out, simplify=false)
end

function _busio(sys::System, busbar)
    (;in=[getproperty(sys, busbar; namespace=false).i_r,
          getproperty(sys, busbar; namespace=false).i_i],
     out=[getproperty(sys, busbar; namespace=false).u_r,
          getproperty(sys, busbar; namespace=false).u_i])
end


####
#### Network level Line representation
####
"""
    compile_line(sys::System; verbose=false, name=getname(sys), assume_io_coupling=false, check=true, kwargs...)

Create an EdgeModel from a `System` that satisfies the line model interface.

# Arguments
- `sys::System`: The system must satisfy the line model interface (see [`islinemodel`](@ref))
- `verbose::Bool=false`: Enable verbose output during creation
- `name`: Name for the line (defaults to system name)
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
function compile_line(sys::System; verbose=false, name=getname(sys), assume_io_coupling=false, check=true, kwargs...)
    if !islinemodel(sys)
        msg = "The system must satisfy the line model interface!"
        if isbranchmodel(sys)
            msg *= " $(get_name(sys)) satisfies the branch interface, did you mean to use `MTKLine`?"
        end
        throw(ArgumentError(msg))
    end
    io = _lineio(sys, :src, :dst)
    edgef = EdgeModel(sys, io.srcin, io.dstin, io.srcout, io.dstout; verbose, name, assume_io_coupling, check)
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
    Line(args...; kwargs...) => compile_line(args...; kwargs...)

DEPRECATED! Use [`compile_line`](@ref) instead.
"""
function Line(args...; kwargs...)
    @warn "`Line(...)` is deprecated in favor of `compile_line(...)`" maxlog=1
    compile_line(args...; kwargs...)
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
