####
#### Network level Bus representation
####

"""
    Bus(sys::System; verbose=false, name=getname(sys), kwargs...)

Create a VertexModel from an `System` that satisfies the bus model interface.

# Arguments
- `sys::System`: The system must satisfy the bus model interface (see [`isbusmodel`](@ref))
- `verbose::Bool=false`: Enable verbose output during creation
- `name`: Name for the bus (defaults to system name)
- `kwargs...`: Additional keyword arguments passed to the Bus constructor

# Returns
- A [`VertexModel`](@extref NetworkDynamics.VertexModel-Tuple{}) representing the bus

```

                                          ╔═════════════════════════╗
                                          ║ VertexModel (compiled)  ║
    ┌────────────────────┐      Network   ║  ┌────────────────────┐ ║
    │MTKBus   ┌─────────┐│     interface  ║  │MTKBus   ┌─────────┐│ ║
    │        ┌┤Generator││                ║  │        ┌┤Generator││ ║
    │┌──────┐│└─────────┘│      current ────→│┌──────┐│└─────────┘│ ║
Bus(││BusBar├o           │) =>            ║  ││BusBar├o           │ ║
    │└──────┘│┌────┐     │      voltage ←────│└──────┘│┌────┐     │ ║
    │        └┤Load│     │                ║  │        └┤Load│     │ ║
    │         └────┘     │                ║  │         └────┘     │ ║
    └────────────────────┘                ║  └────────────────────┘ ║
                                          ╚═════════════════════════╝
```

See also: [`MTKBus`](@ref)
"""
function Bus(sys::System; verbose=false, name=getname(sys), kwargs...)
    if !isbusmodel(sys)
        msg = "The system must satisfy the bus model interface!"
        if isinjectormodel(sys)
            msg *= " $(get_name(sys)) satisfies the component interface, did you mean to use `MTKBus`?"
        end
        throw(ArgumentError(msg))
    end
    io = _busio(sys, :busbar)
    vertexf = VertexModel(sys, io.in, io.out; verbose, name)
    Bus(vertexf; copy=false, kwargs...)
end
"""
    Bus(template::VertexModel; copy=true, vidx=nothing, pf=nothing, name=template.name, pairs...)

Similar to the `Bus` constructor, but takes a pre-compiled `VertexModel`. It copies the VertexModel
and applies the keyword arguments. This is usefull when you want to create new bus models based on a template.
"""
function Bus(template::VertexModel; copy=true, vidx=nothing, pf=nothing, name=template.name, pairs...)
    vertexf = copy ? Base.copy(template) : template
    if name != template.name
        vertexf = VertexModel(vertexf; name, allow_output_sym_clash=true)
    end

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
    Line(sys::System; verbose=false, name=getname(sys), kwargs...)

Create an EdgeModel from a `System` that satisfies the line model interface.

# Arguments
- `sys::System`: The system must satisfy the line model interface (see [`islinemodel`](@ref))
- `verbose::Bool=false`: Enable verbose output during creation
- `name`: Name for the line (defaults to system name)
- `kwargs...`: Additional keyword arguments passed to the Line constructor

# Returns
- An [`EdgeModel`](@extref NetworkDynamics.EdgeModel-Tuple{}) representing the line


```

                                             ╔══════════════════════════════╗
                                             ║ EdgeModel (compiled)         ║
     ┌─────────────────────────────┐     src ║ ┌──────────────────────────┐ ║ dst
     │MTKLine   ┌───────┐          │  vertex ║ │MTKLine   ┌────┐          │ ║ vertex
     │         ┌┤BranchA├┐         │         ║ │         ┌┤    ├┐         │ ║
     │┌───────┐│└───────┘│┌───────┐│     u ───→│┌───────┐│└────┘│┌───────┐│←─── u
Line(││LineEnd├o         o┤LineEnd││) =>     ║ ││LineEnd├o      o┤LineEnd││ ║
     │└───────┘│┌───────┐│└───────┘│     i ←───│└───────┘│┌────┐│└───────┘│───→ i
     │  :src   └┤BranchB├┘  :dst   │         ║ │         └┤    ├┘         │ ║
     │          └───────┘          │         ║ │          └────┘          │ ║
     └─────────────────────────────┘         ║ └──────────────────────────┘ ║
                                             ╚══════════════════════════════╝
```


See also: [`MTKLine`](@ref)
"""
function Line(sys::System; verbose=false, name=getname(sys), kwargs...)
    if !islinemodel(sys)
        msg = "The system must satisfy the line model interface!"
        if isbranchmodel(sys)
            msg *= " $(get_name(sys)) satisfies the branch interface, did you mean to use `MTKLine`?"
        end
        throw(ArgumentError(msg))
    end
    io = _lineio(sys, :src, :dst)
    edgef = EdgeModel(sys, io.srcin, io.dstin, io.srcout, io.dstout; verbose, name)
    Line(edgef; copy=false, kwargs...)
end
function Line(edgef::EdgeModel; copy=true, src=nothing, dst=nothing, name=edgef.name, pairs...)
    if copy
        edgef = Base.copy(edgef)
    end
    if name != edgef.name
        edgef = EdgeModel(edgef; name, allow_output_sym_clash=true)
    end

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
