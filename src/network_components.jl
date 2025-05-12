####
#### Network level Bus representation
####
function Bus(sys::ODESystem; verbose=false, name=getname(sys), kwargs...)
    if !isbusmodel(sys)
        msg = "The system must satisfy the bus model interface!"
        if iscomponentmodel(sys)
            msg *= " $(get_name(sys)) satisfies the component interface, did you mean to use `MTKBus`?"
        end
        throw(ArgumentError(msg))
    end
    io = _busio(sys, :busbar)
    vertexf = VertexModel(sys, io.in, io.out; verbose, name)
    Bus(vertexf; copy=false, kwargs...)
end
function Bus(vertexf::VertexModel; copy=true, vidx=nothing, pf=nothing, name=vertexf.name, pairs...)
    if copy
        vertexf = Base.copy(vertexf)
    end
    if name != vertexf.name
        vertexf = VertexModel(vertexf; name, allow_output_sym_clash=true)
    end

    if !isnothing(vidx)
        set_graphelement!(vertexf, vidx)
    end
    if !isnothing(pf)
        if ispfmodel(pf)
            set_metadata!(vertexf, :pfmodel, pf)
        else
            throw(ArgumentError("The provided pf model is invalid!"))
        end
    end
    for (v, d) in pairs
        set_default!(vertexf, v, d)
    end
    vertexf
end

function simplify_mtkbus(sys::ODESystem; busbar=:busbar)
    @argcheck isbusmodel(sys) "The system must satisfy the bus model interface!"
    io = _busio(sys, busbar)
    structural_simplify(sys, (io.in, io.out))[1]
end

function _busio(sys::ODESystem, busbar)
    (;in=[getproperty(sys, busbar; namespace=false).i_r,
          getproperty(sys, busbar; namespace=false).i_i],
     out=[getproperty(sys, busbar; namespace=false).u_r,
          getproperty(sys, busbar; namespace=false).u_i])
end


####
#### Network level Line representation
####
function Line(sys::ODESystem; verbose=false, name=getname(sys), kwargs...)
    if !islinemodel(sys)
        msg = "The system must satisfy the ine model interface!"
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

function simplify_mtkline(sys::ODESystem; src=:src, dst=:dst)
    @argcheck islinemodel(sys) "The system must satisfy the lie model interface!"
    io = _lineio(sys, src, dst)
    in = vcat(io.srcin, io.dstin)
    out = vcat(io.srcout, io.dstout)
    structural_simplify(sys, (in, out))[1]
end
function _lineio(sys::ODESystem, src, dst)
    (;srcin=[getproperty(sys, src; namespace=false).u_r,
             getproperty(sys, src; namespace=false).u_i],
     dstin=[getproperty(sys, dst; namespace=false).u_r,
            getproperty(sys, dst; namespace=false).u_i],
     srcout=[getproperty(sys, src; namespace=false).i_r,
             getproperty(sys, src; namespace=false).i_i],
     dstout=[getproperty(sys, dst; namespace=false).i_r,
             getproperty(sys, dst; namespace=false).i_i,])
end
