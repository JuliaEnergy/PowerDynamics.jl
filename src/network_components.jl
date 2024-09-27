####
#### Network level Bus representation
####
struct Bus{VF<:VertexFunction}
    compf::VF
end
function Bus(sys::ODESystem; verbose=false)
    if !isbusmodel(sys)
        msg = "The system must satisfy the bus model interface!"
        if iscomponentmodel(sys)
            msg *= " $(get_name(sys)) satisfies the component interface, did you mean to use `MTKBus`?"
        end
        throw(ArgumentError(msg))
    end
    io = _busio(sys, :busbar)
    vertexf = ODEVertex(sys, io.in, io.out; verbose)
    Bus(vertexf)
end

function simplify_mtkbus(sys; busbar=:busbar)
    @argcheck isbusmodel(sys) "The system must satisfy the bus model interface!"
    io = _busio(sys, busbar)
    structural_simplify(sys, (io.in, io.out))[1]
end

function _busio(sys, busbar)
    (;in=[getproperty(sys, busbar; namespace=false).i_r,
          getproperty(sys, busbar; namespace=false).i_i],
     out=[getproperty(sys, busbar; namespace=false).u_r,
          getproperty(sys, busbar; namespace=false).u_i])
end


####
#### Network level Line representation
####
struct Line{EF<:EdgeFunction}
    compf::EF
end
function Line(sys::ODESystem; verbose=false)
    if !islinemodel(sys)
        msg = "The system must satisfy the ine model interface!"
        if isbranchmodel(sys)
            msg *= " $(get_name(sys)) satisfies the branch interface, did you mean to use `MTKLine`?"
        end
        throw(ArgumentError(msg))
    end
    io = _lineio(sys, :src, :dst)
    edgef = StaticEdge(sys, io.srcin, io.dstin, io.out, Fiducial(); verbose)
    Line(edgef)
end

function simplify_mtkline(sys; src=:src, dst=:dst)
    @argcheck islinemodel(sys) "The system must satisfy the lie model interface!"
    io = _lineio(sys, src, dst)
    in = vcat(io.srcin, io.dstin)
    structural_simplify(sys, (in, io.out))[1]
end
function _lineio(sys, src, dst)
    (;srcin=[getproperty(sys, src; namespace=false).u_r,
             getproperty(sys, src; namespace=false).u_i],
     dstin=[getproperty(sys, dst; namespace=false).u_r,
            getproperty(sys, dst; namespace=false).u_i],
     out=[getproperty(sys, dst; namespace=false).i_r,
          getproperty(sys, dst; namespace=false).i_i,
          getproperty(sys, src; namespace=false).i_r,
          getproperty(sys, src; namespace=false).i_i])
end
