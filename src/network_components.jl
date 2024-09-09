####
#### Network level Bus representation
####
struct Bus{VF<:VertexFunction}
    vertexf::VF
end
function Bus(sys::ODESystem)
    if !isbusmodel(sys)
        msg = "The system must satisfy the bus model interface!"
        if iscomponentmodel(sys)
            msg *= " $(sys.name) satisfies the component interface, did you mean to use `BusModel`?"
        end
        throw(ArgumentError("The system must satisfy the bus model interface!"))
    end
    io = _busio(sys, busbar)
    vertexf = ODEVertex(sys, io.in, io.out)
end

function simplify_busmodel(sys; busbar=:busbar)
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
    edgef::EF
end
function Line(sys::ODESystem)
    if !islinemodel(sys)
        msg = "The system must satisfy the ine model interface!"
        if isbranchmodel(sys)
            msg *= " $(sys.name) satisfies the branch interface, did you mean to use `LineModel`?"
        end
        throw(ArgumentError("The system must satisfy the bus model interface!"))
    end
    io = _lineio(sys, src, dst)
    StaticEdge(sys, io.srcin, io.dstin, io.out, Fiducial())
end

function simplify_linemodel(sys; src=:src, dst=:dst)
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
