
@doc doc"""
```julia
CompositeNode(CurrentNodes)
```

A composite node consisting of current sources provided. This assumes that the
current nodes are implemented as

```julia
du = i - I_node
```

with MassMatrix.

In principle we could build a composite node that incorporates power sources,
current sources, and one voltage regulator.
"""

struct _IPU end
struct _IU end
struct _PU end
struct _IP end
struct _I end
struct _P end


CompositionTypes = Dict(["IPU" => _IPU,
    "IU" => _IU,
    "PU" => _PU,
    "I" => _I,
    "P" => _P,
    "IP" => _IP])

@__doc__ mutable struct CompositeNode{T}
    CurrentNodes # constant current nodes
    VoltageNode # constant voltage nodes
    PowerNodes # constant power nodes
    internaldims::Array{Int, 1}
    total_dim::Int
    idxs::Array{Array{Int,1},1}
    symbols::Array{Symbol, 1}
    function CompositeNode(;CurrentNodes=nothing, VoltageNode=nothing, PowerNodes=nothing)# PowerNodes, VoltageNode)

        @assert any([CurrentNodes != nothing, PowerNodes != nothing])

        @assert !(typeof(VoltageNode) <: AbstractArray)

        symbols = [:u_r, :u_i]
        internaldims = Array{Int,1}()
        mix = ""

        if CurrentNodes != nothing
            mix *= "I"
            for current_node in CurrentNodes
                append!(symbols, symbolsof(current_node)[3:end])
                append!(internaldims, dimension(current_node) - 2)
            end
        end


        if PowerNodes != nothing
            mix *= "P"
            for power_node in PowerNodes
                append!(symbols, symbolsof(power_node)[3:end])
                append!(internaldims, dimension(power_node) - 2)
            end
        end


        if VoltageNode != nothing
            mix *= "U"
            append!(symbols, symbolsof(VoltageNode)[3:end])
            append!(internaldims, dimension(VoltageNode) - 2)
        end

        idxs = Array{Array{Int,1},1}([])
        offset = 2
        for idim in internaldims
            println([[[1,2]; (offset + 1):(offset + idim)]])
            append!(idxs, [[[1,2]; (offset + 1):(offset + idim)]])
            offset += idim
        end

        total_dim = 2 + sum(internaldims)
        CT = CompositionTypes[mix]

        new{CT}(CurrentNodes,
        VoltageNode,
        PowerNodes,
        internaldims,
        total_dim,
        idxs,
        symbols)
    end
end


function construct_vertex(cn::CompositeNode{T}) where T <: Union{_IPU, _IU, _PU}

    # first collect the information
    current_vertices = [construct_vertex(n) for n in cn.CurrentNodes]
    power_vertices = [construct_vertex(n) for n in cn.PowerNodes]
    voltage_vertex = construct_vertex(cn.VoltageNode)

    all_vertices = vcat(current_nodes, power_nodes, voltage_vertex)

    cfs = [n.f! for n in current_vertices]
    pfs = [n.f! for n in power_vertices]

    mass_matrix = vcat(voltage_vertex.mass_matrix[1:2], [cn.mass_matrix[3:end] for cn in all_vertices]...)

    num_cv = length(current_vertices)
    num_pv = length(power_vertices)

    function rhs!(dx, x, e_s, e_d, p, t)
        i = total_current(e_s, e_d)
        u = x[1] + x[2] * im
        # p = u*conj(i)

        i_injected = 0.0 + 0.0im
        p_injected = 0.0 + 0.0im

        @views begin

            # constant current
            for (cf, idx) in zip(cfs, cn.idxs[1:num_cv])
                cf(dx[idx], x[idx], e_s, e_d, p, t)
                # The current nodes are assumed to have du = i - I_node
                i_injected += i - complex(dx[1], dx[2])
            end

            for (pf, idx) in zip(pfs, cn.idxs[num_cv+1:num_cv+num_pv])
                pf(dx[idx], x[idx], e_s, e_d, p, t)
                # The current nodes are assumed to have du = p - P_node
                p_injected += p - complex(dx[1], dx[2])
            end

            total_i = i_injected + p_injected / u

            voltage_vertex.f!(dx[cn.idxs[end]], x[cn.idxs[end]], [e_s;[real(total_i), imag(total_i), 0., 0.]], e_d, p, t)

        end # views

    end

    ODEVertex(f! = rhs!, dim=cn.total_dim, mass_matrix=mass_matrix, sym=cn.symbols)
end



function construct_vertex(cn::CompositeNode{T}) where T <: Union{_IP, _I, _P}

    # first collect the information
    current_vertices = [construct_vertex(n) for n in cn.CurrentNodes]
    power_vertices = [construct_vertex(n) for n in cn.PowerNodes]

    all_vertices = vcat(current_nodes, power_nodes)

    cfs = [n.f! for n in current_vertices]
    pfs = [n.f! for n in power_vertices]

    mass_matrix = vcat([0,0], [cn.mass_matrix[3:end] for cn in all_vertices]...)

    num_cv = length(current_vertices)
    num_pv = length(power_vertices)

    function rhs!(dx, x, e_s, e_d, p, t)
        i = total_current(e_s, e_d)
        u = x[1] + x[2] * im
        # p = u*conj(i)

        i_injected = 0.0 + 0.0im
        p_injected = 0.0 + 0.0im

        @views begin

            # constant current
            for (cf, idx) in zip(cfs, cn.idxs[1:num_cv])
                cf(dx[idx], x[idx], e_s, e_d, p, t)
                # The current nodes are assumed to have du = i - I_node
                i_injected += i - complex(dx[1], dx[2])
            end

            for (pf, idx) in zip(pfs, cn.idxs[num_cv+1:num_cv+num_pv])
                pf(dx[idx], x[idx], e_s, e_d, p, t)
                # The current nodes are assumed to have du = p - P_node
                p_injected += p - complex(dx[1], dx[2])
            end

            total_p = u * conj(i_injected - i) + p_injected # Check signs

            dx[1] = real(total_p)
            dx[2] = imag(total_p)

        end # views

    end

    ODEVertex(f! = rhs!, dim=cn.total_dim, mass_matrix=mass_matrix, sym=cn.symbols)
end

symbolsof(cn::CompositeNode) = cn.symbols

dimension(cn::CompositeNode) = cn.total_dim

export CompositeNode
