using SparseArrays

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

@__doc__ struct CompositeNode{T}
    CurrentNodes::T
    idims::Array{Int, 1}
    total_dim::Int
    idxs::Array{Array{Int,1},1}
    symbols::Array{Symbol, 1}
    function CompositeNode(CurrentNodes)
        symbols = [:u_r, :u_i]
        idims = Array{Int,1}()

        for current_node in CurrentNodes
            append!(symbols, symbolsof(current_node)[3:end])
            append!(idims, dimension(current_node) - 2)
        end

        idxs = Array{Array{Int,1},1}([])
        offset = 2
        for idim in idims
            append!(idxs, [[[1,2]; (offset+1):(offset+idim)]])
            offset += idim
        end

        total_dim = 2 + sum(idims)
        new{typeof(CurrentNodes)}(CurrentNodes,
        idims,
        total_dim,
        idxs,
        symbols)

    end
end

function construct_vertex(cn::CompositeNode{T}) where T

    CN_Type = typeof(cn)

    current_nodes = [construct_vertex(cnode) for cnode in cn.CurrentNodes]

    cfs = [cn.f! for cn in current_nodes]

    mass_matrices = [cn.mass_matrix for cn in current_nodes]

    mass_matrix = spzeros(cn.total_dim, cn.total_dim)

    for (i, mm) in enumerate(mass_matrices)
        mass_matrix[cn.idxs[i], cn.idxs[i]] .= mm
    end

    mass_matrix[1, :] .= 0.
    mass_matrix[2, :] .= 0.
    mass_matrix[:, 1] .= 0.
    mass_matrix[:, 2] .= 0.

    function rhs!(dx, x, e_s, e_d, p, t)
        i = total_current(e_s, e_d)
        u = x[1] + x[2] * im

        i_injected = 0.0 + 0.0im

        @views begin
        for (cf, idx) in zip(cfs, cn.idxs)
            cf(dx[idx], x[idx], e_s, e_d, p, t)
            # The current nodes are assumed to have du = i - I_r
            i_injected += i - complex(dx[1], dx[2])
        end
        end

        du = i - i_injected
        dx[1] = real(du)
        dx[2] = imag(du)
    end

    ODEVertex(f! = rhs!, dim=cn.total_dim, mass_matrix=mass_matrix, sym=cn.symbols)
end

symbolsof(cn) = cn.symbols

dimension(cn) = cn.total_dim

export CompositeNode
