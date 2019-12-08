@doc doc"""
```Julia
CompositeNode(CurrentNodes)
```

A composite node of a selection of current sources.

"""
@__doc__ struct CompositeNode{T}
        CurrentNodes::T
    end
function construct_vertex(cn::CompositeNode{T})
    symbols = [:u_r, :u_i]
    idims = Array{Int,1}()

    for current_node in cn.CurrentNodes
        append!(symbols, symbolsof(typeof(current_node))[3:end])
        append!(idims, dimension(typeof(current_node)) - 2)
    end

    idxs = []
    offset = 2
    for idim in idims
        append!(idxs, [[1:2]; (offset+1):(offset+idim)])
        offset += idim
    end

    total_dim = 2 + sum(idims)

    @eval symbolsof(::CompositeNode{T}) = symbols

    @eval dimension(::CompositeNode{T}) = total_dim

    current_nodes = [construct_vertex(cnode) for cnode in cn.CurrentNodes]

    function rhs!(dx, x, e_s, e_d, p, t)
        i = total_current(e_s, e_d)
        u = x[1] + x[2] * im

        i_injected = 0. + im 0.

        @view begin
        for cn, idx in zip(current_nodes, idxs)
            cn(dx[idx], x[idx], e_s, e_p, p, t)
            # The current nodes are assumed to have du = i - I_r
            i_injected += i - du
        end
        end

        du = i - i_injected
        dx[1] = real(du)
        dx[2] = imag(du)
    end

    ODEVertex(f!=rhs!, dim=total_dim, mass_matrix=????, sym=symbols)
end
