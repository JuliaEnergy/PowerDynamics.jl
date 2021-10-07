@doc doc"""
```Julia
FluctuationNode(; P,Q)
```

"""
struct FluctuationNode <: AbstractNode
    P
    Q
end

function construct_vertex(fn::FluctuationNode)

    sym = symbolsof(fn)
    dim = dimension(fn)
    mass_matrix = zeros(Int64,2,2) |> Diagonal

    P = fn.P; Q = fn.Q

    function rhs!(dx, x, edges, p, t)
        i = total_current(edges)
        u = x[1] + x[2] * im
        s = u*conj(i)
        dx[1] = P(t) - real(s)
        dx[2] = Q(t) - imag(s)
        nothing
    end
    
    ODEVertex(f! = rhs!, dim=dim, mass_matrix=mass_matrix, sym=sym)

end

symbolsof(fn::FluctuationNode) = [:u_r, :u_i]

dimension(fn::FluctuationNode) = 2

export FluctuationNode