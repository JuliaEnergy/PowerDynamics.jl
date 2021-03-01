function total_current(edges)
    # ND convention: all flows entering a node are positive
    current = 0.0im
    @inbounds for e in edges
        current += e[1] + e[2]*im
    end
    current
end
