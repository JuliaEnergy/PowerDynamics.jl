function total_current(edges)
    # Keeping with the convention of negative sign for outging current
    current = 0.0im
    for e in edges
        current += e[1] + e[2]*im
    end
    current
end
