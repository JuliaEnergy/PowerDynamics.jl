function total_current(e_s, e_d)
    # Keeping with the convention of negative sign for outging current
    current = 0.0im
    for e in e_s
        current += e[1] + e[2]*im
    end
    for e in e_d
        current -= e[1] + e[2]*im
    end
    current
end
