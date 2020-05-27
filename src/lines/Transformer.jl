"""
```Julia
    Transformer(from, to, y, t_ratio)
```

assuming an ideal transformer in series with an admittance.
The representation uses the Π model.

# Mathematical Representation
The voltage transforms as:
```math
    u_{to} = t_{ratio} u_{from}
```

# Arguments

- `from` : start node
- `to` : end node
- `y`: transformer admittance
- `t_ratio`: transformation ration

# Assumptions

The admittance is here taken to be on the high-voltage side.
"""
@Line Transformer(from, to, y, t_ratio) begin

    # If current is flowing away from the source, it is negative at the source.
    voltage_vector = [source_voltage,destination_voltage]
    Y = PiModel(y, 0, 0, t_ratio, 1)
    current_vector = Y * voltage_vector
end

export Transformer
