"""
```Julia
    Transformer(from, to, Y, T_ratio)
```

assuming an ideal transformer in series with an admittance.
The representation uses the Π model.

# Mathematical Representation
The voltage transforms as:
```math
    u_{to} = T_{ratio} u_{from}
```

# Arguments

- `from` : start node
- `to` : end node
- `Y`: transformer admittance
- `T_ratio`: transformation ration

# Assumptions

The admittance is here taken to be on the high-voltage side.
"""
@Line Transformer(from, to, Y, T_ratio) begin

    # If current is flowing away from the source, it is negative at the source.
    voltage_vector = [source_voltage,destination_voltage]
    Π = PiModel(Y, 0, 0, T_ratio, 1)
    current_vector = Π * voltage_vector
end

export Transformer
