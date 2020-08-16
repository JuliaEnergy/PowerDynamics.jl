"""
```Julia
    PiModelLine(from, to, y, y_shunt_km, y_shunt_mk)
```
A line modelled according to the PI-Model.

See also the Chapter 2 in
  Göran Andersson, _Power System Analysis_, Lecture 227-0526-00, ITET ETH Zurich, 2012

# Arguments
- `from` : node `k`
- `to` : node `m`
- `y`: admittance of line between `k` and `m`
- `y_shunt_km`: shunt admittance at the end connected to node `k`
- `y_shunt_mk`: shunt admittance at the end connected to node `m`

# Assumptions:
- the line admittance is symmetric
"""
@Line PiModelLine(from, to, y, y_shunt_km, y_shunt_mk) begin
    Y = PiModel(y, y_shunt_km, y_shunt_mk, 1, 1)
end begin
    # If current is flowing away from the source, it is negative at the source.
    voltage_vector = [source_voltage,destination_voltage]
    current_vector = Y * voltage_vector
end

export PiModelLine
