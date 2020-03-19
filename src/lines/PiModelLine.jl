"""
```Julia
PiModelLine(from, to, Y, Y_shunt_km, Y_shunt_mk)
```
A line modelled according to the PI-Model.
"""
@Line PiModelLine(from, to, Y, Y_shunt_km, Y_shunt_mk) begin
    Y_Pi = PiModel(Y, Y_shunt_km, Y_shunt_mk, 1, 1)
end begin
    # If current is flowing away from the source, it is negative at the source.
    voltage_vector = [source_voltage,destination_voltage]
    current_vector = Y_Pi * voltage_vector
end

export PiModelLine
