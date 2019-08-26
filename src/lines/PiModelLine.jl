"""
A line modelled according to the PI-Model.
"""
@Line PiModelLine(from, to, y, y_shunt_km, y_shunt_mk) begin
    Y = PiModel(y, y_shunt_km, y_shunt_mk, 1, 1)
end begin
    # If current is flowing away from the source, it is negative at the source.
    voltage_vector = [source_voltage,destination_voltage]
    current_vector = Y * voltage_vector
end

export PiModelLine
