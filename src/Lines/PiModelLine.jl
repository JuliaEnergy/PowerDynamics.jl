@Line PiModelLine(y, y_shunt_km, y_shunt_mk) begin
    # If current is flowing away from the source, it is negative at the source.
    voltage_vector = [source_voltage,destination_voltage]
    Y = PiModel(y, y_shunt_km, y_shunt_mk, 1, 1)
    current_vector = Y * voltage_vector
end
