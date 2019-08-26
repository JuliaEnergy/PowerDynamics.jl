"""
A transformer representation uses the Π model,
assuming an ideal transformer in series with an admittance.
The admittance is here taken to be on the high-voltage side.
"""
@Line Transformer(y, t_ratio) begin

    # If current is flowing away from the source, it is negative at the source.
    voltage_vector = [source_voltage,destination_voltage]
    # Π[:, [k, m]] ./ pu # normalise to per unit admittance
    Y = PiModel(y, 0, 0, t_ratio, 1)
    current_vector = Y * voltage_vector
end


export Transformer
