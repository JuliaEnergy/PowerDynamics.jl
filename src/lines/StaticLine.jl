"""
```Julia
StaticLine(from, to, Y)
```
A static model that represents a line with an admittance Y.
"""
@Line StaticLine(from, to, Y) begin
    # If current is flowing away from the source, it is negative at the source.
    # the current flowing in and out of the line is the same: I_mk=I_km
    # hence, the current vector becomes only one complex current
    complex_current = Y * (destination_voltage - source_voltage)
    current_vector = [complex_current,complex_current]
end

export StaticLine
