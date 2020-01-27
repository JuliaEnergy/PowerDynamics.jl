"""
```Julia
ConnectorLine(from, to)
```
A dummy line that represents a static line with an admittance of Y=0.
For numeric reasons the line is not exactly zero but Z=1e-8+1e-8*j
"""
@Line ConnectorLine(from, to) begin
    # If current is flowing away from the source, it is negative at the source.
    # the current flowing in and out of the line is the same: I_mk=I_km
    # hence, the current vector becomes only one complex current
    Y = (1/(1e-10+1im*1e-10))# TODO
    complex_current = Y * (destination_voltage - source_voltage)
    current_vector = [complex_current,complex_current]
end

export ConnectorLine
