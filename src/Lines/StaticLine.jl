@Line StaticLine(Y) begin
    # If current is flowing away from the source, it is negative at the source.
    # the current flowing in and out of the line is the same: I_mk=I_km
    # hence, the current vector becomes only one complex current
    complex_current = Y * (destination_voltage - source_voltage)
    current_vector = [complex_current,complex_current]
end

##### Macro Generated functions should look like:
#function rhs!()
#    source_voltage = v_s[1] + v_s[2]*im
#    destination_voltage = v_d[1] + v_d[2]*im
#    # If current is flowing away from the source, it is negative at the source.
#    complex_current = admittance * (destination_voltage - source_voltage)
#    e[1] = real(complex_current)
#    e[2] = imag(complex_current)
#end

#function construct_line(sl:: StaticLine)
#    return StaticEdge(f! = rhs!,dim = 2)
#end

export StaticLine
