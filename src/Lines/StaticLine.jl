

@Line StaticLine(complex_admittance_edge)
    #source_voltage = u_s
    #destination_voltage = u_d
    complex_current = admittance * (destination_voltage - source_voltage)
end

# TODO Macro shall generate this function
function (cae::complex_admittance_edge!)(e,v_s,v_d,p,t)
    source_voltage = v_s[1] + v_s[2]*im
    destination_voltage = v_d[1] + v_d[2]*im
    # If current is flowing away from the source, it is negative at the source.
    complex_current = cae.admittance * (destination_voltage - source_voltage)
    e[1] = real(complex_current)
    e[2] = imag(complex_current)
    nothing
end
