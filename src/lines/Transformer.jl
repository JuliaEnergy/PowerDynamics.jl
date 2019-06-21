@Line Transformer(y, t_ratio) begin
    """
    This transformer representation uses the Π model,
    assuming an ideal transformer in series with an admittance.
    The admittance is here taken to be on the high-voltage side.
    """
    # If current is flowing away from the source, it is negative at the source.
    voltage_vector = [source_voltage,destination_voltage]
    # Π[:, [k, m]] ./ pu # normalise to per unit admittance
    Y = PiModel(y, 0, 0, t_ratio, 1)
    current_vector = Y * voltage_vector
end


export Transformer

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

#=
function pi_model(y, y_shunt_km, y_shunt_mk, t_km, t_mk)
    """
    Implementation of the unified Π model.
    See also the Chapter 2 in
      Göran Andersson, _Power System Analysis_, Lecture 227-0526-00, ITET ETH Zurich, 2012
    Assumptions:
    * the line admittance is symmetric
    """
    Π = spzeros(Complex{Float64}, 2, 2)
    Π[1, 1] = abs2(t_km) * (y + y_shunt_km)
    Π[1, 2] = - conj(t_km) * t_mk * y
    Π[2, 1] = - conj(t_mk) * t_km * y
    Π[2, 2] = abs2(t_mk) * (y + y_shunt_mk)
    Π
end

export PiModelLine
=#
