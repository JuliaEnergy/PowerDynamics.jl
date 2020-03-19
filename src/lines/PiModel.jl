function PiModel(Y, Y_shunt_km, Y_shunt_mk, T_km, T_mk)
    """
    Implementation of the unified branch model with our sign conventions.
    See also the Chapter 2 in
      Göran Andersson, _Power System Analysis_, Lecture 227-0526-00, ITET ETH Zurich, 2012
    Assumptions:
    * the line admittance is symmetric
    """
    Π = zeros(Complex{Float64}, 2, 2)
    Π[1, 1] = - abs2(T_km) * (Y + Y_shunt_km) # Our sign convention is opposite for the source of the edge
    Π[1, 2] = conj(T_km) * T_mk * Y # Our sign convention is opposite for the source of the edge
    Π[2, 1] = - conj(T_mk) * T_km * Y
    Π[2, 2] = abs2(T_mk) * (Y + Y_shunt_mk)
    Π
end
