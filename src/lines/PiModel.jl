function PiModel(y, y_shunt_km, y_shunt_mk, t_km, t_mk)
    """
    Implementation of the unified branch model.
    See also the Chapter 2 in
      Göran Andersson, _Power System Analysis_, Lecture 227-0526-00, ITET ETH Zurich, 2012
    Assumptions:
    * the line admittance is symmetric
    """
    Π = zeros(Complex{Float64}, 2, 2)
    Π[1, 1] = abs2(t_km) * (y + y_shunt_km)
    Π[1, 2] = - conj(t_km) * t_mk * y
    Π[2, 1] = - conj(t_mk) * t_km * y
    Π[2, 2] = abs2(t_mk) * (y + y_shunt_mk)
    Π
end
