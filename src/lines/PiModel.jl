"""
```Julia
    PiModel(y, y_shunt_km, y_shunt_mk, t_km, t_mk)
```
Implementation of the unified branch model with our sign conventions.
See also the Chapter 2 in
  Göran Andersson, _Power System Analysis_, Lecture 227-0526-00, ITET ETH Zurich, 2012
Assumptions:
* the line admittance is symmetric

# Arguments
- `y`: line admittance
- `y_shunt_km`: shunt admittance at the end connected to node `k`
- `y_shunt_mk`: shunt admittance at the end connected to node `m`
- `t_km`: transformer ratio at the end connected to node `k`
- `t_mk`: transformer ratio at the end connected to node `m`
"""
function PiModel(y, y_shunt_km, y_shunt_mk, t_km, t_mk)
    Π = zeros(Complex{Float64}, 2, 2)
    Π[1, 1] = - abs2(t_km) * (y + y_shunt_km) # Our sign convention is opposite for the source of the edge
    Π[1, 2] = conj(t_km) * t_mk * y # Our sign convention is opposite for the source of the edge
    Π[2, 1] = - conj(t_mk) * t_km * y
    Π[2, 2] = abs2(t_mk) * (y + y_shunt_mk)
    Π
end
