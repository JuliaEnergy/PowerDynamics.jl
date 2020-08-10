# (C) 2018 Potsdam Institute for Climate Impact Research, authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)

@doc doc"""
```Julia
    VoltageDependentLoad(;P, Q, U, q, l)
```

A node type that locally fixes the active (``P``) and reactive power (``Q``) output of the node.

# Keyword Arguments
- `P`: active power demand
- `Q`: reactive power demand
- `U` : the voltage set point
- `q` : relative share of quadratic voltage dependence 
- `l` : relative share of linear voltage dependence 

# Mathematical Representation
Using `VoltageDependentLoad` for node ``a`` applies the equation
```math
0 = S_a ( q (u_a/U_a)^2 + l (u_a/U_a) + 1 - q - l) - u_a \cdot i_a^*.
```
"""
@DynamicNode VoltageDependentLoad(P, Q, U, q, l) begin
    MassMatrix()
end  begin
    @assert 0 <= q <= 1
    @assert 0 <= l <= 1
end [] begin
    s = u * conj(i)
    u_rel = u / U
    du = complex(P, Q) * (q * u_rel^2 + l * u_rel + 1 - q - l) - s
end

export VoltageDependentLoad