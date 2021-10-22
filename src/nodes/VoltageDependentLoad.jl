# (C) 2018 Potsdam Institute for Climate Impact Research, authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)

@doc doc"""
```Julia
    VoltageDependentLoad(;P, Q, U, A, B)
```

A node type that locally fixes the active (``P``) and reactive power (``Q``) output of the node.

# Keyword Arguments
- `P`: active power demand
- `Q`: reactive power demand
- `U` : the voltage set point
- `A` : relative share of quadratic voltage dependence 
- `B` : relative share of linear voltage dependence 

# Mathematical Representation
Using `VoltageDependentLoad` for node ``a`` applies the equation
```math
0 = S_a ( A \cdot (u_a/U_a)^2 + B \cdot (u_a/U_a) + 1 - A - B) - u_a \cdot i_a^*.
```
"""
@DynamicNode VoltageDependentLoad(P, Q, U, A, B) begin
    MassMatrix()
end  begin
    @assert 0 <= A <= 1
    @assert 0 <= B <= 1
    @assert isreal(U)
end [] begin
    s = u * conj(i)
    u_rel = abs(u) / U
    du = complex(P, Q) * (A * u_rel^2 + B * u_rel + 1 - A - B) - s
end

export VoltageDependentLoad
