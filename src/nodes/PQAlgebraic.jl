# (C) 2018 Potsdam Institute for Climate Impact Research, authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)

@doc doc"""
```Julia
    PQAlgebraic(; P,Q)
```

A node type that locally fixes the active (``P``) and reactive power (``Q``) output of the node.

# Keyword Arguments
- `P`: active power set point
- `Q`: reactive power set point

# Mathematical Representation
Using `PQAlgebraic` for node ``a`` applies the equation
```math
0 = (P_a + jQ_a) - u_a \cdot i_a^*.
```
"""
@DynamicNode PQAlgebraic(P,Q) begin
    MassMatrix()
end  begin
    # no prep
end [] begin
    s = u*conj(i)
    du = complex(P, Q) - s
end

export PQAlgebraic
