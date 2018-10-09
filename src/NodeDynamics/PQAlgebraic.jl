# (C) 2018 Potsdam Institute for Climate Impact Research, authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)

@doc doc"""
```Julia
PQAlgebraic(;S)
```

A node type that locally fixes the active (``P``) and reactive power (``Q``) output of the node.

# Keyword Arguments
- `S = P + Q*im`: the complex power output

# Mathematical Representation
Using `PQAlgebraic` for node ``a`` applies the equation
```math
0 = S_a - u_a \cdot i_a^*.
```
"""
@DynamicNode PQAlgebraic(S) <: OrdinaryNodeDynamicsWithMass(m_u=false, m_int=no_internal_masses)  begin
end [] begin
    s = u*conj(i)
    du = S - s
end

export PQAlgebraic
