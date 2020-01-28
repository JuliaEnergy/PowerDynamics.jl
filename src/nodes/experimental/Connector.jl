# (C) 2018 Potsdam Institute for Climate Impact Research, authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)

@doc doc"""
```Julia
Connector(;)
```

A dummy node type that represents a connector node without any dynamics or load.

It locally fixes the active (``P``) and reactive power (``Q``) to zero but due to numeric
reasons, P and Q are set to -1e-8.

# Keyword Arguments
- none

# Mathematical Representation
Using `Connector` for node ``a`` applies the equation
```math
0 = (-1e-8-1e-8*1im) - u_a \cdot i_a^*.
```
"""
@DynamicNode Connector() begin
    MassMatrix()
end  begin
    # no prep
end [] begin
    s = u*conj(i)
    du = s-(-0-1im)*1e-10
end

export Connector
