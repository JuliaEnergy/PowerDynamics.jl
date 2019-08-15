# (C) 2018 Potsdam Institute for Climate Impact Research, authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)

export RLLoad

@doc doc"""
```Julia
RLLoad(R,L,\omega_r)
```

A node type that represents the RL load model
# Keyword Arguments

- `R`: resistance in [?]
- `L`: inductance in [?]
- `\omega_r`: reference frequency of the system


# Mathematical Representation
```math
    dfrac{di}{dt} = 1/L(i \omega_r L i + R i - u)
```
"""
@DynamicNode RLNode(L, R, ω_r) begin
    MassMatrix(m_u = false, m_int = [true, true])
end  begin
    Y = (im * ω_r * L + R)/L
    L_inv = 1/L
end [[i_r, di_r], [i_i, di_i]] begin
    i_int = complex(i_r, i_i)
    du = i - i_int

    di_int = Y * i_int - L_inv * u
    di_r = real(di_int)
    di_i = imag(di_int)
end
