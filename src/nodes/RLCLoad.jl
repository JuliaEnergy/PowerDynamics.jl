# (C) 2018 Potsdam Institute for Climate Impact Research, authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)

@doc doc"""
```Julia
RLCLoad(R,L,C)
```
A node type that represents the RLC load model according to
"Power Systems Electromagnetic Transients Simulation",
Neville Watson and Jos Arrillaga, IET 2007, p.59, eq. (3.47)

# Keyword Arguments
- `R`: resistance in [?]
- `L`: inductance in [?]
- `C`: capacitance in [?]


# Mathematical Representation
```math
	\dfrac{du_C}{dt} = \frac{1}{C}i_L(t)\\
    \dfrac{di_L}{dt} = -\frac{R}{L} i_L(t)+\frac{1}{L} u(t)
```


"""
@DynamicNode RLCLoad(R,L,C)  begin
MassMatrix(m_int = [true, true,true])
end begin
    @assert R > 0 "Resistance should be >0"
    @assert L > 0 "Inductance should be >0"
    @assert C > 0 "Capacitance should be >0"

end [[u_C, du_C],[i_L,di_L],[ω,dω]] begin
    du_C = 1/C *i_L + 1*im*ω*u_c
    di_L = R/L*i_L+1/L*u + 1*im*ω*i_L
    du = i_L - i # = enforcing current conservation
    dω = 0
end

export RLCLoad
