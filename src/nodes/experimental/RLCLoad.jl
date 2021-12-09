# (C) 2018 Potsdam Institute for Climate Impact Research, authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)

@doc doc"""
```Julia
RLCLoad(R,L,C)
```
EXPERIMENTAL
A node type that represents the RLC load model according to
"Power Systems Electromagnetic Transients Simulation",
Neville Watson and Jos Arrillaga, IET 2007, p.59, eq. (3.47)

# Keyword Arguments
- `R`: resistance
- `L`: inductance
- `C`: capacitance


# Mathematical Representation
```math
	\dfrac{du_C}{dt} = \frac{1}{C}i_L(t)\\
    \dfrac{di_L}{dt} = -\frac{R}{L} i_L(t)+\frac{1}{L} u(t)
```
"""
@DynamicNode RLCLoad(R,L,C)  begin
MassMatrix(m_int = [true, true,true,true,true])
end begin
    @assert R > 0 "Resistance should be >0"
    @assert L > 0 "Inductance should be >0"
    @assert C > 0 "Capacitance should be >0"

end [[u_Cr, du_Cr],[u_Ci, du_Ci],[i_Lr,di_Lr],[i_Li,di_Li],[ω,dω]] begin
    du_C = 1/C *(i_Lr+1im*i_Li) + 1*im*ω*(u_Cr+1im*u_Ci)
    du_Cr = real(du_C)
    du_Ci = imag(du_C)
    di_L = -R/L*(i_Lr+1im*i_Li)+1/L*u -1/L*(u_Cr+1im*u_Ci)+ 1*im*ω*(i_Lr+1im*i_Li)
    di_Lr = real(di_L)
    di_Li = imag(di_L)
    du = (i_Lr+1im*i_Li) - i # = enforcing current conservation
    dω = 0
end

export RLCLoad
