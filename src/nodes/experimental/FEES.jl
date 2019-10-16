@doc doc"""
```Julia
FESS()
```

This implementation of a Virtual Synchronous Machine is undertaken according to TU Clausthal.

Additionally to ``u``, the global network-side voltage, it has the internal dynamic variables:
* `i_{dm}`
* `i_{qm}`
* `\omega_m`
* `v_{qm}`
* `i^*_{qm}`???
* `i_{dg}`
* `i_{qg}`
* `v_{dg}`
* `i^*_{dg}`???
* `ω_r`: frequeny deviation in rad/s.

# Keyword arguments are:
* `L_g/L_m`:
* `R_g/R_m`:
* `E`
* `P_n`
* `\omega^*_m`
* `K_{m1}/K_{m2}`
* `K_{m3}/K_{m4}`
* `K_{g1}/K_{g2}`
* `K_{g3}/K_{g4}`

# Mathematical Representation
Using `FESS` applies the equations for the stator circuit:
```math
    i_d = 1/R_d*(u_d-\frac{dΨ_d}{dt}+ω*Ψ_q)\\
    i_q = 1/R_q*(u_q+\frac{dΨ_q}{dt}-ω*Ψ_d)
```
"""
@DynamicNode FESS() begin
end  begin
end [[θ,dθ],[Ψ_d,dΨ_d],[Ψ_q,dΨ_q],[Ψ_D,dΨ_D],[Ψ_Q,dΨ_Q],[Ψ_e,dΨ_e],[ω,dω]] begin
    
end

export FESS
