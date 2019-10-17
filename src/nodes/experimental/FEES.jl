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
end [[θ_PLL,dθ_PLL],[i_dm,di_dm],[i_qm,di_qm],[ω_m,dω_m],[v_qm,dv_qm],[i_qm_ref,di_qm_ref],[i_dg,di_dg],[i_qg,di_qg],[u_dc,du_dc],[v_dg,dv_dg],[i_dg_ref,di_dg_ref]] begin
    v_x = real(u)
    v_y = imag(u)
    v_d = v_x

    dθ_PLL = -i_qg*k_PLL

    # simplifications
    v_dg = L_g*ω*i_qg-u_dc+E
    v_qg = -L_g*ω*i_dg-u_dc
    v_dm = L_m*P_n*ω_m*i_qm - u_dc
    v_qm = -L_m*P_n*ω_m*i_dm - P_n*Ψ_m-u_dc

    #machine-side converter equations:
    di_dm = 1/L_m*(-R_m*i_dm+v_dm)
    di_qm = 1/L_m*(-R_m*i_qm+v_qm)
    dω_m = 3*P_n/(2*J)*Ψ_m*i_qm
    # PI inner current controller
    dv_qm = K_m1*(di_qm_ref-di_qm)+K_m2*(i_qm_ref-i_qm)
    # PI outer speed controller
    di_qm_ref = K_m3*dω_m+K_m4*(ω_m-ω_m_ref)

    # grid-side converter
    di_dg = 1/L_g*(-R_g*i_dg+v_dg)
    di_qg = 1/L_g*(-R_g*i_qm+v_qg)
    duc = 3/(2*C)*(i_dg+i_dq)
    # PI inner current controller
    dv_dg = K_g1*(di_dg_ref-di_dg)+K_g2*(i_dg_ref-i_dg)
    # PI outer speed controller
    di_dg_ref = K_g3*(du_dc)+K_m4*(u_dc_ref-u_dc)
end

export FESS
