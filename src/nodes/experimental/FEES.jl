@doc doc"""
```Julia
FESS()
```

This implementation of a Virtual Synchronous Machine is undertaken according to TU Clausthal.

Additionally to ``u``, the global network-side voltage, it has the internal dynamic variables:
* `\theta`
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
@DynamicNode FESS(u_dc_ref,C,J,ω_m_ref,k_PLL,f,L_g,L_m,P_n,Ψ_m,R_m,R_g,K_m1,K_m2,K_m3,K_m4,K_g1,K_g2,K_g3,K_g4) begin
    MassMatrix(m_u = false,m_int = [true,true,true,true,true,true,true,true,true,true,true,true])
end  begin
    @assert k_PLL>=0
    @assert f>=0
    @assert C>=0
    @assert L_g >=0
    @assert L_m >=0
    @assert P_n >0
    @assert R_m >=0
    @assert R_g >=0
    @assert K_m1 >=0
    @assert K_m2 >=0
    @assert K_m3 >=0
    @assert K_m4 >=0
    @assert K_g1 >=0
    @assert K_g2 >=0
    @assert K_g3 >=0
    @assert K_g4 >=0
end [[θ_PLL,dθ_PLL],[i_dm,di_dm],[i_qm,di_qm],[ω_m,dω_m],[v_qm,dv_qm],[i_qm_ref,di_qm_ref],[i_dg,di_dg],[i_qg,di_qg],[u_dc,du_dc],[u_dg,du_dg],[u_qg,du_qg],[i_dg_ref,di_dg_ref]] begin

    i_dm_ref=0#??
    i_qg_ref=0#???

    u_x = real(u)
    u_y = imag(u)

    v_d = u_x*cos(θ_PLL)+u_y*sin(θ_PLL)
    v_q = u_x*sin(θ_PLL)-u_y*cos(θ_PLL)

    dθ_PLL = -v_q*k_PLL
    ω = (1+dθ_PLL)*2*π*f

    # simplifications
    u_dg = L_g*ω*i_qg+v_d#-v_dg
    u_qg = -L_g*ω*i_dg+v_q#-v_qg
    v_dm = L_m*P_n*ω_m*i_qm +v_d
    v_qm = -L_m*P_n*ω_m*i_dm - P_n*Ψ_m +v_q

    #machine-side converter equations:
    di_dm = 1/L_m*(-R_m*i_dm+v_dm)
    di_qm = 1/L_m*(-R_m*i_qm+v_qm)
    dω_m = 3*P_n/(2*J)*Ψ_m*i_qm
    # PI outer speed controller
    di_qm_ref = K_m3*dω_m+K_m4*(ω_m-ω_m_ref)
    # PI inner current controller
    dv_qm = K_m1*(di_qm_ref-di_qm)+K_m2*(i_qm_ref-i_qm) # _ref?
    dv_dm = K_m1*(-di_dm)+K_m2*(i_dm_ref-i_dm) #_ref?

    # grid-side converter
    di_dg = 1/L_g*(-R_g*i_dg+u_dg)
    di_qg = 1/L_g*(-R_g*i_qg+u_qg)
    du_dc = 3/(2*C)*(i_dg+i_qg)
    # PI outer speed controller
    di_dg_ref = K_g3*(-du_dc)+K_m4*(u_dc_ref-u_dc)
    # PI inner current controller
    du_dg = K_g1*(di_dg_ref-di_dg)+K_g2*(i_dg_ref-i_dg)
    du_qg = K_g1*(-di_qg)+K_g2*(i_qg_ref-i_qg)


    i_x = i_dg*cos(θ_PLL) + i_qg*sin(θ_PLL)
    i_y = i_dg*sin(θ_PLL) - i_qg*cos(θ_PLL)
    du = i-(i_x+1im*i_y)
end

export FESS
