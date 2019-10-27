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
@DynamicNode FESS(u_dc_ref,B,C,J,ω_m_ref,k_PLL,f,L_g,L_m,P_n,Ψ_m,R_m,R_g,K_m1,K_m2,K_m3,K_m4,K_g1,K_g2,K_g3,K_g4) begin
    MassMatrix(m_u = false,m_int = [true,true,true,true,true,true,true,true,true,true,true,true,true,true])
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
    @assert ω_m_ref>=0
end [[θ_m,dθ_m],[θ_g,dθ_g],[ω_m,dω_m],[i_dm,di_dm],[i_qm,di_qm],[e_I_idm,de_I_idm],[e_I_iqm,de_I_iqm],[e_Iω,de_Iω],[i_dg,di_dg],[i_qg,di_qg],[e_I_idg,de_I_idg],[e_I_iqg,de_I_iqg],[u_Iudc,du_Iudc],[u_dc,du_dc]] begin

    u_m = u*exp(-1im*θ_m)
    dθ_m = ω_m-ω_m_ref #v_q*k_PLL # TODO check: v_d is increased with decreasing θ which means i_qg is decreased with increasing θ
    #ω = dθ_PLL# TODO or is it not the frequency deviation but (1+dθ_PLL)*2*π*f?

    u_g = u*exp(-1im*θ_g)
    dθ_g = imag(u_g)*k_PLL
    ω=dθ_g

    i_dm_ref=0#??
    i_qg_ref=0#???

    # simplifications
    v_g = -1im*L_g*ω*(i_dg+1im*i_qg)+u_m-u_g
    v_m = -1im*L_m*P_n*ω_m*(i_dm+1im*i_qm)+u_m-1im*P_n*Ψ_m*ω_m


    #machine-side converter equations:
    di_dm = 1/L_m*(-R_m*i_dm+real(v_m))
    di_qm = 1/L_m*(-R_m*i_qm+imag(v_m))

    dω_m = 1/J*(Ψ_m*i_qm-D*ω_m)
    # PI outer speed controller
    e_ω = ω_m-ω_m_ref
    de_Iω = e_ω
    i_qm_ref = K_m3*e_ω+K_m4*e_Iω
    # PI inner current controller
    e_idm = i_dm_ref-i_dm
    de_I_idm = e_idm
    e_iqm = i_qm_ref-i_qm
    de_I_iqm = e_iqm
    v_m = (K_m1*e_idm + K_m2*e_I_idm)+1im*(K_m1*e_iqm + K_m2*e_I_iqm) # _ref?

    # grid-side converter
    di_dg = 1/L_g*(-R_g*i_dg+ω*L_g*i_qg+real(u_m)-real(u_g))
    di_qg = 1/L_g*(-R_g*i_q-ω*L_g*i_dg+imag(u_m)-imag(u_g))

    du_dc = 3/(2*C)*abs(i_dg+1im*i_dg)
    # PI outer speed controller
    e_udc = u_dc_ref-u_dc
    du_Iudc = e_udc
    i_dg_ref = K_g3*e_udc+K_m4*u_Iudc
    # PI inner current controller
    e_iqg = i_qg_ref-i_qg
    de_I_iqg = e_iqg
    e_idg = i_dg_ref-i_dg
    de_I_idg = e_idg
    u_g = (K_g1*e_idg+K_g2*e_I_idg)+1im*(K_g1*e_iqg+K_g2*e_I_iqg)

    println("i: ",i)
    println("i_dg: ",i_dg)

    du = i-(i_dg+1im*i_qg)*exp(1im*θ_g)
end

export FESS

#function PIControl(e)
#    deI = e
#    [e,eI]
#end
