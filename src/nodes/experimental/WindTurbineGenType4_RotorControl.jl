@doc doc"""
```Julia
WindTurbineGenType4_RotorControl(;T_L,ω_0,T_H,K_ω,K_PPLL,K_IPLL,Q_ref,C,J,P,ω_rref,u_dcref,K_Q,K_v,K_pv,K_iv,K_ptrq,K_itrq) begin
```

```
"""

@DynamicNode WindTurbineGenType4_RotorControl(T_L,ω_0,T_H,K_ω,K_PPLL,K_IPLL,Q_ref,C,J,P,ω_rref,u_dcref,K_Q,K_v,K_pv,K_iv,K_ptrq,K_itrq) begin
    MassMatrix(m_u = false,m_int = [true,true,true,true,true,true,true,true,true,true])
end  begin
    @assert J>=0
    @assert K_pv>=0
    @assert K_iv>=0
    @assert K_ptrq >=0
    @assert K_itrq >=0
    @assert C>=0
    @assert K_PPLL>=0
    @assert K_IPLL>=0
    @assert K_ω>=0
end [[θ_PLL,dθ_PLL],[e_Idθ,de_Idθ],[ω,dω],[e_IP,de_IP],[z,dz],[e_IV,de_IV],[u_dc,du_dc],[i_q,di_q],[u_tref,du_tref],[ω_r,dω_r]] begin
    # transformation from global to local dq-coordinates with θ_PLL
    u_dq = u*exp(-1im*θ_PLL)
    v_d = real(u_dq)
    v_q = imag(u_dq)

    # PLL dynamics for frequency measurement
    de_Idθ,ω_PLL=PIControl(e_Idθ,-v_q,K_PPLL,K_IPLL)
    dθ_PLL=ω_PLL

    # PI speed control:
    de_IP,t_eref=PIControl(e_IP,ω_rref-ω_PLL-ω_r,K_ptrq,K_itrq)

    # inertial respons emulation
    dω=1/T_L*(ω_PLL-ω)
    Δt_eref=1/(T_H*ω_0)*(-z*ω_0-ω)
    dz = Δt_eref
    t_e =(t_eref+Δt_eref*K_ω)

    # turbine dynamics
    t_m = P/ω_r
    p_in = t_e*ω_r
    dω_r = 1/J*(t_m-t_e)#-D*ω_r)

    # DC voltage control
    de_IV,i_d=PIControl(e_IV,(u_dcref-u_dc),K_pv,K_iv)
    i_dq = i_d+1im*i_q
    s_e = u_dq*conj(i_dq)# s=(v_d*i_d+v_q*i_q)+j(v_q*i_d-v_d*i_q)
    p_e = real(s_e)
    q_e = imag(s_e)
    du_dc = 1/C*(p_in-p_e)

    # reactive power control
    u_t=abs(u)
    du_tref = K_Q*(Q_ref-q_e)
    di_q = K_v*(u_tref-u_t)

    du = i - i_dq*exp(1im*θ_PLL)
end

export WindTurbineGenType4_RotorControl
