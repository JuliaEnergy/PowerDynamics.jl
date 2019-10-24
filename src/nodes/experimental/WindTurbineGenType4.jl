@doc doc"""
```Julia
WindTurbineGenType4(;I_n,k_PLL,f,f_s,T_m,k_P,τ_ω)
```

```
"""
PI_control(e,u)

@DynamicNode WindTurbineGenType4(K_PLL,Q_ref,C,J,P,ω_rref,u_dcref,K_Q,K_v,K_g1,K_g2,K_r1,K_r2) begin
    MassMatrix(m_u = false,m_int = [true,true,true,true,true,true, true])
end  begin
    @assert J>=0
    @assert K_g1>=0
    @assert K_g2>=0
    @assert K_r1 >=0
    @assert K_r2 >=0
    @assert ω_rref>0
    @assert C>=0
end [[θ,dθ],[t_e,dt_e],[u_dc,du_dc],[i_dref,di_dref],[i_qref,di_qref],[u_tref,du_tref],[ω_r,dω_r]] begin
    u_dq = u*exp(-1im*θ)
    v_d = real(u_dq)
    v_q = imag(u_dq)
    dθ = v_q*K_PLL

    t_m = P/ω_r
    dω_r = 1/J*(t_m-t_e)

    dt_e = K_r1*dω_r+K_r2*(ω_r-ω_rref)
    p_in = t_e*ω_r
    s_e = u*conj(i)# s=(v_d*i_d+v_q*i_q)+j(v_q*i_d-v_d*i_q)
    p_e = real(s_e)
    q_e = imag(s_e)

    du_dc = 1/C*(p_in-p_e)

    u_t=v_d#TODO: discuss!!!
    di_dref = K_g1*du_dc + K_g2*(u_dc-u_dcref)
    du_tref = K_Q*(Q_ref-q_e)
    di_qref = K_v*(u_tref-u_t)

    du = i - (i_dref+1im*i_qref)*exp(1im*θ)
end

export WindTurbineGenType4
