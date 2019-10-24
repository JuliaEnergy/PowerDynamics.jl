@doc doc"""
```Julia
WindTurbineGenType4(;I_n,k_PLL,f,f_s,T_m,k_P,τ_ω)
```

```
"""
function PI_control(e,u,K_P,K_I)
    @assert K_P>=0
    @assert K_I>=0
    de_I=e
    u = K_P*e+K_I*e_I
    return [de_I,u]
end

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
end [[θ,dθ],[e_IP,de_IP],[e_IV,de_IV],[u_dc,du_dc],[i_q,di_q],[u_tref,du_tref],[ω_r,dω_r]] begin
    u_dq = u*exp(-1im*θ)
    v_d = real(u_dq)
    v_q = imag(u_dq)
    dθ = v_q*K_PLL

    s_e = u*conj(i)# s=(v_d*i_d+v_q*i_q)+j(v_q*i_d-v_d*i_q)
    p_e = real(s_e)
    q_e = imag(s_e)

    # PI speed control:
    e_p = -(ω_rref-ω_r)
    de_Ip = e_p
    t_e = K_r1*e_p+K_r2*e_IP
    p_in =t_e/ω_r

    t_m = P/ω_r
    dω_r = 1/J*(t_m-t_e)

    du_dc = 1/C*(p_in-p_e)

    u_t=v_d#TODO: discuss!!!

    # PI voltage control:
    e_V = -(u_dcref-u_dc)
    de_IV = e_V
    i_d = K_g1*e_V + K_g2*de_IV

    du_tref = K_Q*(Q_ref-q_e)
    di_q = K_v*(u_tref-u_t)
    println("v_d: ",v_d)

    du = i - (i_d+1im*i_q)*exp(1im*θ)
end

export WindTurbineGenType4
