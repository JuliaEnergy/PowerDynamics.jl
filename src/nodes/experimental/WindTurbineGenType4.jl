@doc doc"""
```Julia
WindTurbineGenType4(;I_n,k_PLL,f,f_s,T_m,k_P,τ_ω)
```

```
"""

@DynamicNode WindTurbineGenType4(D,K_PLL,Q_ref,C,J,P,ω_rref,K_Q,K_v,K_g1,K_g2,K_r1,K_r2) begin
    MassMatrix(m_u = false,m_int = [true,true,true,true,true,true,true,true,true,true])
end  begin
    @assert J>=0
    @assert K_g1>=0
    @assert K_g2>=0
    @assert K_r1 >=0
    @assert K_r2 >=0
    @assert ω_rref>=0
    @assert C>=0
    @assert K_PLL>=0
end [[θ_PLL,dθ_PLL],[P_inertia,dP_inertia],[e_Idθ,de_Idθ],[e_IP,de_IP],[u_dc,du_dc],[ω_r,dω_r],[i_q,di_q],[u_tref,du_tref],[e_Iu,de_Iu],[e_IV,de_IV]] begin#
    function PI_control(e_I,e,K_P,K_I)
        @assert K_P>=0
        @assert K_I>=0
        de_I=e
        u = K_P*e+K_I*e_I
        return [de_I,u]
    end
    u_dq = u*exp(-1im*θ_PLL)
    v_d = real(u_dq)
    v_q = imag(u_dq)
    #i_q = 0.


    de_Idθ,ω_PLL=PI_control(e_Idθ,v_q,K_PLL,K_PLL)#v_q*K_PLL
    #println("ω_PLL",ω_PLL)

    dθ_PLL=ω_PLL
    τ_L=1
    T_H =5.5
    K_f=10
    f=50
    f_g = ω_PLL#*2π*50
    #dω = 1/τ_L*(-ω + f_g)
    #dz = 1/T_H*(-z+ω)
    #t_ω=K*dz
    #t_ω = (-ω*K_f-t_Iω/T_H)
    #dt_Iω = t_ω
    #println("t_ω: ",t_ω)
    #println("dω: ",dω)
    #println("ω_PLL",ω_PLL)
    #println("ω:",ω)

    #println("dθ_PLL: ", dθ_PLL)

    #ω = ω_r-1-dθ_PLL
    #de_IV,i_d=PI_control(e_IV,(u_dcref-u_dc),K_g1,K_g2)


    # speed control:
    de_IP,t_e=PI_control(e_IP,(ω_rref-ω_PLL-ω_r),K_r1,K_r2)

    f_trigger=-0.3
    f_m = ω_PLL*2π
    if f_m<f_trigger
        dP_inertia = (f_trigger-f_m)*0.1*P
    else dP_inertia = 0.
    end

    t_m = P/ω_r
    #t_e_ref=t_e+dP_inertia/ω_r
    #de_IV,i_d=PI_control(e_IV,(t_e_ref-t_e),K_g1,K_g2)
    p_in = t_e*ω_r

    # DC voltage control:
    #di_d = (u_dc_ref-u_dc)
    println("dP_inertia",dP_inertia)
    t_e_ref = (P+dP_inertia)/ω_r
    println("t_e:",t_e)
    println("t_e_ref:",t_e_ref)
    println("ω_PLL",ω_PLL)
    de_Iu,u_dc_ref = PI_control(e_Iu,(t_e_ref-t_e),1,1)
    de_IV,i_d=PI_control(e_IV,(u_dc_ref-u_dc),K_g1,K_g2)
    i_dq = (i_d+1im*i_q)

    s_e = u_dq*conj(i_dq)# s=(v_d*i_d+v_q*i_q)+j(v_q*i_d-v_d*i_q)
    p_e = real(s_e)
    q_e = imag(s_e)

    du_dc = 1/(u_dc*C)*(p_in-p_e)

    #de_IV,i_d=PI_control(e_IV,(u_dcref-u_d,c),K_g1,K_g2)


    dω_r = 1/J*(t_m-t_e)#-D*ω_r)

    # reactive power control
    u_t=abs(u)
    du_tref = K_Q*(Q_ref-q_e)
    di_q = K_v*(u_tref-u_t)

    du = i - i_dq*exp(1im*θ_PLL)
end

export WindTurbineGenType4
