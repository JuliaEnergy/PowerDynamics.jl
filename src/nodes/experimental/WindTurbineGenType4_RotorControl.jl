@doc doc"""
```Julia
WindTurbineGenType4(;I_n,k_PLL,f,f_s,T_m,k_P,τ_ω)
```

```
"""

@DynamicNode WindTurbineGenType4_inertia(ΔP_max,K_PLL,Q_ref,C,J,P,ω_rref,u_dcref,K_Q,K_v,K_g1,K_g2,K_r1,K_r2) begin
    MassMatrix(m_u = false,m_int = [true,true,true,true,true,true,true])
end  begin
    @assert J>=0
    @assert K_g1>=0
    @assert K_g2>=0
    @assert K_r1 >=0
    @assert K_r2 >=0
    @assert ω_rref>=0
    @assert C>=0
    @assert K_PLL>=0
end [[θ_PLL,dθ_PLL],[e_Idθ,de_Idθ],[e_IP,de_IP],[e_IV,de_IV],[u_dc,du_dc],[i_q,di_q],[u_tref,du_tref],[ω_r,dω_r]] begin
    u_dq = u*exp(-1im*θ_PLL)
    v_d = real(u_dq)
    v_q = imag(u_dq)

    de_Idθ,ω_PLL=PI_control(e_Idθ,v_q,K_PLL,1)
    dθ_PLL=ω_PLL

    dθ_r = ω_r
    #println("dθ_PLL: ", dθ_PLL)

    if abs(f_m)>f_trigger
        ΔP =-k_P*(f_m-sign(f_m)*f_trigger)*P
        t_e_ref = (P + ΔP)/ω_r
        #ΔP =-k_P*ω*P
        #dP_inertia = (f_trigger-f_m)*0.1*P
        if abs(ΔP)>(ΔP_max)*P
            ΔP=-sign(f_m)*(ΔP_max)*P
            t_e_ref = (P + ΔP)/ω_r
        end
    else
        ΔP=0
        t_e_ref = P/ω_r
    end
    println("ΔP:",ΔP)

    #u_dq = u*exp(-1im*θ)
    #dθ = v_d*K_PLL
    #ω = dθ_PLL
    #println("ω: ",ω)

    s_e = u*conj(i)# s=(v_d*i_d+v_q*i_q)+j(v_q*i_d-v_d*i_q)
    p_e = real(s_e)
    q_e = imag(s_e)

    # PI speed control:
    de_IP,t_e=PI_control(e_IP,ω_rref-ω_r,K_r1,K_r2)
    #println("t_e:",t_e)
    #e_P = (ω_rref-ω_r)
    #de_IP = e_P
    #println("e_IP: ",e_IP)
    #t_e = K_r1*e_P+K_r2*e_IP
    #p_in +=ΔP
    p_in =(t_e+Δt)/ω_r
    #println("(ω_rref-ω_r): ",(ω_rref-ω_r))
    #println("e_IP: ",e_IP)
    #println("ω_r:",ω_r)
    #println(p_e," ",abs(u)," ",abs(i),";")
    #println("i: ",i)

    t_m = P/ω_r
    dω_r = 1/J*(t_m-t_e)#-D*ω_r)
    #println("t_m-t_e: ",t_m-t_e)
    #println("p_e: ",p_e)
    #println("p_in: ",p_in)

    du_dc = 1/C*(p_in-p_e)
    #println("p_e: ",p_e)

    u_t=abs(u)#TODO: discuss!!!

    # PI voltage control:
    de_IV,i_d=PI_control(e_IV,(u_dcref-u_dc),K_g1,K_g2)
    #e_V = (u_dcref-u_dc)
    #de_IV = e_V
    #i_d = K_g1*e_V + K_g2*de_IV
    #println("u_tref: ",u_tref)
    #println("u_dc: ",u_dc)

    du_tref = K_Q*(Q_ref-q_e)
    di_q = K_v*(u_tref-u_t)
    #println("v_d: ",v_d)
    #println("v_q: ",v_q)
    #println("abs(u)",abs(u))

    du = i - (i_d+1im*i_q)*exp(1im*θ_PLL)
end

export WindTurbineGenType4_inertia
