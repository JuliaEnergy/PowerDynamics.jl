@doc doc"""
```Julia
WindTurbineGenType4(;I_n,k_PLL,f,f_s,T_m,k_P,τ_ω)
```

```
"""

@DynamicNode WindTurbineGenType4(ΔP_max,K_P,K_PLL,P) begin
    MassMatrix(m_u = false,m_int = [true,true,true,true])
end  begin
    @assert K_PLL>=0
    @assert K_P>=0
end [[θ_PLL,dθ_PLL],[ω,dω],[e_Idθ,de_Idθ],[i_d,di_d]] begin
    u_dq = u*exp(-1im*θ_PLL)
    v_d = real(u_dq)
    v_q = imag(u_dq)
    u_dc_ref=1

    de_Idθ,ω_PLL=PIControl(e_Idθ,v_q,K_PLL,1)
    dθ_PLL=ω_PLL
    dω = (-ω + ω_PLL)

    f_m = ω_PLL*2π

    # speed control:

    ω_r=1

    f_trigger=0.1
    if abs(f_m)>f_trigger
        ΔP =-(f_m-sign(f_m)*f_trigger)*P
        p_in_ref = (P + ΔP)#/ω_r
        if abs(ΔP)>(ΔP_max)*P
            ΔP=-sign(f_m)*(ΔP_max)*P
            p_in_ref = (P + ΔP)
        end
    else
        ΔP=0
        p_in_ref = P
    end

    u_dc_ref=1
    i_q=0
    i_dq = (i_d+1im*i_q)
    s_e = u_dq*conj(i_dq)# s=(v_d*i_d+v_q*i_q)+j(v_q*i_d-v_d*i_q)
    p_e = real(s_e)
    q_e = imag(s_e)
    di_d = p_in_ref/v_d-i_d

    du = i - i_dq*exp(1im*θ_PLL)
end

export WindTurbineGenType4
