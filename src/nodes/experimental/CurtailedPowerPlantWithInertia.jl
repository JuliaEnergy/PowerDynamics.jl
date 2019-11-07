@doc doc"""
```Julia
CurtailedPowerPlantWithInertia(;)
```
"""
@DynamicNode CurtailedPowerPlantWithInertia(option1,P,K_ω,K_PPLL,K_IPLL,T_d,T_f) begin
    MassMatrix(m_u = false,m_int = [true,true,true,true])
end  begin
    @assert K_PPLL>=0
end [[θ_PLL,dθ_PLL],[e_Idθ,de_Idθ],[z,dz],[y,dy]] begin
    function PI_control(e_I,e,K_P,K_I)
        @assert K_P>=0
        @assert K_I>=0
        de_I=e
        u = K_P*e+K_I*e_I
        return [de_I,u]
    end

    i_dq = i*exp(-1im*θ_PLL)
    u_dq = u*exp(-1im*θ_PLL)
    v_d = real(u_dq)
    v_q = imag(u_dq)
    i_d = real(i_dq)

    #PLL dynamics for frequency measurement
    de_Idθ,ω_PLL=PI_control(e_Idθ,v_q,K_PPLL,K_IPLL)
    dθ_PLL=ω_PLL

    ω = 1/T_d*(ω_PLL*K_ω-z)
    dz = ω

    if option1
        di_d_ref = 1/T_H*(i_d_ref-ω)
        deId,v_d = PI_control(eId,i_d_ref-i_d,K_Pid,K_Iid)
        deIq,u_q = PI_control(eIq,i_q_ref-i_q,K_Pid,K_Iid)
        du = i - i_dq*exp(1im*θ_PLL)
    else
        p = real(u_dq * conj(i_dq))
        P_I = 1/T_f*(ω-y)
        dy = P_I
        P_ref= P - P_I
        println("i_d",i_d)
        println("p",p)
        du = p - P_ref
    end
end

export CurtailedPowerPlantWithInertia
