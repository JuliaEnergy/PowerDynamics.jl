@doc doc"""
```Julia
CurtailedPowerPlantWithInertia(;)
```
"""
@DynamicNode CurtailedPowerPlantWithInertia(PICFunction,P,ω_0,T_AI,K_PPLL,K_IPLL,T_d,T_f,K_PV,K_IV) begin
    MassMatrix(m_u = false,m_int = [true,true,true,true,true,true])#true,true,true,true
end  begin
    @assert K_PPLL>=0
end [[θ_PLL,dθ_PLL],[e_Iω,de_Iω],[ω,dω],[y,dy],[e_Iiq,de_Iiq],[e_Iid,de_Iid]] begin#,,[i_d,di_d],[i_q,di_q],[v_q,dv_q],[v_d,dv_d]
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
    i_q = imag(i_dq)

    i_q_ref=0
    #PLL dynamics for frequency measurement
    e_ω= -v_q
    de_Iω= e_ω
    ω_PLL=K_PPLL*e_ω+K_IPLL*e_Iω
    dθ_PLL=ω_PLL


    dω = 1/T_d*(ω_PLL-ω)

    i_dI = 1/(T_f*ω_0)*(-ω-y*ω_0)
    dy = i_dI
    i_d_ref= P/v_d+T_AI*i_dI#(P-i_q*v_q)/v_d+i_dI

    # current controller
    if PICFunction
        de_Iid,v_d_ref = PIControl(e_Iid,i_d_ref-i_d,K_PV,K_IV)
        de_Iiq,v_q_ref = PIControl(e_Iiq,i_q_ref-i_q,K_PV,K_IV)
    else
        e_id = i_d_ref-i_d
        de_Iid= e_id
        v_d_ref=K_PV*e_id+K_IV*e_Iid

        e_iq = i_q_ref-i_q
        de_Iiq= e_iq
        v_q_ref=K_PV*e_iq+K_IV*e_Iiq
    end

    du = u - (v_d_ref+1im*v_q_ref)*exp(1im*θ_PLL)
    #0 = i - (i_d_ref+1im*i_q_ref)*exp(1im*θ_PLL)
end

export CurtailedPowerPlantWithInertia
