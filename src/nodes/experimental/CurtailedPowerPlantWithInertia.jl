@doc doc"""
```Julia
CurtailedPowerPlantWithInertia(;)
```
"""
@DynamicNode CurtailedPowerPlantWithInertia(P,ω_0,T_AI,K_PPLL,K_IPLL,T_d,T_f,K_PV,K_IV) begin
    MassMatrix(m_u = false,m_int = [true,true,true,true,true,true])
end  begin
    @assert K_PPLL>=0
    @assert K_IPLL>=0
    @assert K_PV>=0
    @assert K_IV>=0
    @assert T_d>=0
    @assert T_f>=0
    @assert T_AI>=0
end [[θ_PLL,dθ_PLL],[e_Iω,de_Iω],[ω,dω],[y,dy],[e_Iiq,de_Iiq],[e_Iid,de_Iid]] begin
    # transformation from global to local dq-coordinates
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

    # inertia emulation with df/dt response
    dω = 1/T_d*(ω_PLL-ω) # low pass filter
    i_dI = 1/(T_f*ω_0)*(-ω-y*ω_0) # high pass filter
    dy = i_dI
    i_d_ref= P/v_d+T_AI*i_dI

    # current controller
    de_Iid,v_d_ref = PIControl(e_Iid,i_d_ref-i_d,K_PV,K_IV)
    de_Iiq,v_q_ref = PIControl(e_Iiq,i_q_ref-i_q,K_PV,K_IV)

    du = u - (v_d_ref+1im*v_q_ref)*exp(1im*θ_PLL)
end

export CurtailedPowerPlantWithInertia
