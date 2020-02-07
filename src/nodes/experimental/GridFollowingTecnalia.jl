@DynamicNode GridFollowingTecnalia(τ_u,ω_ini,K_pω,K_iω,K_ω,K_v,ω_r,V_r,P,Q_r) begin
    MassMatrix(m_u = false, m_int = [true,true,true,true])
end begin
    @assert τ_u >= 0
    @assert ω_ini >= 0
    @assert K_pω >= 0
    @assert K_iω >= 0
    @assert K_ω >= 0
    @assert K_v >= 0
end [[u_fil_d,du_fil_d],[u_fil_q,du_fil_q],[e_Iω,de_Iω],[θ,dθ]] begin

    u_dq = exp(-im*θ)*u
    du_fil_d = 1/τ_u*(-u_fil_d + real(u_dq))
    du_fil_q = 1/τ_u*(-u_fil_q + imag(u_dq))
    u_fil_dq = u_fil_d + im*u_fil_q

    e_ω = -u_fil_q
    de_Iω= e_ω
    ω = ω_ini - K_pω*e_ω - K_iω*e_Iω
    dθ = ω

    p = -K_ω*(ω-ω_r) + P
    q = -K_v*(abs(u_fil_dq)-V_r) + Q_r
    i_fil_dq = (p-im*q)/u_fil_dq
    du = i - i_fil_dq*exp(im*θ)
end

export GridFollowingTecnalia
