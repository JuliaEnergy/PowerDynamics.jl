@DynamicNode GridFormingTecnalia(τ_U, τ_I, τ_P, τ_Q, n_P, n_Q, k_P, k_Q, P, Q, V_r, R_f, X_f) begin
    MassMatrix(m_u = false, m_int = [true,true,true,true,true])
end begin
    @assert τ_U >= 0
    @assert τ_I >= 0
    @assert τ_P >= 0
    @assert τ_Q >= 0
    @assert n_P >= 0
    @assert n_Q >= 0
    @assert k_P >= 0
    @assert k_Q >= 0
    @assert V_r >= 0
    @assert R_f >= 0
    @assert X_f >= 0
end [[u_fil,du_fil],[i_fil,di_fil],[p,dp],[q,dq],[θ,dθ],[ω, dω],[v,dv]] begin

    u_dq = exp(-1im*θ)*u
    i_dq = exp(-1im*θ)*i

    du_fil = 1/τ_U*(-u_fil + u_dq)
    di_fil = 1/τ_I*(-i_fil + i_dq)

    # p = real(u_fil * conj(i_fil))
    # q = imag(u_fil * conj(i_fil))
    dp = real(u_fil * conj(di_fil) + du_fil * conj(i_fil))
    dq = imag(u_fil * conj(di_fil) + du_fil * conj(i_fil))

    dθ = ω
    dω = -1/τ_P*ω + k_P/τ_P*(P-p) - k_P/n_P*dp
    dv = 1/τ_Q*(V_r - v) + k_Q/τ_Q*(Q-q)) - k_Q/n_Q*dq

    v_out = v - R_f*(i_dq - i_fil) - 1im*X_f*i_fil
    du = u - v_out*exp(1im*θ)
end
