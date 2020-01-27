@DynamicNode GridFormingTecnalia(ω_r,τ_U, τ_I, τ_P, τ_Q, n_P, n_Q, K_P, K_Q, P, Q, V_r, R_f, X_f) begin
    MassMatrix(m_u = false, m_int = [true,true,true,true,true,true,true,true,true])
end begin
    @assert τ_U >= 0
    @assert τ_I >= 0
    @assert τ_P >= 0
    @assert τ_Q >= 0
    @assert n_P >= 0
    @assert n_Q >= 0
    @assert K_P >= 0
    @assert K_Q >= 0
    @assert V_r >= 0
    @assert R_f >= 0
    @assert X_f >= 0
end [[u_fil_r,du_fil_r],[u_fil_i,du_fil_i],[i_fil_r,di_fil_r],[i_fil_i,di_fil_i],[p,dp],[q,dq],[θ,dθ],[ω, dω],[v,dv]] begin
    # TODO: adjust P_ref to compensate for losses and islanding
    u_dq = exp(-1im*θ)*u
    i_dq = exp(-1im*θ)*i

    du_fil_r = 1/τ_U*(-u_fil_r + real(u_dq))
    du_fil_i = 1/τ_U*(-u_fil_i + imag(u_dq))

    di_fil_r = 1/τ_I*(-i_fil_r + real(i_dq))
    di_fil_i = 1/τ_I*(-i_fil_i + imag(i_dq))

    u_fil = u_fil_r +1im*u_fil_i
    i_fil = i_fil_r +1im*i_fil_i
    du_fil = du_fil_r +1im*du_fil_i
    di_fil = di_fil_r +1im*di_fil_i

    # p = real(u_fil * conj(i_fil)) = u_filr*i_filr-u_fili*i_fili
    # q = imag(u_fil * conj(i_fil)) = -u_fil_r*i_fil_i+u_fil_i*i_fil_r
    #dp = real(u_fil * conj(di_fil) + du_fil * conj(i_fil))
    dp =  du_fil_r*i_fil_r+u_fil_r*di_fil_r-du_fil_i*i_fil_i-du_fil_i*di_fil_i

    #dq = imag(u_fil * conj(di_fil) + du_fil * conj(i_fil))
    dq = -du_fil_r*i_fil_i-u_fil_r*di_fil_i+du_fil_i*i_fil_r+u_fil_i*di_fil_r

    dθ = ω-ω_r
    dω = 1/τ_P*(ω_r-ω) + K_P/τ_P*(P-p) - K_P/n_P*dp
    #println(ω)
    #println("p:",p)
    #println("P:",P)
    dv = 1/τ_Q*(V_r - v) + K_Q/τ_Q*(Q-q) - K_Q/n_Q*dq
    #println(dv)

    v_out = v - R_f*(i_dq - i_fil) - 1im*X_f*i_fil
    du = u - v_out*exp(1im*θ)
end

export GridFormingTecnalia
