@DynamicNode GridFormingTecnalia(ω_r,τ_U, τ_I, τ_P, τ_Q, n_P, n_Q, K_P, K_Q, P, Q, V_r, R_f, X_f) begin
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
end [[u_fil_r,du_fil_r],[u_fil_i,du_fil_i],[i_fil_r,di_fil_r],[i_fil_i,di_fil_i],[ω, dω]] begin

    du_fil_r = 1/τ_U*(-u_fil_r + real(u))
    du_fil_i = 1/τ_U*(-u_fil_i + imag(u))

    di_fil_r = 1/τ_I*(-i_fil_r + real(i))
    di_fil_i = 1/τ_I*(-i_fil_i + imag(i))

    u_fil = u_fil_r + 1im*u_fil_i
    i_fil = i_fil_r + 1im*i_fil_i

    p = real(u_fil * conj(i_fil))
    dp =  du_fil_r*i_fil_r+u_fil_r*di_fil_r+du_fil_i*i_fil_i+du_fil_i*di_fil_i

    q = imag(u_fil * conj(i_fil))
    dq = -du_fil_r*i_fil_i-u_fil_r*di_fil_i+du_fil_i*i_fil_r+u_fil_i*di_fil_r

    dϕ = ω-ω_r
    dω = 1/τ_P*(ω_r-ω) + K_P/τ_P*(P-p) - K_P/n_P*dp
    v = abs(u) 
    #v = abs(u_fil)
    dv = 1/τ_Q*(V_r-v) + K_Q/τ_Q*(Q-q) - K_Q/n_Q*dq
    du = u * 1im * dϕ + dv*(u/v)
end

export GridFormingTecnalia
