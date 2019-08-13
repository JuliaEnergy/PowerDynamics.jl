@doc doc"""
```Julia
PVInverterWithFrequencyControl(I_n,k_PLL,f,f_s,T_m,k_P)

TODO docs
```

"""
@DynamicNode PVInverterWithFrequencyControl(I_n,k_PLL,f,f_s,T_m,k_P) begin
    MassMatrix(m_u = true,m_int = [true,true,true,false])
end  begin
    # no prep
end [[θ_PLL,dθ_PLL],[v_xm,dv_xm],[v_ym,dv_ym],[P,dP]] begin
    v_x = real(u)
    v_y = imag(u)
    p = real(u * conj(i))
    I_P = I_n
    I_Q = 0

    dv_xm = 1/T_m*(v_x-v_xm)
    dv_ym = 1/T_m*(v_y-v_ym)

    v_d = v_xm*cos(θ_PLL)+v_ym*sin(θ_PLL)
    v_q = -v_xm*sin(θ_PLL)+v_ym*cos(θ_PLL)

    dθ_PLL = -v_q*(k_PLL)
    if isnan(dθ_PLL)
        stop()
    end

    f_m = (1+dθ_PLL)*f
    P = I_P*v_d -k_P*(f_m-f_s)/f*p
    dP = p-P
    #println(P)
    println(v_q)
    i_x = I_P*cos(θ_PLL) + I_Q*sin(θ_PLL)
    i_y = I_P*sin(θ_PLL) - I_Q*cos(θ_PLL)
    du = u*1im*dθ_PLL
end

export PVInverterWithFrequencyControl
