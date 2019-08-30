@doc doc"""
```Julia
PVInverterWithFrequencyControl(I_n,k_PLL,f,f_s,T_m,k_P)

TODO docs
```

"""
@DynamicNode PVInverterWithFrequencyControl(I_n,k_PLL,f,f_s,T_m,k_P) begin
    MassMatrix(m_u = false,m_int = [true,true,true,true])
end  begin
    @assert k_PLL>=0
    @assert f>=0
    @assert f_s>=0
    @assert T_m >=0
    @assert k_P >=0
end [[θ_PLL,dθ_PLL],[v_xm,dv_xm],[v_ym,dv_ym],[P,dP]] begin
    v_x = real(u)
    v_y = imag(u)
    p = real(u * conj(i))

    dv_xm = 1/T_m*(v_x-v_xm)
    dv_ym = 1/T_m*(v_y-v_ym)

    println("abs(u): ",abs(u))
    println("v_xm: ",v_xm)
    println("v_ym: ",v_ym)


    v_d = v_xm*cos(θ_PLL)+v_ym*sin(θ_PLL)
    v_q = -v_xm*sin(θ_PLL)+v_ym*cos(θ_PLL)

    println("v_d: ",v_d)
    println("v_q: ",v_q)

    dθ_PLL = -v_q*(k_PLL)
    if isnan(dθ_PLL)
        stop()
    end


    f_m = (1+dθ_PLL)*f
    ω = dθ_PLL*f*2π
    println("P: ",P)
    println("p: ",p)
    #P=p

    println("f_m: ",f_m)
    dP =-k_P*(f_m-f_s)/f*(v_d*I_n)
    I_P = P/v_d
    println("I_P: ",I_P)
    I_Q = 0.
    i_x = I_P*cos(θ_PLL) - I_Q*sin(θ_PLL)
    i_y = I_P*sin(θ_PLL) + I_Q*cos(θ_PLL)
    println("abs(i_x+1im*i_y): ",i_x+1im*i_y)
    println("i: ",i)
    du = i-(i_x+1im*i_y)
    #du = u*1im*ω
end

export PVInverterWithFrequencyControl
