@doc doc"""
```Julia
VISMA(;P_m,Ω,u_e,R_e,L_e,R_d,R_q,L_d,L_q,R_D,R_Q,L_D,L_Q,M_Dd,M_Qq,M_ed,M_eD,Z_p,H)
```

This implementation of a Virtual Synchronous Machine is undertaken according to TU Clausthal.

Additionally to ``u``, the global network-side voltage, it has the internal dynamic variables:
* `θ`:
* `Ψ_d`:
* `Ψ_q`:
* `Ψ_D`:
* `Ψ_Q`:
* `Ψ_e`:
* `ω_r`: frequeny deviation in rad/s.

# Keyword arguments are:
* `R_e`:
* `L_e`:
* `R_d`:
* `R_q`:
* `L_d`:
* `L_q`:
* `R_D`:
* `R_Q`:
* `L_D`:
* `L_Q`:
* `M_{Dd}`:
* `M_{Qq}`:
* `M_{ed}`:
* `M_{eD}`:
* `Z_{p}`:
* `J`:

# Mathematical Representation
Using `VISMA` applies the equations for the stator circuit:
```math
    i_d = 1/R_d*(u_d-\frac{dΨ_d}{dt}+ω*Ψ_q)\\
    i_q = 1/R_q*(u_q+\frac{dΨ_q}{dt}-ω*Ψ_d)
```
equal to:
```math
    \frac{dΨ_d}{dt} =  u_d - i_d*R_d +ω*Ψ_q\\
    \frac{dΨ_q}{dt} = -u_q + i_q*R_q +ω*Ψ_d
```

where
```math
    Ψ_d = L_d*i_d+M_{Dd}*i_D+M_{ed}*i_e\\
    Ψ_q = L_q*i_q+M_{Qq}*i_Q
```
The local d-q-corrdinates are transformed from

```math
    (u_d,u_q) = T*(u_a,u_b,u_c)\\
    (i_a,i_b,i_c) = T^{-1}*(i_d,i_q)
```
Exciter- and damper:
```math
    i_e = 1/R_e*(u_e-\frac{dΨ_e}{dt})\\
    Ψ_e = L_e*i_e+ M_{ed}*i_d+ M_{eD}*i_D\\
    i_D = -\frac{1}{R_D}*\frac{dΨ_D}{dt}\\
    i_Q = -\frac{1}{R_Q}*frac{dΨ_Q}{dt}\\
    Ψ_D = L_D*i_D+M_{Dd}*i_d+M_{eD}*i_e\\
    Ψ_Q = L_Q*i_Q+M_{Qq}*i_q
```
equal to:
```math
    \frac{dΨ_e}{dt}=u_e-i_e*R_e\\
    Ψ_e = L_e*i_e+ M_{ed}*i_d+ M_{eD}*i_D\\
    \frac{dΨ_D}{dt}= -i_D*R_D\\
    \frac{dΨ_Q}{dt}= -i_Q*R_Q\\
    Ψ_D = L_D*i_D+M_{Dd}*i_d+M_{eD}*i_e\\
    Ψ_Q = L_Q*i_Q+M_{Qq}*i_q
```

Rotor Mechanics:
```math
    P_{el} = \frac{3}{2}*Z_p*(Ψ_d*i_q-Ψ_q*i_d)*Ω \\
    \frac{dω}{dt}=\frac{1}{H}(P_{el}_P_m)\\
    \frac{dθ}{dt}=ω
```

"""
@DynamicNode VISMA(P_m,Ω,u_e,R_e,L_e,R_d,R_q,L_d,L_q,R_D,R_Q,L_D,L_Q,M_Dd,M_Qq,M_ed,M_eD,Z_p,H)   begin
    MassMatrix(m_u = false,m_int=[true,true,true,true,true,true,true])
end   begin
    @assert Ω>0
    @assert R_e>=0
    @assert L_e>=0
    @assert R_d >=0
    @assert R_q >=0
    @assert L_d >=0
    @assert L_q >=0
    @assert R_D >=0
    @assert R_Q >=0
    @assert L_D >=0
    @assert L_Q >=0
    @assert M_Dd>=0
    @assert M_Qq>=0
    @assert M_ed>=0
    @assert M_eD>=0
    @assert Z_p>0
    @assert H>0
end [[θ,dθ],[Ψ_d,dΨ_d],[Ψ_q,dΨ_q],[Ψ_D,dΨ_D],[Ψ_Q,dΨ_Q],[Ψ_e,dΨ_e],[ω,dω]] begin
    e_c = 1im*u*exp(-1im*θ)
    #p = real(u * conj(i))
    u_d = real(e_c)
    u_q = imag(e_c)

    # introducing helper parameters
    δ = L_e+M_eD^2/L_D
    α = L_d*L_D-M_Dd^2
    β = (M_eD*M_Dd-L_D*M_ed)/δ
    γ = -M_ed-M_eD*M_Dd/L_D

    i_d = (Ψ_d*L_D-Ψ_D*M_Dd+(Ψ_e-M_eD/L_D*Ψ_D)*β)/(α-γ*β)
    i_q=(Ψ_q*L_Q-Ψ_Q*M_Qq)/(L_q*L_Q -M_Qq^2)

    i_e = (Ψ_e - M_ed*i_d- M_eD/L_D*Ψ_D-M_eD*M_Dd/L_D*i_d)/δ

    i_D=1/L_D*(Ψ_D -M_Dd*i_d-M_eD*i_e)
    i_Q = 1/L_Q*(Ψ_Q -M_Qq*i_q)

    dΨ_d =  u_d - i_d*R_d +ω*Ψ_q
    dΨ_q = -u_q + i_q*R_q +ω*Ψ_d

    Ψ_d = L_d*i_d+M_Dd*i_D+M_ed*i_e
    Ψ_q = L_q*i_q+M_Qq*i_Q

    dΨ_e=u_e-i_e*R_e
    #Ψ_e = L_e*i_e+ M_ed*i_d+ M_eD*i_D
    #i_e =1/L_e*(Ψ_e - M_ed*i_d- M_eD*i_D)

    dΨ_D= -i_D*R_D
    dΨ_Q= -i_Q*R_Q

    #Ψ_D = L_D*i_D+M_Dd*i_d+M_eD*i_e
    #Ψ_Q = L_Q*i_Q+M_Qq*i_q

    P_el = 3/2*Z_p*(Ψ_d*i_q-Ψ_q*i_d)
    dω=Ω/H*(P_el-P_m)
    dθ=ω
    println(u)
    #du = -1im*de_c*exp(1im*θ)+ u*1im*ω
    du = (i_d+1im*i_q)*exp(1im*θ)- 1im*i
end

export VISMA
