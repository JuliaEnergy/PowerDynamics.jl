@doc doc"""
```Julia
VSIMinimal(;τ_P,τ_Q,K_P,K_Q,E_r,P,Q)
```

A node type that applies the frequency and voltage droop control to control the frequency and
voltage dynamics. The following implementation is taken from
J. Schiffer et. al. , Automatica 50 (2014) 2457–2469.

Additionally to ``u``, it has the internal dynamic variable ``\omega`` representing
the frequency of the rotator relative to the grid frequency ``\Omega``, i.e. the
real frequency ``\omega_r`` of the rotator is given as ``\omega_r = \Omega + \omega``.

# Keyword Arguments
- `τ_p`: time constant active power measurement
- `τ_Q`: time constant reactive power measurement
- `K_P`: droop constant frequency droop
- `K_Q`: droop constant voltage droop
- `V_r`: reference/ desired voltage
- `P`: active (real) power infeed
- `Q`: reactive (imag) power infeed


# Mathematical Representation
Using `VSIMinimal` for node ``a`` applies the equations
```math
\dot{\phi}_a=\omega_a\\
 \dot{\omega}_a=\frac{1}{\tau_{P,a}}[-\omega_a-K_{P,a} (\Re\left(u_a \cdot i_a^*\right)-P_{ref,a})]\\
\tau_Q\dot{v}_a=-v_a+V_{ref}-K_{Q,a} (\Im\left(u_a \cdot i_a^*\right)-Q_{ref,a})\\
 \dot{u}_a=\dot{v_a}e^{j\phi}+j\omega_a u_a
```
```
"""
@DynamicNode VSIMinimal(τ_P,τ_Q,K_P,K_Q,V_r,P,Q) <: OrdinaryNodeDynamics() begin
    @assert τ_P > 0 "time constant active power measurement should be >0"
    @assert τ_Q > 0 "time constant reactive power measurement should be >0"
    @assert K_Q > 0 "reactive power droop constant should be >0"
    @assert K_P > 0 "active power droop constant reactive power measurement should be >0"
end [[ω, dω]] begin
    p = real(u * conj(i))
    q = imag(u * conj(i))
    dϕ = ω
    v = abs(u)
    dv = 1/τ_Q*(-v + V_r- K_Q *(q-Q))
    du = u * 1im * dϕ + dv*(u/v)
    dω = 1/τ_P*(-ω-K_P*(p-P))
end

export VSIMinimal
