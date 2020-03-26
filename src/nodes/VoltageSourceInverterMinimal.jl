@doc doc"""
```Julia
VSIMinimal(;τ_P,τ_Q,K_P,K_Q,E_r,P,Q)
```

A node type that applies the frequency and voltage droop control to control the frequency and
voltage dynamics. Implemented according to  Schiffer et. al., "Conditions for stability of
droop-controlled inverter-based microgrids", Automatica, 2014.

`VSIMinimal` models an inverters as AC voltage source which means the amplitude and
frequency can defined by the designer (often called grid-forming inverter mode).
The frequency and voltage regulation is assumed to be instantaneous.
In addition simple proportional controllers are implemented for frequency and voltage such that the frequency `ω` and
voltage amplitudes `v` of the inverters are modified depending on the deviations (with respect to a desired value) of the active and reactive powers, respectively.
it is assumed that active and reactive power are measured via low pass filters with time constant `τ_P` and `τ_Q`, respectively.
`VSIMinimal` can be derived from `VSIVoltagePT1` by assuming an instantaneous voltage regulation without delay.

Additionally to ``u``, it has the internal dynamic variable ``ω`` representing
the frequency of the inverter frequency relative to the grid frequency ``Ω=2π50``Hz, i.e. the
real frequency ``ω_r`` of the inverter frequency is given as ``ω_r = Ω + ω``.


# Keyword Arguments
- `τ_p`: time constant active power measurement
- `τ_Q`: time constant reactive power measurement
- `K_P`: droop constant frequency droop
- `K_Q`: droop constant voltage droop
- `V_r`: reference/ desired voltage
- `P`: active (real) power infeed
- `Q`: reactive (imag) power infeed


# Mathematical Representation
Using `VSIMinimal` for node ``a`` (according to J. Schiffer et. al., eq. (7)) gives the equations
```math
\dot{\phi}_a=\omega_a\\
 \dot{\omega}_a=\frac{1}{\tau_{P,a}}[-\omega_a-K_{P,a} (\Re\left(u_a \cdot i_a^*\right)-P_{a})]\\
\tau_Q\dot{v}_a=-v_a+V_{r}-K_{Q,a} (\Im\left(u_a \cdot i_a^*\right)-Q_{a})\\
 \dot{u}_a=\dot{v_a}e^{j\phi}+j\omega_a u_a
```
```
"""
@DynamicNode VSIMinimal(τ_P,τ_Q,K_P,K_Q,V_r,P,Q) begin
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
