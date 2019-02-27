@doc doc"""
```Julia
VSIVoltagePT1(;τ_v,τ_P,τ_Q,K_P,K_Q,E_r,P,Q)
```

A node type that applies the frequency and voltage droop control to control the frequency and
voltage dynamics.

`VSIVoltagePT1` models an inverters as AC voltage source which means the amplitude and
frequency can defined by the designer (often called grid-forming inverter mode).
The frequency regulation is assumed to be instantaneous, but the voltage control happens with a delay
`τ_v` that is represented by a first order filter.
In addition simple proportional controllers are implemented for frequency and voltage such that the frequency `ω` and
voltage amplitudes `v` of the inverters are modified depending on the deviations (with respect to a desired value) of the active and reactive powers, respectively.
ift is assumed that active and reactive power are measured via low pass fileters with time constant `τ_P` and `τ_Q`, respectively.

Hence, additionally to `u`, it has the internal dynamic variables
* `ω` representing the frequency of the rotator relative to the grid frequency ``Ω``, i.e. the real frequency ``ω_r`` of the rotator is given as ``ω_r = Ω + ω``.
* `q_m` is the measured reactive power at the grid connection point.

# Keyword Arguments
- `τ_v`: time constant voltage control delay
- `τ_P`: time constant active power measurement
- `τ_Q`: time constant reactive power measurement
- `K_P`: droop constant frequency droop
- `K_Q`: droop constant voltage droop
- `V_r`: reference/ desired voltage
- `P`: active (real) power infeed
- `Q`: reactive (imag) power infeed


# Mathematical Representation
Using `VSIVoltagePT1` for node ``a`` (according to J. Schiffer et. al., eq. (6)) gives the equations
```math
\dot{\phi}_a=\omega_a\\
 \dot{\omega}_a=\frac{1}{\tau_{P,a}}[-\omega_a-K_{P,a} (\Re\left(u_a \cdot i_a^*\right)-P_{ref,a})]\\
 \tau_v\dot{v}_{a}=-v_a+V_{ref}-K_{Q,a}(q_{m,a}-Q_{ref,a})\\
 \tau_Q \dot{q}_{m,a}=-q_{m,a}+\Im\left(u_a \cdot i_a^*\right)\\
 \dot{u}_a=\dot{v_a}e^{j\phi}+j\omega_a u_a\\
```
In general ``τ_V ≪ τ_P``, assuming ``τ_V = 0`` would then lead to `VSIMinimal`.

"""
@DynamicNode VSIVoltagePT1(τ_v,τ_P,τ_Q,K_P,K_Q,V_r,P,Q) <: OrdinaryNodeDynamics() begin
    @assert τ_v > 0 "time constant voltage droop delay should be >0"
    @assert τ_P > 0 "time constant active power measurement should be >0"
    @assert τ_Q > 0 "time constant reactive power measurement should be >0"
    @assert K_Q > 0 "reactive power droop constant should be >0"
    @assert K_P > 0 "active power droop constant reactive power measurement should be >0"
end [[ω, dω],[q_m,dq_m]] begin
    p = real(u * conj(i))
    q = imag(u * conj(i))
    dϕ = ω
    v = abs(u)
    dv = 1/τ_v*(-v+V_r - K_Q*(q_m-Q))
    dq_m = 1/τ_Q*(q-q_m)
    du = u * 1im * dϕ + dv*(u/v)
    dω = 1/τ_P*(-ω-K_P*(p-P))
end

export VSIVoltagePT1
