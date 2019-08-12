@doc doc"""
```Julia
VSIMinimal(;τ_P,τ_Q,K_P,K_Q,E_r,P,Q)
```

A node type that applies the frequency and voltage droop control to control the frequency and
voltage dynamics.

`VSIMinimal` models an inverters as AC voltage source which means the amplitude and
frequency can defined by the designer (often called grid-forming inverter mode).
The frequency and voltage regulation is assumed to be instantaneous.
In addition simple proportional controllers are implemented for frequency and voltage such that the frequency `ω` and
voltage amplitudes `v` of the inverters are modified depending on the deviations (with respect to a desired value) of the active and reactive powers, respectively.
ift is assumed that active and reactive power are measured via low pass fileters with time constant `τ_P` and `τ_Q`, respectively.
`VSIMinimal` can be derived from `VSIVoltagePT1` by assuming an instantaneous voltage regulation without delay.

Additionally to ``u``, it has the internal dynamic variable ``ω`` representing
the frequency of the rotator relative to the grid frequency ``Ω``, i.e. the
real frequency ``ω_r`` of the rotator is given as ``ω_r = Ω + ω``.


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
@DynamicNode PVInverterWithFrequencyControl(I_n,k_PLL,f,f_s,T_m,k_P) begin
    MassMatrix()
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

    v_q = v_xm*cos(θ_PLL)+v_ym*sin(θ_PLL)
    v_d = -v_xm*sin(θ_PLL)+v_ym*cos(θ_PLL)

    dθ_PLL = v_q*(-k_PLL)

    f_m = (1+dθ_PLL)*f
    P = I_P*v_d
    dP = -k_P*(f_m-f_s)/f*p

    i_x = I_P*cos(θ_PLL) + I_Q*sin(θ_PLL)
    i_y = I_P*sin(θ_PLL) - I_Q*cos(θ_PLL)

    du = i-(i_x+1im*i_y)
end

export PVInverterWithFrequencyControl
