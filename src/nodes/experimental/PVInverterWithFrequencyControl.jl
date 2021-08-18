@doc doc"""
```Julia
PVInverterWithFrequencyControl(;I_n,k_PLL,f,f_s,T_m,k_P,τ_ω)
```

This implementation of a generic inverter model is built with standard components according to the report on "Modelling of Inverter-Based
Generation for Power System Dynamic Studies" of the joint CIGRE working group.

Additionally to ``u``, the global network-side voltage, it has the internal dynamic variables:
* `θ_PLL`: phase determined by the PLL
* `v_xm`: x-component of measured grid-side voltage
* `v_ym`: y-component of measured grid-side voltage
* `P`: power infeed of power plant
* `ω`: frequeny deviation in rad/s.

# Keyword arguments are:
* `I_n`: the nominal current of the PV plant
* `k_PLL`: the PLL constant,
* `f`: the nominal frequency (50Hz usually)
* `f_s`: the set point at which droop control is triggered
* `T_m`: the time constant of the low pass filter for measuring the voltage locally at the inverter
* `k_P`: the droop control constant
* `τ_ω` time constant of the frequeny filter

# Mathematical Representation
Using `PVInverterWithFrequencyControl` applies the equations
```math
    v_x = \Re(u)\\
    v_y = \Im(u)
```

The network-side power is
```math
p = \Re(u \cdot i*)
```

The nominal current of the PV inverter, ``I_n``, is completely active current (we do not consider voltage regulation so far):
```math
    I_{P,max} = I_n \\
    I_{Q,max} = 0
```

```math
    \frac{dv_{xm}}{dt} = 1/T_m(v_x-v_{xm})\\
    \frac{dv_{ym}}{dt} = 1/T_m(v_y-v_{ym})
```

Since PowerDynamics is working with phasor units, this model has two d-q-systems (for the power plant and for the grid), $I_p/v_d$ and $I_q/v_q$ are the local coordinates and $i_x/v_x$ and $i_y/v_y$ are the global (grid) coordinates.
```math
    u_dq = (v_{xm}+j v_{ym})e^{-j\theta_{PLL}}
    v_d = \Re(u_dq)\\
    v_q = \Im(u_dq)
```

The local coordinates are chosen such that ``v_q=0``:
```math
    \dot{\theta}_{PLL} = v_q(k_{PLL})
```
The frequency deviation, ``\dot{\theta}_{PLL}``,  is obtained thanks to the PLL controller of the units. Therefore, the measured frequency in Hz is given by
```math
        f_m = (1+\dot{\theta}_{PLL})f.
```

An additional filter is added to avoid a too fast PV reaction leading to unwanted oscillations of the active current during a short-term fault:
```math
    \frac{d\omega}{dt} = 1/\tau_\omega(-\omega + \dot{\theta}_{PLL}\cdot2\pi f)
```
With the active power determined by the droop control the current `I_P`equals:
```math
    I_P=P/v_d
```

Implementing the frequency dead band for overfrequnency control with:
```math
    \frac{dP}{dt} =-k_P \cdot d\omega \cdot P_{ext} \ , \text{if ($f_m>f_s)$}\\
    \frac{dP}{dt}  = 0 \ , \text{else}
```
where ``f_s`` is the frequency at which the active power output starts decreasing and ``k_P`` is the droop control constant (in percentage of the rated power ``P_{ext}``):
```math
    P_{ext} = \Re(I_{P,max}\cdot(v_d+jv_q))= I_N v_d.
```

```math
    0\cdot \dot{u} =i - (I_P+ j I_Q)e^{j\theta_{PLL}}
```
"""
@DynamicNode PVInverterWithFrequencyControl(I_n,k_PLL,f,f_s,T_m,k_P,τ_ω) begin
    MassMatrix(m_u = false,m_int = [true,true,true,true,true])
end  begin
    @assert k_PLL>=0
    @assert f>=0
    @assert f_s>=0
    @assert T_m >=0
    @assert k_P >=0
    @assert τ_ω >0
end [[θ_PLL,dθ_PLL],[v_xm,dv_xm],[v_ym,dv_ym],[P,dP],[ω,dω]] begin
    p = real(u * conj(i))

    v_x = real(u)
    v_y = imag(u)

    dv_xm = 1/T_m*(v_x-v_xm)
    dv_ym = 1/T_m*(v_y-v_ym)

    u_dq = (v_xm+1im*v_ym)*exp(-1im*θ_PLL)
    v_d = real(u_dq)
    v_q = imag(u_dq)
    #v_d = v_xm*cos(θ_PLL)+v_ym*sin(θ_PLL)
    #v_q = v_xm*sin(θ_PLL)-v_ym*cos(θ_PLL)

    dθ_PLL = v_q*(k_PLL)

    f_m = (1+dθ_PLL)*f

    # adding filter to avoid a too fast PV reaction leading to unwanted oscillations of the active current during a short
    # term fault since PV frequency regulation are generally expected to occur within seconds.
    dω = 1/τ_ω*(-ω + dθ_PLL*f*2*π)

    if (f_m>f_s)
        dP = k_P*dω*(v_d*I_n)
        #dP =k_P*(f_m-f_s)/f*(v_d*I_n)
    else
        dP=0
    end
    #if (P<v_d*I_n)
    I_P = P/v_d
    #else
    #    I_P = I_n
    #end
    I_Q = 0.
    du = i-(I_P+1im*I_Q)*exp(1im*θ_PLL)
end

export PVInverterWithFrequencyControl
