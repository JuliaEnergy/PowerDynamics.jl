# (C) 2018 Potsdam Institute for Climate Impact Research, authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)

@doc doc"""
```Julia
SwingEq(;H, P, D, Ω)
```

A node type that applies the swing equation to the frequency/angle dynamics and
keeps the voltage magnitude as is.

Additionally to ``u``, it has the internal dynamic variable ``\omega`` representing
the frequency of the rotator relative to the grid frequency ``\Omega``, i.e. the
real frequency ``\omega_r`` of the rotator is given as ``\omega_r = \Omega + \omega``.

# Keyword Arguments
- `H`: inertia
- `P`: active (real) power output
- `D`: damping coefficient
- `Ω`: rated frequency of the power grid, often 50Hz
- `T_d_dash`: time constant of d-axis
- `T_q_dash`: time constant of q-axis
- `X_d_dash`: transient reactance of d-axis
- `X_q_dash`: transient reactance of q-axis
- `X_d`: reactance of d-axis
- `X_d`: reactance of q-axis

# Mathematical Representation
Using `FourthEq` for node ``a`` applies the equations
```math
\frac{du_a}{dt} = i u_a  \omega_a, \\
\frac{H}{2\pi\Omega}\frac{d\omega_a}{dt} = P_a - D_a\omega_a - \Re\left(u_a \cdot i_a^*\right),
\frac{du}{dt} = \frac{1}{T_q_dash} \left(- \Re(u) + Ω_H(X_q - X_q_dash) \Im(i_c)\right)
+ j \frac{1}{T_d_dash} \left(- \Im(u) - Ω_H(X_d - X_d_dash) \cdot \Re(i_c) + E_f\right)

```
"""
@DynamicNode FourthEq(H, P, D, Ω, E_f, T_d_dash ,T_q_dash ,X_q_dash ,X_d_dash,X_d, X_q) <: OrdinaryNodeDynamics() begin
    @assert H > 0 "inertia (H) should be >0"
    @assert D > 0 "damping (D) should be >0"
    @assert T_d_dash > 0 "time constant of d-axis (T_d_dash) should be >0"
    @assert T_q_dash > 0 "time constant of q-axis (T_q_dash) should be >0"
    @assert X_d_dash > 0 "transient reactance of d-axis (X_d_dash) should be >0"
    @assert X_q_dash > 0 "transient reactance of q-axis (X_q_dash) should be >0"
    @assert X_d > 0 "reactance of d-axis (X_d_dash) should be >0"
    @assert X_q > 0 "reactance of q-axis (X_q_dash) should be >0"

    Ω_H = Ω * 2pi / H #    norm = 2 * H / (2 * np.pi * 50)  # normalize the parameters as done for coupling_const, input_power, damping_const

end [[ω, dω]] begin
    p = real(u * conj(i))

    dϕ = ω # dϕ is only a temp variable that Julia should optimize out
    #du = u * im * dϕ
    du_d = (1 / T_q_dash)* (- real(u) + Ω_H*(X_q - X_q_dash)* imag(i))
    du_q = (1 / T_d_dash)* (- imag(u) - Ω_H*(X_d - X_d_dash) * real(i) + E_f)
    du = du_d + 1im*du_q
    dω = (P - D*ω - p)*Ω_H
end

export FourthEq
