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
    u = -je_c e^{j\theta} = -j(e_d + je_q)e^{j\theta}\\
    e_c= e_d + je_q = jue^{-j\theta}\\
    i  = -ji'e^{j\theta} = -j(i_d+ j i_q )e^{j\theta} = Y^L \cdot (u) \\
    i_c= i_d + ji_q = jie^{-j\theta}\\
    p = \Re (i^* u)
```
The fourth-order equations read (according to Sauer, p. 140, eqs. (6110)-(6114)) and p. 35 eqs(3.90)-(3.91)
```math
    \frac{d\theta}{dt} = \omega \\
     \frac{d\omega}{dt} = P-D\omega - p -(x'_q-x'_d)i_d i_q\\
    \frac{d e_q}{dt} = \frac{1}{T'_d} (- e_q - (x_d - x'_d) i_{d}+ e_f) \\
    \frac{d e_d}{dt} = \frac{1}{T'_q} (- e_d + (x_q - x'_q) i_{q})  \\
```
With the PowerDynamics.jl \time{naming conventions} of $i$ and $u$ they read as
```math
   \dot u = \frac{d}{dt}(-j e_c e^{j\theta})=-j(\dot e_d + j\dot e_q)e^{j\theta} + uj\omega
```
"""
@DynamicNode FourthEq(H, P, D, Ω, E_f, T_d_dash ,T_q_dash ,X_q_dash ,X_d_dash,X_d, X_q) <: OrdinaryNodeDynamics() begin
    @assert H > 0 "inertia (H) should be >0"
    @assert D > 0 "damping (D) should be >0"
    @assert T_d_dash > 0 "time constant of d-axis (T_d_dash) should be >0"
    @assert T_q_dash > 0 "time constant of q-axis (T_q_dash) should be >0"
    @assert X_d_dash >= 0 "transient reactance of d-axis (X_d_dash) should be >=0"
    @assert X_q_dash >= 0 "transient reactance of q-axis (X_q_dash) should be >=0"
    @assert X_d >= 0 "reactance of d-axis (X_d_dash) should be >=0"
    @assert X_q >= 0 "reactance of q-axis (X_q_dash) should be >=0"

    Ω_H = Ω * 2pi / H #    norm = 2 * H / (2 * np.pi * 50)  # normalize the parameters as done for coupling_const, input_power, damping_const

end [[θ,dθ],[ω, dω]] begin
    i_c = 1im*i*exp(-1im*θ)
    e_c = 1im*u*exp(-1im*θ)
    p = real(u * conj(i))
    e_d = real(e_c)
    e_q = imag(e_c)
    i_d = real(i_c)
    i_q = imag(i_c)

    dθ = ω
    de_d = (1 / T_q_dash)* (- e_d + (X_q - X_q_dash)* i_q)
    de_q = (1 / T_d_dash)* (- e_q - (X_d - X_d_dash) * i_d + E_f)
    de_c = de_d + 1im*de_q
    du = -1im*de_c*exp(1im*θ)+ u*1im*ω
    dω = (P - D*ω - p- (X_q_dash - X_d_dash)*i_d* i_q)*Ω_H
end

export FourthEq
