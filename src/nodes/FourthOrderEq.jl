# (C) 2018 Potsdam Institute for Climate Impact Research, authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)

@doc doc"""
```Julia
FourthEq(H, P, D, Ω, E_f, T_d_dash ,T_q_dash ,X_q_dash ,X_d_dash,X_d, X_q)
```

A node type that applies the 4th-order synchronous machine model with frequency/angle and voltage dynamics,
which is implemented according to P. Sauer, "Power System Dynamics and Stability".
For an illustration of a synchronous machine schematic see P. Sauer, Fig. 3.1 on p. 25.

Usually the swing equation (`SwingEq`) is used for short time periods to analyze the transient behavior
of generators in a power grid, the so-called first swing. The 4th-order model  also takes the back reaction
of the power flow onto the voltage into account. This has the effect that the angle of the voltage as seen by the power grid,
and the angle of the rotating mass are no longer the same but become dynamically coupled.

Additionally to ``u``, it has the internal dynamic variables
* ``ω`` representing the frequency of the rotator relative to the grid frequency ``Ω``, i.e. the real frequency ``ω_r`` of the rotator is given as ``\omega_r = \Omega + \omega`` and
* ``θ`` representing the relative angle of the rotor with respect to the voltage angle ``ϕ``.

# Keyword Arguments
- `H`: shaft inertia constant (given in [s]), defined according to P. Sauer, p. 33, eq. (3.60)
- `P`: active (real) power output, also called the mechanical torque applied to the shaft, given in [pu]
- `D`: damping coefficient (given in [s], see P. Sauer, eq. (5.156) where the damping torque is equal `Dω`)
- `Ω`: rated frequency of the power grid, often ``2π⋅50Hz``
- `T_d_dash`: time constant of d-axis, given in [s], see P. Sauer, chapter 3.7, p. 54 for a general explanation on time constants
- `T_q_dash`: time constant of q-axis, given in [s]
- `X_d_dash`: transient reactance of d-axis, given in [pu]
- `X_q_dash`: transient reactance of q-axis, given in [pu]
- `X_d`: reactance of d-, given in [pu]
- `X_d`: reactance of q-axis, given in [pu]
- `E_f`: scaled field voltage, which, if set equal to 1.0 pu, gives 1.0 pu open-circuit terminal voltage. The physical device that provides the value of `E_f` is called the exciter (according to P. Sauer, p. 65)

# Mathematical Representation
Using `FourthEq` for node ``a`` applies the equations
```math
    u = -je_c e^{j\theta} = -j(e_d + je_q)e^{j\theta}\\
    e_c= e_d + je_q = jue^{-j\theta}\\
    i  = -ji'e^{j\theta} = -j(i_d+ j i_q )e^{j\theta} = Y^L \cdot u \\
    i_c= i_d + ji_q = jie^{-j\theta}\\
    p = \Re (i^* u)
```
where complex voltage and current are described in a co-rotating frame with axes labeled d and q.

The fourth-order equations read (according to P. Sauer, "Power System Dynamics and Stability", p. 140, eqs. (6110)-(6114)) and p. 35 eqs(3.90)-(3.91)
```math
    \frac{d\theta}{dt} = \omega \\
     \frac{d\omega}{dt} = (P-D\omega - p -(X'_q-X'_d)i_d i_q)Ω_H\\
    \frac{d e_q}{dt} = \frac{1}{T'_d} (- e_q - (X_d - X'_d) i_{d}+ E_f) \\
    \frac{d e_d}{dt} = \frac{1}{T'_q} (- e_d + (X_q - X'_q) i_{q})  \\
```
The equations for frequency and phase represent energy conservation and phase shift.
The dynamic equations for the complex voltage show the relationship between the dynamicy of flux linkages and currents which must reflect a
conservative coupling field.

With the PowerDynamics.jl naming conventions of ``i`` and ``u`` they read as
```math
   \dot u = \frac{d}{dt}(-j e_c e^{j\theta})=-j(\dot e_d + j\dot e_q)e^{j\theta} + uj\omega
```
"""
@DynamicNode FourthOrderEq(H, P, D, Ω, E_f, T_d_dash ,T_q_dash ,X_q_dash ,X_d_dash,X_d, X_q) begin
    @assert H > 0 "inertia (H) should be >0"
    @assert D >= 0 "damping (D) should be >=0"
    @assert T_d_dash > 0 "time constant of d-axis (T_d_dash) should be >0"
    @assert T_q_dash > 0 "time constant of q-axis (T_q_dash) should be >0"
    @assert X_d_dash >= 0 "transient reactance of d-axis (X_d_dash) should be >=0"
    @assert X_q_dash >= 0 "transient reactance of q-axis (X_q_dash) should be >=0"
    @assert X_d >= 0 "reactance of d-axis (X_d_dash) should be >=0"
    @assert X_q >= 0 "reactance of q-axis (X_q_dash) should be >=0"

    Ω_H = Ω / (2*H)
    
end [[θ,dθ],[ω, dω]] begin # u= u_r+1im*u_i = v*exp(1im*ϕ)
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

export FourthOrderEq
