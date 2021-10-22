# (C) 2018 Potsdam Institute for Climate Impact Research, authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)

@doc doc"""
```Julia
SwingEq(;H, P, D, Ω)
```

A node type that applies the swing equation to the frequency/angle dynamics and
keeps the voltage magnitude as is. In the following, we followed the implementation
of the 2nd-order Synchronous Machine Model according to  Sauer et. al.
"Power system dynamics and stability", 1998.

Additionally to ``u``, it has the internal dynamic variable ``\omega`` representing
the frequency of the rotator relative to the grid frequency ``\Omega``, i.e. the
real frequency ``\omega_r`` of the rotator is given as ``\omega_r = \Omega + \omega``.

# Keyword Arguments
- `H`: inertia constant (given in [s]), defined according to P. Sauer, p. 33, eq. (3.60)
- `P`: active (real) power output, also called the mechanical torque applied to the shaft, given in [pu]
- `D`: damping coefficient, (given in [s], see P. Sauer, eq. (5.156) where the damping torque is equal `Dω`)
- `Ω`: rated frequency in [1/s] of the power grid, often ``2\pi⋅50``Hz

# Mathematical Representation
Using `SwingEq` for node ``a`` applies the equations
```math
\frac{du_a}{dt} = j \omega_a u_a, \\
\frac{2H}{\Omega}\frac{d\omega_a}{dt} = P_a - D_a\omega_a - \Re\left(u_a \cdot i_a^*\right),
```
which is equivalent to
```math
\frac{d\phi_a}{dt} = \omega_a, \\
v = v(t=0) = \text{const.} \\
\frac{2H}{\Omega}\frac{d\omega_a}{dt} = P_a - D_a\omega_a - \Re\left(u_a \cdot i_a^*\right),
```
where ``H = \frac{1}{2}\frac{J\Omega^2}{S_b}`` for a two-pole machine accoding to
Sauer et. al. eq. (3.60) on p. 33. `S_b` is the rated three-phase MVA of the power system.

"""
@DynamicNode SwingEq(H, P, D, Ω) begin
    @assert D >= 0 "damping (D) should be >=0"
    @assert H > 0 "inertia (H) should be >0"
    Ω_H = Ω / (2*H)
end [[ω, dω]] begin
    p = real(u * conj(i))
    dϕ = ω # dϕ is only a temp variable that Julia should optimize out
    du = u * im * dϕ
    dω = (P - D*ω - p)*Ω_H
end

export SwingEq

@doc doc"""
```Julia
SwingEqLVS(;H, P, D, Ω, Γ, V)
```

A node type that applies the swing equation to the frequency/angle dynamics and
has a linear voltage stability (LVS) term.

Additionally to ``u``, it has the internal dynamic variable ``\omega`` representing
the frequency of the rotator relative to the grid frequency ``\Omega``, i.e. the
real frequency ``\omega_r`` of the rotator is given as ``\omega_r = \Omega + \omega``.

# Keyword Arguments
- `H`: inertia constant (given in [s]), defined according to P. Sauer, p. 33, eq. (3.60)
- `P`: active (real) power output, also called the mechanical torque applied to the shaft, given in [pu]
- `D`: damping coefficient, (given in [s], see P. Sauer, eq. (5.156) where the damping torque is equal `Dω`)
- `Ω`: rated frequency in [1/s] of the power grid, often ``2\pi⋅50``Hz
- `Γ`: voltage stability coefficient
- `V`: set voltage, usually `1`

# Mathematical Representation
Using `SwingEq` for node ``a`` applies the equations
```math
\frac{du_a}{dt} = j \omega_a u_a - \frac{u}{\|u\|} Γ_a  (v_a - V_a), \\
\frac{2H}{\Omega}\frac{d\omega_a}{dt} = P_a - D_a\omega_a - \Re\left(u_a \cdot i_a^*\right),
```
which is equivalent to
```math
\frac{d\phi_a}{dt} = \omega_a, \\
\frac{dv_a}{dt} = - Γ_a  (v_a - V_a) \\
\frac{2H}{\Omega}\frac{d\omega_a}{dt} = P_a - D_a\omega_a - \Re\left(u_a \cdot i_a^*\right),
```
"""
@DynamicNode SwingEqLVS(H, P, D, Ω, Γ, V) begin
    @assert D >= 0 "damping (D) should be >=0"
    @assert H > 0 "inertia (H) should be >0"
    @assert Γ > 0 "voltage magnitude stability coefficient (Γ) should be > 0"
    Ω_H = Ω / (2*H)
end [[ω, dω]] begin
    p = real(u * conj(i))
    dϕ = ω # dϕ is only a temp variable that Julia should optimize out
    dω = (P - D*ω - p)*Ω_H
    v = abs(u)
    # Linear Voltage Stability (LVS) term
    dv = - Γ * (v - V)
    du =  u/v * dv + u * im * dϕ
end

convert(::Type{SwingEq}, p::SwingEqLVS) = SwingEq(H=p.H, P=p.P, D=p.D, Ω=p.Ω)

export SwingEqLVS
