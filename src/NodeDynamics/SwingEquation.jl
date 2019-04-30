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

# Mathematical Representation
Using `SwingEq` for node ``a`` applies the equations
```math
\frac{du_a}{dt} = i u_a  \omega_a, \\
\frac{H}{2\pi\Omega}\frac{d\omega_a}{dt} = P_a - D_a\omega_a - \Re\left(u_a \cdot i_a^*\right),
```
which is equivalent to
```math
\frac{d\phi_a}{dt} = \omega, \\
v = v(t=0) = \text{const.} \\
\frac{H}{2\pi\Omega}\frac{d\omega_a}{dt} = P_a - D_a\omega_a - \Re\left(u_a \cdot i_a^*\right),
```
"""
@DynamicNode SwingEq(H, P, D, Ω) begin
    @assert D >= 0 "damping (D) should be >=0"
    @assert H > 0 "inertia (H) should be >0"
    Ω_H = Ω * 2pi / H
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
- `H`: inertia
- `P`: active (real) power output
- `D`: damping coefficient
- `Ω`: rated frequency of the power grid, often 50Hz
- `Γ`: voltage stability coefficient
- `V`: set voltage, usually `1`

# Mathematical Representation
Using `SwingEq` for node ``a`` applies the equations
```math
\frac{du_a}{dt} = i u_a \omega - \frac{u}{\|u\|} Γ_a  (v_a - V_a), \\
\frac{H}{2\pi\Omega}\frac{d\omega_a}{dt} = P_a - D_a\omega_a - \Re\left(u_a \cdot i_a^*\right),
```
which is equivalent to
```math
\frac{d\phi_a}{dt} = \omega_a, \\
\frac{dv_a}{dt} = - Γ_a  (v_a - V_a) \\
\frac{H}{2\pi\Omega}\frac{d\omega_a}{dt} = P_a - D_a\omega_a - \Re\left(u_a \cdot i_a^*\right),
```
"""
@DynamicNode SwingEqLVS(H, P, D, Ω, Γ, V) begin
    @assert Γ > 0 "voltage magnitude stability coefficient (Γ) should be > 0"
    # FIXME
    # the following is a bit of a hack and should only be done if you really know what you do
    swing_node_dyn = construct_node_dynamics(convert(SwingEq, par))
    #@assert typeof(swing_node_dyn) === OrdinaryNodeDynamics{SwingEq}
    #@assert swing_node_dyn.n_int == 1
    swing_eq_rhs! = swing_node_dyn.f!
end [[ω, dω]] begin
    v = abs(u)
    # Linear Voltage Stability (LVS) term
    dv = - Γ * (v - V)

    du =  u/v * dv + swing_eq_rhs!(dint, u, i, int, t)
    dω = dint[1]
end
convert(::Type{SwingEq}, p::SwingEqLVS) =
    SwingEq(H=p.H, P=p.P, D=p.D, Ω=p.Ω)

export SwingEqLVS
