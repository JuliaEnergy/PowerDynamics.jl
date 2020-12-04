# (C) 2018 Potsdam Institute for Climate Impact Research, authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)

@doc doc"""
A node type that applies the 4th-order synchronous machine model with frequency/angle and voltage dynamics,
which is implemented according to P. Sauer, "Power System Dynamics and Stability".
For an illustration of a synchronous machine schematic see P. Sauer, Fig. 3.1 on p. 25.
Usually the swing equation (`SwingEq`) is used for short time periods to analyze the transient behavior
of generators in a power grid, the so-called first swing. The 4th-order model  also takes the back reaction
of the power flow onto the voltage into account. This has the effect that the angle of the voltage as seen by the power grid,
and the angle of the rotating mass are no longer the same but become dynamically coupled.

In addition to the 4th-order model the  Type IEEEG1 governor was implemented.
It is used for tandem compound, double reheat steam turbine systems.
The type was implemented using the following resources as the guide line: 
"MatDyn" Copyright ©2009 Stijn Cole and
"Dynamic Models for Turbine-Governors in Power System Studies", IEEE PES, 2013, Page 2 - 2

Additionally to ``u``, it has the internal dynamic variables
* ``ω`` representing the frequency of the rotator relative to the grid frequency ``Ω``, i.e. the real frequency ``ω_r`` of the rotator is given as ``\omega_r = \Omega + \omega`` and
* ``θ`` representing the relative angle of the rotor with respect to the voltage angle ``ϕ``.
* ``Pm`` representing the active power output, also called the mechanical torque applied to the shaft, given in [pu]
* ``x1`` representing an internal variable block input
* ``z`` representing an internal variable block output
* ``P`` representing the power of the servomotor [pu]

# Keyword Arguments
- `H`: shaft inertia constant (given in [s]), defined according to P. Sauer, p. 33, eq. (3.60)
- `D`: damping coefficient (given in [s], see P. Sauer, eq. (5.156) where the damping torque is equal `Dω`)
- `Ω`: rated frequency of the power grid, often ``2π⋅50Hz``
- `T_d_dash`: time constant of d-axis, given in [s], see P. Sauer, chapter 3.7, p. 54 for a general explanation on time constants
- `T_q_dash`: time constant of q-axis, given in [s]
- `X_d_dash`: transient reactance of d-axis, given in [pu]
- `X_q_dash`: transient reactance of q-axis, given in [pu]
- `X_d`: reactance of d-, given in [pu]
- `X_d`: reactance of q-axis, given in [pu]
- `E_f`: scaled field voltage, which, if set equal to 1.0 pu, gives 1.0 pu open-circuit terminal voltage. The physical device that provides the value of `E_f` is called the exciter (according to P. Sauer, p. 65)

# Governer
- `P0`: Reference power
- `Pmax`: Max power limit imposed by valve or gate control [pu]
- `Pmin`: Min power limit imposed by valve or gate control [pu]
- `Pup`: Max main control valve rate of change [pu/s]
- `Pdown`: Min main control valve rate of change [pu/s]
- `T1`: Controller lag compensation [s]
- `T2`: Controller lead compensation [s]
- `T3`: Valve position time constant (servomotor mechanism) [s]
- `K`: Total effective speed-governing system gain [pu]
"""

@DynamicNode FourthOrderEqGovernorIEEEG1(H, D, Ω, E_f, T_d_dash, T_q_dash, X_d_dash, X_q_dash, X_d, X_q, P0, Pmax, Pmin, Pup, Pdown, T1, T2, T3, K) begin
    @assert H > 0 "inertia (H) should be >0"
    @assert D >= 0 "damping (D) should be >=0"
    @assert T_d_dash > 0 "time constant of d-axis (T_d_dash) should be >0"
    @assert X_d_dash >= 0 "transient reactance of d-ax11is (x11_d_dash) should be >=0"
    @assert X_d >= 0 "reactance of d-ax11is (x11_d_dash) should be >=0"
    Ω_H = Ω / (2 * H)
end [[θ,dθ],[ω, dω], [Pm, dPm], [x1, dx1], [z, dz], [P, dP]] begin
    i_c = 1im*i*exp(-1im*θ)
    e_c = 1im*u*exp(-1im*θ)
    p = real(u * conj(i))
    e_d = real(e_c)
    e_q = imag(e_c)
    i_d = real(i_c)
    i_q = imag(i_c)

    # Governor Type IEEEG1
    dx1 = K * (-1 / T1 * x1 + (1 - T2 /  T1) * ω) # Block Input

    dP = (1 / T1) * x1 + (T2 / T1) * ω

    y = (1 / T3) * (P0 - P - Pm)                  # Block Output
    y_temp = y                                    # temorary variable

    # Limiting the valve rate of change
    if y > Pup
        y_temp = Pup
    end

    if y < Pdown
        y_temp = Pdown
    end

    dz = y_temp
    dPm = y_temp

    # Limiting the power imposed by the valve or gate control
    if z > Pmax
        dPm = (1 - Pmax) * dPm
    end

    if z < Pmin
        dPm = (1 - Pmin) * dPm
    end

    # Fourth Order Model
    dθ = ω
    de_d = (1 / T_q_dash) * (- e_d + (X_q - X_q_dash) * i_q)
    de_q = (1 / T_d_dash) * (- e_q - (X_d - X_d_dash) * i_d + E_f)
    de_c = de_d + 1im * de_q
    du = -1im * de_c * exp(1im * θ) + u * 1im * ω
    dω = (Pm - D * ω - p - (X_q_dash - X_d_dash) * i_d * i_q) * Ω_H
end

export FourthOrderEqGovernorIEEEG1
