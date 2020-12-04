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

In addition to the 4th-order model the Type IEEEDC1A exciter was implemented.
The type was implemented using the following resources as the guide line: 
"MatDyn" Copyright ©2009 Stijn Cole and 
"IEEE Recommended Practice for Excitation System Models for Power System Stability Studies", IEEE Power and Energy Society, 2016

Additionally to ``u``, it has the internal dynamic variables
* ``ω`` representing the frequency of the rotator relative to the grid frequency ``Ω``, i.e. the real frequency ``ω_r`` of the rotator is given as ``\omega_r = \Omega + \omega``
* ``θ`` representing the relative angle of the rotor with respect to the voltage angle ``ϕ``.
* ``E_f`` representing the scaled field voltage, The exciter model (here IEEEDC1A) proviedes the value for ``E_f``(according to P. Sauer, p. 65)
* ``U_r`` representing the Voltage regulator output in [pu]
* ``U_f`` representing the excitation system stabilization rate feedback

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


# IEEE DC1A Exciter Model
- `K_e`: Exciter constant related to self-excited field [pu]
- `K_f`: Excitation control system stabilizer gains [pu]
- `K_a`: Voltage Regulator gain [pu]
- `U`:
- `U_ref`: Reference value of the stator terminal voltage [pu]
- `U_ref2`: Reference value of the stator terminal voltage [pu]
- `U_rmax`: Voltage regulator maximum output [pu]
- `U_rmin`: Voltage regulator minimum output [pu]
- `T_a`: Time constant of the voltage regulator [s]
- `T_f`: Excitation control system stabilizer time constant [s]
- `T_e`: Exciter time constant, integration rate associated with exciter control [s]
"""

@DynamicNode FourthOrderEqExciterIEEED1A(H, P, D, Ω, T_d_dash, T_q_dash, X_d_dash, X_q_dash, X_d, X_q, K_e, K_f, K_a, U, U_ref, U_ref2, U_rmax, U_rmin, T_a, T_f, T_e, S_E_max, S_E_tq, V_R_max) begin
    @assert H > 0 "inertia (H) should be >0"
    @assert D >= 0 "damping (D) should be >=0"
    @assert T_d_dash > 0 "time constant of d-axis (T_d_dash) should be >0"
    @assert X_d_dash >= 0 "transient reactance of d-axis (X_d_dash) should be >=0"
    @assert X_d >= 0 "reactance of d-axis (X_d_dash) should be >=0"

    Ω_H = Ω / (2 * H)
end [[θ,dθ],[ω, dω], [E_f, dE_f], [U_f, dU_f], [U_r, dU_r]] begin
    i_c = 1im*i*exp(-1im*θ)
    e_c = 1im*u*exp(-1im*θ)
    p = real(u * conj(i))
    e_d = real(e_c)
    e_q = imag(e_c)
    i_d = real(i_c)
    i_q = imag(i_c)

    # Exciter Model, IEEEDC1A
	E_f_d_max, A_x, B_x = ExciterSaturationEq(K_e, S_E_max, S_E_tq, V_R_max) # Saturation modeling
	U_x = A_x * exp(B_x * E_f)

	dU_r = 1 / T_a * (K_a * (U_ref - U + U_ref2 - U_f) - U_r)
	dU_f = 1 / T_f * (K_f / T_e * (U_r - U_x - K_e * E_f) - U_f)

    if U_r > U_rmax
		println("Ur > Urmax")
		U_r2 = U_rmax

	elseif U_r < U_rmin
		println("Ur < Urmax")
		U_r2 = U_rmin
	else
		U_r2 = U_r
	end

	dE_f = 1 / T_e * ( U_r2 - U_x - K_e * E_f)

    # Fourth Order Model
	dθ = ω
    de_d = (1 / T_q_dash) * (- e_d + (X_q - X_q_dash) * i_q)
    de_q = (1 / T_d_dash) * (- e_q - (X_d - X_d_dash) * i_d + E_f)
    de_c = de_d + 1im * de_q
    du = -1im * de_c * exp(1im*θ) + u * 1im * ω
    dω = (P - D * ω - p - (X_q_dash - X_d_dash) * i_d * i_q) * Ω_H

end

export FourthOrderEqExciterIEEED1A
