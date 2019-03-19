# (C) 2018 Potsdam Institute for Climate Impact Research, authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)

@doc doc"""
```Julia
FourthOrderEqGovernorExciterAVR(H, P, D, Ω, T_d_dash ,T_q_dash ,X_q_dash ,X_d_dash,X_d, X_q, T_e, T_a, T_f, K_e, K_a, K_f, V_ref, R_d, T_sv, T_ch)
```

A node type that applies the 4th-order synchronous machine model with frequency/angle and voltage dynamics, including an Exciter, Automatic Voltage Regulator and Governor
which is implemented according to P. Sauer, "Power System Dynamics and Stability".
For an illustration of a synchronous machine schematic see P. Sauer, Fig. 3.1 on p. 25.

Exciter and Automatic Voltage Regulator:
The equations for the systems that balance the AC synchronous machine voltage level by increasing or decreasing the exciter DC voltage. Note, within this model, the transient reactance in the d-axis
of the generator needs to be included into the nodal admittance matrix. As the bus of this generator node type is constructed to be an internal generator bus.

Governor:
The prime mover provides the mechanism for controlling the synchronous machine speed and, hence, terminal voltage frequency.

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
- `Ω`: rated frequency of the power grid, often 50Hz
- `T_d_dash`: time constant of d-axis, given in [s], see P. Sauer, chapter 3.7, p. 54 for a general explanation on time constants
- `T_q_dash`: time constant of q-axis, given in [s]
- `X_d_dash`: transient reactance of d-axis, given in [pu]
- `X_q_dash`: transient reactance of q-axis, given in [pu]
- `X_d`: reactance of d-, given in [pu]
- `X_d`: reactance of q-axis, given in [pu]

- 'T_e' : Exciter time constant, integration rate associated with exciter control [s]
- 'T_a' : Maximum voltage regulator output [s]
- 'T_f' : Excitation control system stabilizer time constant [s]
- 'K_e' : Exciter constant related to self-excited field [pu]
- 'K_a' : Voltage Regulator gain [pu]
- 'K_f' : Excitation control system stabilizer gains [pu]
- 'V_ref' : Reference voltage for the AVR [pu]

- 'R_d' : Speed regulation [2πdroop/ω_s]
- 'T_sv' : Steam Valve time constant [s]
- 'T_ch' : Steam Chest time constant [s]


# Mathematical Representation Synchronous Machine
Using `FourthEq` for node ``a`` applies the equations
```math
    u = -je_c e^{j\theta} = -j(e_d + je_q)e^{j\theta}\\
    e_c= e_d + je_q = jue^{-j\theta}\\
    i  = -ji'e^{j\theta} = -j(i_d+ j i_q )e^{j\theta} = Y^L \cdot u \\
    i_c= i_d + ji_q = jie^{-j\theta}\\
    p = \Re (i^* u) \\
```

where complex voltage and current are described in a co-rotating frame with axes labeled d and q.

The fourth-order equations read (according to P. Sauer, "Power System Dynamics and Stability", p. 140, eqs. (6110)-(6114)) and p. 35 eqs(3.90)-(3.91)
```math
    \frac{d\theta}{dt} = \omega \\
    \frac{d\omega}{dt} = (P-D\omega - p -(X'_q-X'_d)i_d i_q)Ω_H\\
    \frac{d e_q}{dt} = \frac{1}{T'_d} (- e_q - (X_d - X'_d) i_{d}+ E_f) \\
    \frac{d e_d}{dt} = \frac{1}{T'_q} (- e_d + (X_q - X'_q) i_{q}) \\
```

# Exciter and AVR equations
```math
	u_{terminal} = e'_c - j X'_d i \\
	S_{e}(e_{fd}) = 0.098e^{0.55 e_{fd}} (according to P. Sauer, p. 70) \\
	\dfrac{dR_f}{dt} = \dfrac{1}{T_f} (-R_f + \dfrac{K_f}{T_f} e_f) \\
	\dfrac{dv_r}{dt} = \dfrac{1}{T_a} (-v_r + (K_a R_f) -\dfrac{K_a K_f}{T_f}e_{fd} + K_a (V_{ref} - abs(u_{terminal}))) \\
	\dfrac{de_{fd}}{dt} = \dfrac{1}{T_e} (-K_e + S_{e}(e_{fd})e_{fd} + v_r) \\
```

# Governor equations
```math
    \dfrac{dP_m}{dt} = \dfrac{1}{T_{ch}} (-P_m + P_{sv}) \\
    Assumption: T_m = P_m \\
    \dfrac{dP_{sv}}{dt} = \dfrac{1}{T_{sv}} (-P_{sv} + P_c -\dfrac{1}{R_d} (\dfrac{\omega}{\omega_s} - 1)) \\
```

The equations for frequency and phase represent energy conservation and phase shift.
The dynamic equations for the complex voltage show the relationship between the dynamicy of flux linkages and currents which must reflect a
conservative coupling field.

With the PowerDynamics.jl naming conventions of ``i`` and ``u`` they read as
```math
   \dot u = \frac{d}{dt}(-j e_c e^{j\theta})=-j(\dot e_d + j\dot e_q)e^{j\theta} + uj\omega \\
```
"""
@DynamicNode FourthOrderEqGovernorExciterAVR(H, P, D, Ω, T_d_dash ,T_q_dash ,X_q_dash ,X_d_dash,X_d, X_q, T_e, T_a, T_f, K_e, K_a, K_f, V_ref, R_d, T_sv, T_ch) <: OrdinaryNodeDynamics() begin # S_e_fd,
    @assert H > 0 "inertia (H) should be >0"
    @assert D >= 0 "damping (D) should be >=0"
    @assert T_d_dash > 0 "time constant of d-axis (T_d_dash) should be >0"
    @assert T_q_dash > 0 "time constant of q-axis (T_q_dash) should be >0"
    @assert X_d_dash >= 0 "transient reactance of d-axis (X_d_dash) should be >=0"
    @assert X_q_dash >= 0 "transient reactance of q-axis (X_q_dash) should be >=0"
    @assert X_d >= 0 "reactance of d-axis (X_d_dash) should be >=0"
    @assert X_q >= 0 "reactance of q-axis (X_q_dash) should be >=0"
    @assert T_e > 0 "Exciter time constant"
    @assert T_a > 0 "Maximum voltage regulator output"
    @assert T_f > 0 "Excitation control system stabilizer time constant"
    @assert K_a > 0 "Voltage Regulator gain should be >0"
    @assert K_f > 0 "Excitation control system stabilizer gains should be >0"
    @assert V_ref > 0 "Reference voltage for the AVR should be >0"
    @assert R_d > 0 "Speed regulation should be >0"
    @assert T_sv > 0 "Steam Valve time constant should be >0"
    @assert T_ch > 0 "T_ch' : Steam Chest time constant should be >0"

    Ω_H = (Ω * 2pi) / H #    norm = 2 * H / (2 * np.pi * 50)  # normalize the parameters as done for coupling_const, input_power, damping_const

end [[θ, dθ],[ω, dω],[e_f, de_f],[v_r, dv_r],[r_f,dr_f],[P_sv, dP_sv],[P_m, dP_m]] begin
    ##### Model R1 #####
    i_c = 1im*i*exp(-1im*θ)
    e_c = 1im*u*exp(-1im*θ)
    p = real(u * conj(i))
    e_d = real(e_c)
    e_q = imag(e_c)
    i_d = real(i_c)
    i_q = imag(i_c)

    #Exciter
    V_mes = e_c - 1im*X_d_dash*i_c
    dr_f = (1 / T_f) * (- r_f + ((K_f/T_f) * e_f))
    dv_r = (1 / T_a) * (- v_r + (K_a * r_f) - ((K_a * K_f)/T_f)*e_f + K_a*(V_ref - abs(V_mes)))
    de_f = (1 / T_e) * ((- (K_e + (0.098*exp(0.55*e_f))) * e_f) + v_r) #S_e_fd, Sauer p70

    #Governor
    dP_sv = (1 / T_sv) * (-P_sv + P - (1/R_d)*(((ω+(Ω*2pi))/(Ω*2pi))-1))
    dP_m  = (1 / T_ch) * (-P_m  + P_sv)

    #Synchronous machine
    dθ = ω
    de_d = (1 / T_q_dash) * (- e_d + (X_q - X_q_dash) * i_q)
    de_q = (1 / T_d_dash) * (- e_q - (X_d - X_d_dash) * i_d + e_f)
    de_c = de_d + 1im*de_q
    du = -1im*de_c*exp(1im*θ)+ u*1im*ω
    dω = (P_m - D*ω - p - (X_q_dash - X_d_dash) * i_d * i_q) * Ω_H
end

export FourthOrderEqGovernorExciterAVR
