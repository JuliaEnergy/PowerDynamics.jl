# (C) 2018 Potsdam Institute for Climate Impact Research, authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)

@doc doc"""
According to P. Sauer "Power System Dynamics and Stability" p. 67 the
relationship between v_f_d and i_in1 is nonlinear due to saturation
of the exciter iron.

The saturation is then often described through two characteristic points on the
saturation curve: S_E_max (saturation at maximal excitation voltage)
and S_E0.75(saturation at 0.75* maximal excitation voltage).

These can be used to fit the following function:
S_E(E_f_d) = A_x * exp(B_x * E_f_d)               P.Sauer, P.70, Eq(4.22)

In a steady-state it is also given that:
0 = - (K_e + S_E_max) * E_f_d_max + V_R_max       P.Sauer, P.69, Eq(4.21)
# Inputs:
- `K_e`: Exciter constant related to self-excited field [pu]
- `S_E_max`: Saturation at maximal excitation voltage [pu]
- `S_E_tq`: Saturation at 0.75 * maximal excitation voltage [pu]
- `V_R_max`: Maximal reference voltage for the AVR [pu]

# Outputs:
- `E_f_d_max`: maximal excitation voltage
- `A_x`: Fitting Parameter
- `B_x`: Fitting Parameter
"""

function ExciterSaturationEq(K_e, S_E_max, S_E_tq, V_R_max)
    E_f_d_max = V_R_max / (K_e + S_E_max)
    A_x = S_E_tq^4 / S_E_max^3
    B_x = 4/3 * log(S_E_tq / A_x) / E_f_d_max
    return E_f_d_max, A_x, B_x
end

export ExciterSaturationEq
