# (C) 2018 Potsdam Institute for Climate Impact Research, authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)

@doc doc"""
```Julia
ExponentialRecoveryLoad(P0, Q0, Nps, Npt, Nqs, Nqt, Tp, Tq, V0)
```
A node type that represents the exponential recovery load model. The exponential recovery load model aims to
capture the load restoration characteristics with an exponential recovery process expressed as an
input–output relationship between powers (real and reactive) and voltage.

# Keyword Arguments
- `P0`: Active power load demand [pu]
- `Q0`: Reactive power load demand [pu]
- `Nps`: Steady-state load voltage dependence p-axis [pu]
- `Npt`: Transient load voltage dependence p-axis [pu]
- `Nqs`: Steady-state load voltage dependence q-axis [pu]
- `Nqt`: Transient load voltage dependence q-axis [ pu]
- `Tp`: Load recovery constant p-axis [s]
- `Tq`: Load recovery constant q-axis [s]
- `V0`: Reference grid voltage [pu]

# Mathematical Representation
```math
	\dfrac{dx_p}{dt} = \dfrac{1}{T_p}(-x_p + P_0(\dfrac{|u|}{V_0})^{N_{ps}} - P_0(\dfrac{|u|}{V_0})^{N_{pt}}) \\
    \dfrac{dx_q}{dt} = \dfrac{1}{T_q}(-x_q + Q_0(\dfrac{|u|}{V_0})^{N_{qs}} - Q_0(\dfrac{|u|}{V_0})^{N_{qt}}) \\
    P = x_p + P_0(\dfrac{|u|}{V_0})^{N_{pt}} \\
    Q = x_q - Q_0(\dfrac{|u|}{V_0})^{N_{qt}} \\
```

IEEE TRANSACTIONS ON POWER SYSTEMS, VOL. 21, NO. 3, AUGUST 2006
Measurement-Based Dynamic Load Models: Derivation, Comparison, and Validation
Byoung-Kon Choi, Member, IEEE, Hsiao-Dong Chiang, Fellow, IEEE, Yinhong Li, Hua Li, Member, IEEE, Yung-Tien Chen, Der-Hua Huang, and Mark G. Lauby
"""
@DynamicNode ExponentialRecoveryLoad(P0, Q0, Nps, Npt, Nqs, Nqt, Tp, Tq, V0) begin
    MassMatrix(m_int = [true, true])
end  begin
    @assert V0 > 0 "Nominal Voltage should be >0"
    @assert Tp > 0 "Load recovery constant should be >0"
    @assert Tq > 0 "Load recovery constant should be >0"
    @assert Nps > 0 "Steady-state load voltage dependence p-axis should be >0"
    @assert Npt > 0 "Transient load voltage dependence p-axis should be >0"
    @assert Nqs > 0 "Steady-state load voltage dependence q-axis should be >0"
    @assert Nqt > 0 "Transient load voltage dependence p-axis should be >0"

end [[x_p, dx_p],[x_q, dx_q]] begin
    Pd = real(u*conj(i))
    Qd = imag(u*conj(i))

    dx_p = (1/Tp)*(-x_p + P0*((abs(u)/V0)^Nps) - P0*((abs(u)/V0)^Npt))
    dx_q = (1/Tq)*(-x_q + Q0*((abs(u)/V0)^Nqs) - Q0*((abs(u)/V0)^Nqt))

    du = -Pd + x_p + P0*((abs(u)/V0)^Npt) + im*(-Qd + x_q + Q0*((abs(u)/V0)^Nqt))
end

export ExponentialRecoveryLoad
