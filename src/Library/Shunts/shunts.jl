"""
    StaticShunt(; G, B)

Static (algebraic) shunt element with constant admittance Y = G + jB.

This model represents a linear shunt connected to a bus, drawing current proportional
to the bus voltage: I = Y * V. It has no dynamic states and is evaluated algebraically.

# Parameters
- `G`: Shunt conductance [pu]. Positive values represent resistive losses (real power consumption).
- `B`: Shunt susceptance [pu]. Positive values represent capacitive behavior (leading current),
       negative values represent inductive behavior (lagging current).
"""
@mtkmodel StaticShunt begin
    @parameters begin
        G, [description="Shunt conductance [pu]"]
        B, [description="Shunt susceptance [pu]"]
    end
    @components begin
        terminal = Terminal()
    end
    begin
        Y = G + im*B
        ishunt = -Y * (terminal.u_r + im*terminal.u_i)
    end
    @equations begin
        terminal.i_r ~ real(ishunt)
        terminal.i_i ~ imag(ishunt)
    end
end

"""
    DynamicRCShunt(; R, C, ω0=2π*50)

Dynamic shunt element modeled as a parallel R || C circuit.

This model represents a parallel combination of resistance and capacitance connected to a bus.
The capacitor voltage is a differential state, suitable for:
- Shunt capacitor banks with resistive losses
- Dynamic susceptance for DAE index reduction at current-source buses
- Fast transient behavior of reactive compensation

# Parameters
- `R`: Shunt resistance [pu]
- `C`: Shunt susceptance [pu] at frequency ω0. Related to physical capacitance by C = ω0 * C_actual.
- `ω0`: Frame angular frequency [rad/s]. Default: 2π*50 rad/s.
"""
@mtkmodel DynamicRCShunt begin
    @components begin
        terminal = Terminal()
    end
    @parameters begin
        ω0=2π*50, [description="Frame angular frequency [rad/s]"]
        R, [description="Shunt resistance [pu]"]
        C, [description="Shunt susceptance [pu] (frequency-normalized capacitance)"]
    end
    @variables begin
        V_C_r(t), [guess=1, description="Capacitor voltage real part (dq frame) [pu]"]
        V_C_i(t), [guess=0, description="Capacitor voltage imaginary part (dq frame) [pu]"]
        i_C_r(t), [guess=0, description="Capacitor current real part (dq frame) [pu]"]
        i_C_i(t), [guess=0, description="Capacitor current imaginary part (dq frame) [pu]"]
        i_R_r(t), [guess=0, description="Resistor current real part (dq frame) [pu]"]
        i_R_i(t), [guess=0, description="Resistor current imaginary part (dq frame) [pu]"]
    end
    @equations begin
        # Capacitor dynamics in rotating dq frame (C is susceptance)
        C/ω0 * Dt(V_C_r) ~ -i_C_r + C*V_C_i
        C/ω0 * Dt(V_C_i) ~ -i_C_i - C*V_C_r
        # Resistor current (Ohm's law)
        i_R_r ~ -terminal.u_r / R
        i_R_i ~ -terminal.u_i / R
        # Terminal voltage = capacitor voltage
        terminal.u_r ~ V_C_r
        terminal.u_i ~ V_C_i
        # Grid current: capacitor + resistor (parallel connection)
        terminal.i_r ~ i_C_r + i_R_r
        terminal.i_i ~ i_C_i + i_R_i
    end
end
