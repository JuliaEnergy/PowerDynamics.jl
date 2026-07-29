"""
    StaticShunt(; G, B)

Static (algebraic) shunt element with constant admittance Y = G + jB.

This model represents a linear shunt connected to a bus, drawing current proportional
to the bus voltage: I = Y * V. It has no dynamic states and is evaluated algebraically.

# Parameters
- `G`: Shunt conductance [pu]. Positive values represent resistive losses (real power consumption).
- `B`: Shunt susceptance [pu]. Positive values represent capacitive behavior (leading current),
       negative values represent inductive behavior (lagging current).

$(PowerDynamics.ref_source_file(@__FILE__, @__LINE__))
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
    DynamicCShunt(; C)

Dynamic shunt element modeled as a pure capacitor.

The capacitor voltage is a differential state, suitable for DAE index reduction at
current-source buses and modelling shunt capacitor banks without resistive losses.

# Parameters
- `C`: Shunt susceptance [pu] at frequency `ωbase`. Related to physical capacitance by C = ωbase * C_actual.

$(PowerDynamics.ref_source_file(@__FILE__, @__LINE__))
"""
@mtkmodel DynamicCShunt begin
    @components begin
        terminal = Terminal()
    end
    @parameters begin
        ωbase, [bound_to = :systembase₊ωbase, description="System frequency base [rad/s]"]
        ωframe, [bound_to = :systembase₊ωframe, description="Global dq frame speed [pu]"]
        C, [description="Shunt susceptance [pu] (frequency-normalized capacitance)"]
    end
    @variables begin
        V_C_r(t), [guess=1, description="Capacitor voltage real part (dq frame) [pu]"]
        V_C_i(t), [guess=0, description="Capacitor voltage imaginary part (dq frame) [pu]"]
        i_C_r(t), [guess=0, description="Capacitor current real part (dq frame) [pu]"]
        i_C_i(t), [guess=0, description="Capacitor current imaginary part (dq frame) [pu]"]
    end
    @equations begin
        # Capacitor dynamics in the rotating dq frame (C is susceptance). The cross
        # terms are the transport-theorem jω contribution → frame speed ωframe [pu].
        C/ωbase * Dt(V_C_r) ~ -i_C_r + ωframe*C*V_C_i
        C/ωbase * Dt(V_C_i) ~ -i_C_i - ωframe*C*V_C_r
        # Terminal voltage = capacitor voltage
        terminal.u_r ~ V_C_r
        terminal.u_i ~ V_C_i
        # Grid current: capacitor only
        terminal.i_r ~ i_C_r
        terminal.i_i ~ i_C_i
    end
end

"""
    DynamicParallelRCShunt(; R, C)

Dynamic shunt element modeled as a parallel R ∥ C circuit.

This model represents a parallel combination of resistance and capacitance connected to a bus.
The capacitor voltage is a differential state, suitable for:
- Shunt capacitor banks with resistive losses
- Dynamic susceptance for DAE index reduction at current-source buses
- Fast transient behavior of reactive compensation

# Parameters
- `R`: Shunt resistance [pu]
- `C`: Shunt susceptance [pu] at frequency `ωbase`. Related to physical capacitance by C = ωbase * C_actual.

$(PowerDynamics.ref_source_file(@__FILE__, @__LINE__))
"""
@mtkmodel DynamicParallelRCShunt begin
    @components begin
        terminal = Terminal()
    end
    @parameters begin
        ωbase, [bound_to = :systembase₊ωbase, description="System frequency base [rad/s]"]
        ωframe, [bound_to = :systembase₊ωframe, description="Global dq frame speed [pu]"]
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
        # Capacitor dynamics in the rotating dq frame (C is susceptance). The cross
        # terms are the transport-theorem jω contribution → frame speed ωframe [pu].
        C/ωbase * Dt(V_C_r) ~ -i_C_r + ωframe*C*V_C_i
        C/ωbase * Dt(V_C_i) ~ -i_C_i - ωframe*C*V_C_r
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
