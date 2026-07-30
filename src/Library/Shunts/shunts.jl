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
    DynamicCShunt(; B)

Dynamic shunt element modeled as a pure capacitor.

**Differential state: the capacitor voltage** `V_C_r`/`V_C_i` — suitable for DAE index
reduction at current-source buses and for modelling shunt capacitor banks without
resistive losses.

# Parameters
- `B`: Shunt susceptance [pu], evaluated at `ωbase` (physical capacitance is `B/ωbase` [pu·s]).
       Capacitive, i.e. `B > 0`; an inductive shunt has the reactor current as its state and
       is therefore a different model.

$(PowerDynamics.ref_source_file(@__FILE__, @__LINE__))
"""
@mtkmodel DynamicCShunt begin
    @components begin
        terminal = Terminal()
    end
    @parameters begin
        ωbase, [bound_to = :systembase₊ωbase, description="System frequency base [rad/s]"]
        ωframe, [bound_to = :systembase₊ωframe, description="Global dq frame speed [pu]"]
        B, [description="Shunt susceptance [pu] at ωbase"]
    end
    @variables begin
        V_C_r(t), [guess=1, description="Capacitor voltage real part (dq frame) [pu]"]
        V_C_i(t), [guess=0, description="Capacitor voltage imaginary part (dq frame) [pu]"]
        i_C_r(t), [guess=0, description="Capacitor current real part (dq frame) [pu]"]
        i_C_i(t), [guess=0, description="Capacitor current imaginary part (dq frame) [pu]"]
    end
    @equations begin
        # Capacitor dynamics in the rotating dq frame. The cross terms are the
        # transport-theorem jω contribution → frame speed ωframe [pu].
        B/ωbase * Dt(V_C_r) ~ -i_C_r + ωframe*B*V_C_i
        B/ωbase * Dt(V_C_i) ~ -i_C_i - ωframe*B*V_C_r
        # Terminal voltage = capacitor voltage
        terminal.u_r ~ V_C_r
        terminal.u_i ~ V_C_i
        # Grid current: capacitor only
        terminal.i_r ~ i_C_r
        terminal.i_i ~ i_C_i
    end
end

"""
    DynamicParallelRCShunt(; R, B)

Dynamic shunt element modeled as a parallel R ∥ C circuit.

This model represents a parallel combination of resistance and capacitance connected to a bus.
**Differential state: the capacitor voltage** `V_C_r`/`V_C_i`. Suitable for:
- Shunt capacitor banks with resistive losses
- Dynamic susceptance for DAE index reduction at current-source buses
- Fast transient behavior of reactive compensation

Note that the phasor admittance `Y = 1/R + jB` does not by itself determine the dynamics —
a series R–C shunt has the same `Y` at nominal frequency but different transients. The
topology is fixed by the model choice, after which `R` and `B` are unambiguous.

# Parameters
- `R`: Shunt resistance [pu]
- `B`: Shunt susceptance [pu], evaluated at `ωbase` (physical capacitance is `B/ωbase` [pu·s])

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
        B, [description="Shunt susceptance [pu] at ωbase"]
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
        # Capacitor dynamics in the rotating dq frame. The cross terms are the
        # transport-theorem jω contribution → frame speed ωframe [pu].
        B/ωbase * Dt(V_C_r) ~ -i_C_r + ωframe*B*V_C_i
        B/ωbase * Dt(V_C_i) ~ -i_C_i - ωframe*B*V_C_r
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
