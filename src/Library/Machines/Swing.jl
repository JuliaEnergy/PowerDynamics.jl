"""
    Swing

Simplified swing-equation generator with constant voltage magnitude and active power balance.

The rotor angle `θ` is measured against the **global dq frame**, the rotor speed `ω` is in
pu on `ωbase`:

```math
\\begin{aligned}
\\dot{\\theta} &= \\omega_\\mathrm{base}\\,(\\omega - \\omega_\\mathrm{frame})\\\\
M\\,\\dot{\\omega} &= P_\\mathrm{m} - D\\,(\\omega - \\omega_\\mathrm{set}) - P_\\mathrm{el}
\\end{aligned}
```

Note the two distinct roles the frequency plays here: `ωframe` is the speed of the
reference frame (a gauge, pinned to 1), while `ωset` is the damping setpoint. They
coincide numerically today but are not the same quantity.

# Parameters
- `M`: Inertia ``M = 2H`` [s].
- `D`: Damping coefficient [pu power / pu frequency].
- `ωset`: Frequency setpoint [pu] the damping holds against.
- `Pm`: Mechanical power [pu] (or a `RealInput` if `Pm_input=true`).
- `V`: Terminal voltage magnitude [pu].
- `ωbase`/`ωframe` are inherited from the container's `systembase`; not constructor arguments.

$(PowerDynamics.ref_source_file(@__FILE__, @__LINE__))
"""
@mtkmodel Swing begin
    @structural_parameters begin
        Pm_input=false;
    end
    @components begin
        terminal = Terminal()
        if Pm_input
            Pm = RealInput()
        end
    end
    @variables begin
        ω(t), [guess=1, description="Rotor frequency [pu]"]
        θ(t), [guess=0, description="Rotor angle [rad]"]
        Pel(t), [description="Electrical Power injected into the grid [pu]"]
    end
    @parameters begin
        M=6.0, [description="Inertia M = 2H [s]"]
        D=2.0, [description="Damping [pu power / pu frequency]"]
        V, [guess=1, description="Voltage magnitude [pu]"]
        ωbase, [bound_to = :systembase₊ωbase, description="System frequency base [rad/s]"]
        ωframe, [bound_to = :systembase₊ωframe, description="Global dq frame speed [pu]"]
        ωset=1, [description="Frequency setpoint [pu]"]
        if !Pm_input
            Pm, [guess=1, description="Mechanical Power [pu]"]
        end
    end
    begin
        _Pm = Pm_input ? Pm.u : Pm
    end
    @equations begin
        # θ is measured against the global dq frame → ωframe (a gauge), while the
        # damping term references the commanded speed → ωset (a setpoint).
        Dt(θ) ~ ωbase*(ω - ωframe)
        M*Dt(ω) ~ _Pm - D*(ω - ωset) - Pel

        Pel ~ terminal.u_r*terminal.i_r + terminal.u_i*terminal.i_i
        terminal.u_r ~ V*cos(θ)
        terminal.u_i ~ V*sin(θ)
    end
end
