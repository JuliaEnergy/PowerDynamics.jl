"""
    IdealDroopInverter

Simplified droop-controlled inverter with first-order power filtering and ideal voltage source output.

The droop law sets the inverter frequency `ω` in pu; the voltage angle `δ` is measured
against the **global dq frame**, hence `Dt(δ) ~ ωbase*(ω - ωframe)`. `ωset` (droop
setpoint) and `ωframe` (frame gauge) coincide numerically but are distinct quantities.

$(PowerDynamics.ref_source_file(@__FILE__, @__LINE__))
"""
@mtkmodel IdealDroopInverter begin
    @components begin
        terminal = Terminal()
    end

    @parameters begin
        Pset, [description="Active power setpoint [pu]", guess=1]
        Qset, [description="Reactive power setpoint [pu]", guess=0]
        Vset, [description="Voltage magnitude setpoint [pu]", guess=1]
        ωset=1, [description="Frequency setpoint [pu]"]
        ωbase, [bound_to = :systembase₊ωbase, description="System frequency base [rad/s]"]
        ωframe, [bound_to = :systembase₊ωframe, description="Global dq frame speed [pu]"]
        Kp=0.05, [description="Active power droop coefficient [pu freq / pu power]"]
        Kq=0.05, [description="Reactive power droop coefficient"]
        τ_p = 0.1, [description="Active Power filter time constant [s]"]
        τ_q = 0.1, [description="Reactive Power filter time constant [s]"]
    end

    @variables begin
        Pmeas(t), [description="Measured active power [pu]", guess=1]
        Qmeas(t), [description="Measured reactive power [pu]", guess=0]
        Pfilt(t), [description="Filtered active power [pu]", guess=1]
        Qfilt(t), [description="Filtered reactive power [pu]", guess=1]
        ω(t), [description="Frequency [pu]"]
        δ(t), [description="Voltage angle [rad]", guess=0]
        V(t), [description="Voltage magnitude [pu]"]
    end

    @equations begin
        ## Power measurement from terminal quantities
        Pmeas ~ terminal.u_r*terminal.i_r + terminal.u_i*terminal.i_i
        Qmeas ~ terminal.u_i*terminal.i_r - terminal.u_r*terminal.i_i

        ## First-order low-pass filtering
        τ_p * Dt(Pfilt) ~ Pmeas - Pfilt
        τ_q * Dt(Qfilt) ~ Qmeas - Qfilt

        ## Droop control equations (ωset is a setpoint, not the frame speed)
        ω ~ ωset - Kp * (Pfilt - Pset)  # Frequency decreases with excess power
        V ~ Vset - Kq * (Qfilt - Qset)  # Voltage decreases with excess reactive power

        ## Voltage angle dynamics (δ is measured against the global dq frame → ωframe)
        Dt(δ) ~ ωbase*(ω - ωframe)

        ## Output voltage components
        terminal.u_r ~ V*cos(δ)
        terminal.u_i ~ V*sin(δ)
    end
end;
