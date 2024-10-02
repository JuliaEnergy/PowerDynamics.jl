@mtkmodel Swing begin
    @structural_parameters begin
        ω_ref_input=false;
        Pm_input=false;
    end
    @components begin
        terminal = Terminal()
        if ω_ref_input
            ω_ref = RealInput()
        end
        if Pm_input
            Pm = RealInput()
        end
    end
    @variables begin
        ω(t), [guess=0, description="Rotor frequency"]
        θ(t), [guess=0, description="Rotor angle"]
        Pel(t), [description="Electrical Power injected into the grid"]
    end
    @parameters begin
        M=0.005, [description="Inertia"]
        D=0.0001, [description="Damping"]
        V, [guess=1, description="Voltage magnitude"]
        if !ω_ref_input
            ω_ref=1, [description="Reference frequency"]
        end
        if !Pm_input
            Pm, [guess=1, description="Mechanical Power"]
        end
    end
    begin
        _ω_ref = ω_ref_input ? ω_ref.u : ω_ref
        _Pm = Pm_input ? Pm.u : Pm
    end
    @equations begin
        Dt(θ) ~ ω - _ω_ref
        M*Dt(ω) ~ _Pm - D*(ω - _ω_ref) - Pel

        Pel ~ terminal.u_r*terminal.i_r + terminal.u_i*terminal.i_i
        terminal.u_r ~ V*cos(θ)
        terminal.u_i ~ V*sin(θ)
    end
end
