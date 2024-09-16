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
        ω(t)=0.0, [description="Rotor frequency"]
        θ(t)=0.0, [description="Rotor angle"]
        Pel(t), [description="Electrical Power injected into the grid"]
    end
    @parameters begin
        M=0.005, [description="Inertia"]
        D=0.1, [description="Damping"]
        V=1.0, [description="Voltage magnitude"]
        if !ω_ref_input
            ω_ref=0, [description="Reference frequency"]
        end
        if !Pm_input
            Pm, [description="Mechanical Power"]
        end
    end
    @equations begin
        if ω_ref_input
            Dt(θ) ~ ω - ω_ref.u
        else
            Dt(θ) ~ ω - ω_ref
        end
        if Pm_input
            Dt(ω) ~ 1/M * (Pm.u - D*ω - Pel)
        else
            Dt(ω) ~ 1/M * (Pm - D*ω - Pel)
        end

        Pel ~ terminal.u_r*terminal.i_r + terminal.u_i*terminal.i_i
        terminal.u_r ~ V*cos(θ)
        terminal.u_i ~ V*sin(θ)
    end
end
