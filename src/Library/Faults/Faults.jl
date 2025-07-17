@mtkmodel RXGroundFault begin
    @parameters begin
        R, [description="Resistance to ground during fault [pu]"]
        X, [description="Reactance to ground during fault [pu]"]
        enable=0, [description="Enable fault"]
    end
    @components begin
        terminal = Terminal()
    end
    @variables begin
        P(t), [description="Active Fault Power [pu]"]
        Q(t), [description="Reactive Fault Power [pu]"]
    end
    begin
        Z = R + im*X
        V = terminal.u_r + im*terminal.u_i
        iout = enable*V/Z
    end
    @equations begin
        terminal.i_r ~ simplify(real(iout))
        terminal.i_i ~ simplify(imag(iout))

        # observables
        P ~ terminal.u_r*terminal.i_r + terminal.u_i*terminal.i_i
        Q ~ terminal.u_i*terminal.i_r - terminal.u_r*terminal.i_i
    end
end
