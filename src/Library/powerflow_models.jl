@mtkmodel VδConstraint begin
    @structural_parameters begin
        force_constraint=false
    end
    @components begin
        terminal = Terminal()
    end
    @parameters begin
        V
        δ=0
    end
    @equations begin
        if force_constraint
            @no_simplify terminal.u_r ~ V*cos(δ)
            @no_simplify terminal.u_i ~ V*sin(δ)
        else
            terminal.u_r ~ V*cos(δ)
            terminal.u_i ~ V*sin(δ)
        end
    end
end

@mtkmodel UrUiConstraint begin
    @components begin
        terminal = Terminal()
    end
    @parameters begin
        u_r, [guess=1, description="Real part of fixed terminal voltage"]
        u_i, [guess=0, description="Imaginary part of fixed terminal voltage"]
    end
    @equations begin
        terminal.u_r ~ u_r
        terminal.u_i ~ u_i
    end
end

@mtkmodel PVConstraint begin
    @components begin
        terminal = Terminal()
    end
    @parameters begin
        P
        V
    end
    @equations begin
        P ~ terminal.u_r*terminal.i_r + terminal.u_i*terminal.i_i
        V^2 ~ terminal.u_r^2 + terminal.u_i^2
    end
end

@mtkmodel PQConstraint begin
    @components begin
        terminal = Terminal()
    end
    @parameters begin
        P, [guess=0.1, description="Active Power [pu]"]
        Q, [guess=0, description="Reactive Power [pu]"]
    end
    begin
        S = P + im*Q
        uc = terminal.u_r + im*terminal.u_i
        ic = conj(S/uc)
    end
    @equations begin
        terminal.i_r ~ simplify(real(ic))
        terminal.i_i ~ simplify(imag(ic))
    end
end
