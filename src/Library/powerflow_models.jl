@mtkmodel VδConstraint begin
    @components begin
        terminal = Terminal()
    end
    @parameters begin
        V
        δ=0
    end
    @equations begin
        terminal.u_r ~ V*cos(δ)
        terminal.u_i ~ V*sin(δ)
    end
end

@mtkmodel UrUiConstraint begin
    @components begin
        terminal = Terminal()
    end
    @parameters begin
        u_r
        u_i
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
        P
        Q
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
