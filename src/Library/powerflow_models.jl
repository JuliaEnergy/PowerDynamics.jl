@cached_component function VδConstraint(; name, force_constraint=false)
    @named terminal = Terminal()
    params = @parameters begin
        V
        δ=0
    end
    eqs = if force_constraint
        [@no_simplify(terminal.u_r ~ V*cos(δ)),
         @no_simplify(terminal.u_i ~ V*sin(δ))]
    else
        [terminal.u_r ~ V*cos(δ),
         terminal.u_i ~ V*sin(δ)]
    end
    System(eqs, t, [], params; name, systems=[terminal])
end

@cached_component function UrUiConstraint(; name)
    @named terminal = Terminal()
    params = @parameters begin
        u_r, [guess=1, description="Real part of fixed terminal voltage"]
        u_i, [guess=0, description="Imaginary part of fixed terminal voltage"]
    end
    eqs = [
        terminal.u_r ~ u_r
        terminal.u_i ~ u_i
    ]
    System(eqs, t, [], params; name, systems=[terminal])
end

@cached_component function PVConstraint(; name)
    @named terminal = Terminal()
    params = @parameters begin
        P
        V
    end
    eqs = [
        P ~ terminal.u_r*terminal.i_r + terminal.u_i*terminal.i_i
        V^2 ~ terminal.u_r^2 + terminal.u_i^2
    ]
    System(eqs, t, [], params; name, systems=[terminal])
end

@cached_component function PQConstraint(; name)
    @named terminal = Terminal()
    params = @parameters begin
        P, [guess=0.1, description="Active Power [pu]"]
        Q, [guess=0, description="Reactive Power [pu]"]
    end
    S = P + im*Q
    uc = terminal.u_r + im*terminal.u_i
    ic = conj(S/uc)
    eqs = [
        terminal.i_r ~ simplify(real(ic))
        terminal.i_i ~ simplify(imag(ic))
    ]
    System(eqs, t, [], params; name, systems=[terminal])
end
