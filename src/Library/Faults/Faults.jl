"""
    RXGroundFault

Shunt RX impedance to ground representing a fault, switchable via the `enable` parameter.

$(PowerDynamics.ref_source_file(@__FILE__, @__LINE__))
"""
@cached_component function RXGroundFault(; name)
    @named terminal = Terminal()
    params = @parameters begin
        R, [description="Resistance to ground during fault [pu]"]
        X, [description="Reactance to ground during fault [pu]"]
        enable=0, [description="Enable fault"]
    end
    vars = @variables begin
        P(t), [description="Active Fault Power [pu]"]
        Q(t), [description="Reactive Fault Power [pu]"]
    end
    Z = R + im*X
    V = terminal.u_r + im*terminal.u_i
    iout = enable*V/Z
    eqs = [
        terminal.i_r ~ simplify(real(iout))
        terminal.i_i ~ simplify(imag(iout))

        # observables
        P ~ terminal.u_r*terminal.i_r + terminal.u_i*terminal.i_i
        Q ~ terminal.u_i*terminal.i_r - terminal.u_r*terminal.i_i
    ]
    System(eqs, t, vars, params; name, systems=[terminal])
end
