@mtkmodel PQLoad begin
    @parameters begin
        Pset, [description="Active Power demand"]
        Qset, [description="Reactive Power demand"]
    end
    @components begin
        terminal = Terminal()
    end
    @variables begin
        P(t), [description="Active Power"]
        Q(t), [description="Reactive Power"]
    end
    @equations begin
        P ~ terminal.u_r*terminal.i_r + terminal.u_i*terminal.i_i
        Q ~ terminal.u_i*terminal.i_r - terminal.u_r*terminal.i_i
        # makes it easier for the solver then P ~ Pset and Q ~ Qset
        terminal.i_r ~  simplify(real((Pset + im*Qset)/(terminal.u_r + im*terminal.u_i)))
        terminal.i_i ~ -simplify(imag((Pset + im*Qset)/(terminal.u_r + im*terminal.u_i)))
    end
end

@mtkmodel VoltageDependentLoad begin
    @parameters begin
        Pset, [guess=-1, description="Active Power demand"]
        Qset, [guess=0, description="Reactive Power demand"]
        αP, [description="Active Power exponent"]
        αQ, [description="Reactive Power exponent"]
        Vn, [guess=1, description="Nominal voltage (where real power equals set power)"]
    end
    @components begin
        terminal = Terminal()
    end
    @variables begin
        P(t), [description="Active Power [pu]"]
        Q(t), [description="Reactive Power [pu]"]
    end
    begin
        v = sqrt(terminal.u_r^2 + terminal.u_i^2)
        Pload = Pset * (v/Vn)^αP
        Qload = Qset * (v/Vn)^αQ
        Sload = Pload + im*Qload
        vcomplex = terminal.u_r + im*terminal.u_i
        iout = conj(Sload/vcomplex)
    end
    @equations begin
        terminal.i_r ~ simplify(real(iout))
        terminal.i_i ~ simplify(imag(iout))

        # observables
        P ~ terminal.u_r*terminal.i_r + terminal.u_i*terminal.i_i
        Q ~ terminal.u_i*terminal.i_r - terminal.u_r*terminal.i_i
    end
end

@mtkmodel ConstantYLoad begin
    @structural_parameters begin
        allow_zero_conductance = false
    end
    @components begin
        terminal = Terminal()
    end
    @parameters begin
        B, [guess=1, description="Shunt susceptance [pu]"]
        G, [guess=1, description="Shunt conductance [pu]"]
    end
    @variables begin
        i_r(t), [guess=0.1, description="Current real part [pu]"]
        i_i(t), [guess=0.1, description="Current imaginary part [pu]"]
        P(t), [description="Active Power [pu]"]
        Q(t), [description="Reactive Power [pu]"]
    end
    begin
        Y = G + im*B
        iload = -Y * (terminal.u_r + im*terminal.u_i)
    end
    @equations begin
        if !allow_zero_conductance
            i_r ~ simplify(real(iload))
            i_i ~ simplify(imag(iload))
        else
            @no_simplify i_r ~ simplify(real(iload))
            @no_simplify i_i ~ simplify(imag(iload))
        end
        terminal.i_r ~ i_r
        terminal.i_i ~ i_i

        # observables
        P ~ terminal.u_r*terminal.i_r + terminal.u_i*terminal.i_i
        Q ~ terminal.u_i*terminal.i_r - terminal.u_r*terminal.i_i
    end
end

# TODO: S_b for loads?
@mtkmodel ZIPLoad begin
    @parameters begin
        Pset, [description="Active Power at operation point [pu]"]
        Qset, [description="Reactive Power at operation point [pu]"]
        Vset, [guess=1,description="Voltage at operation point [pu]"]
        KpZ, [description="Active power constant impedance fraction"]
        KqZ, [description="Reactive power constant impedance fraction"]
        KpI, [description="Active power constant current fraction"]
        KqI, [description="Reactive power constant current fraction"]
        KpC=1-KpZ-KpI, [description="Active power constant power fraction"]
        KqC=1-KqZ-KqI, [description="Reactive power constant power fraction"]
    end
    @components begin
        terminal = Terminal()
    end
    @variables begin
        Vrel(t), [description="Relative voltage magnitude"]
        P(t), [description="Active Power"]
        Q(t), [description="Reactive Power"]
    end
    @equations begin
        Vrel ~ sqrt(terminal.u_r^2 + terminal.u_i^2)/Vset
        P ~ Pset*(KpZ*Vrel^2 + KpI*Vrel + KpC)
        Q ~ Qset*(KqZ*Vrel^2 + KqI*Vrel + KqC)
        # formulate equations for i_r and i_i instead
        terminal.i_r ~  simplify(real((P + im*Q)/(terminal.u_r + im*terminal.u_i)))
        terminal.i_i ~ -simplify(imag((P + im*Q)/(terminal.u_r + im*terminal.u_i)))
    end
end

@mtkmodel ConstantCurrentLoad begin
    @parameters begin
        Iset, [guess=-1, description="Current magnitude setpoint"]
        θset, [guess=0, description="Phase offset relative to voltage"]
        ε = 1e-10, [description="Small regularization term to avoid division by zero"]
    end
    @components begin
        terminal = Terminal()
    end
    begin
        umag = sqrt(terminal.u_r^2 + terminal.u_i^2 + ε^2)
    end
    @equations begin
        terminal.i_r ~ Iset * (terminal.u_r/umag * cos(θset) - terminal.u_i/umag * sin(θset))
        terminal.i_i ~ Iset * (terminal.u_i/umag * cos(θset) + terminal.u_r/umag * sin(θset))
    end
end
