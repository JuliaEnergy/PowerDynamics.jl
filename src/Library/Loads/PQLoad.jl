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
        Pset, [description="Active Power demand"]
        Qset, [description="Reactive Power demand"]
        αP, [description="Active Power exponent"]
        αQ, [description="Reactive Power exponent"]
        Vn, [description="Nominal voltage (where real power equals set power)"]
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
    @components begin
        terminal = Terminal()
    end
    @parameters begin
        Pset, [description="Active Power demand [pu]"]
        Qset, [description="Reactive Power demand [pu]"]
        Vset, [guess=1,description="Nominal voltage [pu]"]
    end
    @variables begin
        P(t), [description="Active Power [pu]"]
        Q(t), [description="Reactive Power [pu]"]
    end
    begin
        S = Pset + im*Qset
        Y = conj(S)/Vset^2
        iload = Y * (terminal.u_r + im*terminal.u_i)
    end
    @equations begin
        terminal.i_r ~ simplify(real(iload))
        terminal.i_i ~ simplify(imag(iload))

        # observables
        P ~ terminal.u_r*terminal.i_r + terminal.u_i*terminal.i_i
        Q ~ terminal.u_i*terminal.i_r - terminal.u_r*terminal.i_i
    end

end
