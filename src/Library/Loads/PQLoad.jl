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
