@mtkmodel LineEnd begin
    @components begin
        terminal = Terminal()
    end
    @variables begin
        u_r(t), [description="line end d-voltage", input=true]
        u_i(t), [description="line end q-voltage", input=true]
        i_r(t), [description="line end d-current", output=true]
        i_i(t), [description="line end d-current", output=true]
    end
    @equations begin
        u_r ~  terminal.u_r
        u_i ~  terminal.u_i
        i_r ~ -terminal.i_r
        i_i ~ -terminal.i_i
    end
end

function LineModel(branches...; name=:line)
    if !all(isbranchmodel.(b))
        throw(ArgumentError("All components must satisfy the branch model interface!"))
    end
    systems = @named begin
        src = LineEnd()
        dst = LineEnd()
    end

    eqs = [[connect(src.terminal, branch.src) for branch in branches]...,
           [connect(dst.terminal, branch.dst) for branch in branches]...]

    ODESystem(eqs, t; systems=[ systems..., first.(branches)...], name)
end
