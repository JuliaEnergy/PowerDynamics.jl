@mtkmodel LineEnd begin
    @components begin
        terminal = Terminal()
    end
    @variables begin
        u_r(t), [description="line end d-voltage", input=true]
        u_i(t), [description="line end q-voltage", input=true]
        i_r(t), [description="line end d-current", output=true]
        i_i(t), [description="line end d-current", output=true]
        P(t), [description="line end active power"]
        Q(t), [description="line end reactive power"]
        u_mag(t), [description="line end voltage magnitude"]
        u_arg(t), [description="line end voltage argument"]
        i_mag(t), [description="line end current magnitude"]
        i_arg(t), [description="line end current argument"]
    end
    @equations begin
        u_r ~  terminal.u_r
        u_i ~  terminal.u_i
        i_r ~ -terminal.i_r
        i_i ~ -terminal.i_i
        # observed equations
        P ~ u_r * i_r + u_i * i_i
        Q ~ u_i * i_r - u_r * i_i
        u_mag ~ sqrt(u_r^2 + u_i^2)
        u_arg ~ atan(u_i, u_r)
        i_mag ~ sqrt(i_r^2 + i_i^2)
        i_arg ~ atan(i_i, i_r)
    end
end

function MTKLine(branches...; name=:line)
    if !all(isbranchmodel.(branches))
        throw(ArgumentError("All components must satisfy the branch model interface!"))
    end
    systems = @named begin
        src = LineEnd()
        dst = LineEnd()
    end

    eqs = [[connect(src.terminal, branch.src) for branch in branches]...,
           [connect(dst.terminal, branch.dst) for branch in branches]...]

    ODESystem(eqs, t; systems=[ systems..., branches...], name)
end
