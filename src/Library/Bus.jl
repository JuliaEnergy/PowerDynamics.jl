@mtkmodel BusBar begin
    @components begin
        terminal = Terminal()
    end
    @variables begin
        u_r(t), [description="bus d-voltage", output=true]
        u_i(t), [description="bus q-voltage", output=true]
        i_r(t), [description="bus d-current", input=true]
        i_i(t), [description="bus d-current", input=true]
    end
    @equations begin
        u_r ~ terminal.u_r
        u_i ~ terminal.u_i
        i_r ~ terminal.i_r
        i_i ~ terminal.i_i
    end
end

function Bus(_injectors...; name=:bus)
    injectors = Tuple[]
    for injector in _injectors
        if injector isa Tuple
            push!(injectors, injector)
        else
            push!(injectors, (injector, injector.terminal))
        end
    end
    @named busbar = BusBar()
    eqs = [connect(busbar.terminal, inj[2]) for inj in injectors]
    ODESystem(eqs, t; systems=[busbar, first.(injectors)...], name)
end
