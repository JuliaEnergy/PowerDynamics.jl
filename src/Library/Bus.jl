@mtkmodel BusBase begin
    @variables begin
        u_r(t)=1, [description="bus d-voltage", output=true]
        u_i(t)=0, [description="bus q-voltage", output=true]
        i_r(t), [description="bus d-current", input=true]
        i_i(t), [description="bus d-current", input=true]
        P(t), [description="bus active power"]
        Q(t), [description="bus reactive power"]
    end
    @equations begin
        P ~ u_r * i_r + u_i * i_i
        Q ~ u_i * i_r - u_r * i_i
    end
end

@mtkmodel BusBar begin
    @extend BusBase()
    @components begin
        terminal = Terminal()
    end
    @equations begin
        u_r ~ terminal.u_r
        u_i ~ terminal.u_i
        i_r ~ -terminal.i_r
        i_i ~ -terminal.i_i
    end
end

@mtkmodel SlackAlgebraic begin
    @components begin
        busbar = BusBase()
    end
    @parameters begin
        u_set_r=1, [description="bus d-voltage setpoint"]
        u_set_i=0, [description="bus q-voltage setpoint"]
    end
    @equations begin
        busbar.u_r ~ u_set_r
        busbar.u_i ~ u_set_i
    end
end

@mtkmodel SlackDifferential begin
    @parameters begin
        u_init_r=1, [description="bus d-voltage initial value"]
        u_init_i=0, [description="bus q-voltage initial value"]
    end
    @components begin
        busbar = BusBase(;u_r=u_init_r, u_i=u_init_i)
    end
    @equations begin
        Dt(busbar.u_r) ~ 0
        Dt(busbar.u_i) ~ 0
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
