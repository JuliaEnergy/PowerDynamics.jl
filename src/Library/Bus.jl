@mtkmodel BusBase begin
    @variables begin
        u_r(t)=1, [description="bus d-voltage", output=true]
        u_i(t)=0, [description="bus q-voltage", output=true]
        i_r(t), [description="bus d-current (flowing into bus)", input=true]
        i_i(t), [description="bus d-current (flowing into bus)", input=true]
        P(t), [description="bus active power (flowing into network)"]
        Q(t), [description="bus reactive power (flowing into network)"]
        u_mag(t), [description="bus voltage magnitude"]
        u_arg(t), [description="bus voltage argument"]
        i_mag(t), [description="bus current magnitude"]
        i_arg(t), [description="bus current argument"]
        # ω(t), [description="bus angular frequency"]
    end
    @equations begin
        #observed equations
        # attension: flipped sign in P and Q, flow direction opposite to i
        P ~ u_r * (-i_r) + u_i * (-i_i)
        Q ~ u_i * (-i_r) - u_r * (-i_i)
        u_mag ~ sqrt(u_r^2 + u_i^2)
        u_arg ~ atan(u_i, u_r)
        i_mag ~ sqrt(i_r^2 + i_i^2)
        i_arg ~ atan(i_i, i_r)
        # ω ~ Dt(u_arg) # this can lead to Dt(i_r) and Dt(i_i) in the rhs of the equations
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
        i_r ~ terminal.i_r
        i_i ~ terminal.i_i
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

function MTKBus(injectors...; name=:bus)
    if !all(iscomponentmodel.(injectors))
        throw(ArgumentError("All components must satisfy the bus component model interface!"))
    end
    @named busbar = BusBar()
    eqs = [connect(busbar.terminal, inj.terminal) for inj in injectors]
    ODESystem(eqs, t; systems=[busbar, injectors...], name)
end


# this is a hack to convice MTK that i_r and i_i do depend on u_r and u_i
# we tell MTK to not further resolve, which makes it accept
# the curent constraint as a valid constraint for ur/ui
_to_zero(x) = 0.0
ModelingToolkit.@register_symbolic _to_zero(x)::Float64
@mtkmodel KirchoffBus begin
    @components begin
        busbar = BusBase()
    end
    @equations begin
        busbar.i_r ~ _to_zero(busbar.u_r)
        busbar.i_i ~ _to_zero(busbar.u_i)
    end
end

MTKBus(; name=:bus) = KirchoffBus(; name)

"""
    CompositeInjector(systems, eqs; name=Symbol(join(getname.(systems), "_")))

Create a injector object which contains several subsystems. Every subsystem which has a `terminal` will be connected
to a newly created terminal of the composite injector. Can contain further systems, such as controllers, with the
additional connection equations `eqs`.

For example, one could create a composite injector with three subsustens:
- a generator,
- a controller, and
- a load;
which is augmented with 2 equations
- one for the measurments, and
- one for the actuation.

The returned `CompositeInjector` system has the following structure:
It will automaticially create a new terminal `t` (thus satisfing the [Injector
Interface](@ref)) which will be connected to the terminals of the subsystems
which have a terminal (machine and load in this case).

        ┌────────────────────────────────────┐
        │ CompositeInjector                  │
        │              ╭───→───╮ measurements│
        │    ┌─────────┴─┐   ┌─┴──────────┐  │
    (t) │  o─┤ Generator │   │ Controller │  │
     o──┼──┤ └─────────┬─┘   └─┬──────────┘  │
        │  │           ╰───←───╯ actuation   │
        │  │ ┌──────┐                        │
        │  o─┤ Load │                        │
        │    └──────┘                        │
        └────────────────────────────────────┘
"""
function CompositeInjector(systems, eqs; name=Symbol(join(ModelingToolkit.getname.(systems), "_")))
    @named terminal = Terminal()
    ivs = ModelingToolkit.get_iv.(systems)
    @assert allequal(ivs) "Systems have different independent variables! $ivs"
    iv = first(ivs)
    termeqs = [connect(sys.terminal, terminal) for sys in systems if iscomponentmodel(sys)]
    ODESystem(vcat(termeqs,eqs), iv; systems=vcat(terminal, systems), name)
end

function CompositeInjector(systems...; kwargs...)
    systems = collect(systems) # tuple -> vector
    with_terminal = filter(iscomponentmodel, systems)
    without_terminal = filter(x->!iscomponentmodel(x), systems)

    outputs = mapreduce(vcat, systems) do sys
        subouts = filter(ModelingToolkit.get_systems(sys)) do subsys
            contains(repr(ModelingToolkit.get_gui_metadata(subsys).type), "Blocks.RealOutput")
        end
        subout_names = getname.(subouts)
        subouts_resolved = getproperty.(Ref(sys), subout_names)
        subout_names .=> subouts_resolved
    end
    inputs = mapreduce(vcat, systems) do sys
        subins = filter(ModelingToolkit.get_systems(sys)) do subsys
            contains(repr(ModelingToolkit.get_gui_metadata(subsys).type), "Blocks.RealInput")
        end
        subin_names = getname.(subins)
        subins_resolved = getproperty.(Ref(sys), subin_names)
        subin_names .=> subins_resolved
    end
    out_dict = Dict(outputs...)
    in_dict = Dict(inputs...)

    # FIXME: hard coded for now
    outnames = collect(keys(out_dict))
    eqs = []
    for (iname, isys) in in_dict
        out = findmatch(iname, outnames)
        push!(eqs, connect(out_dict[out], isys))
    end
    CompositeInjector(systems, eqs; kwargs...)
end

function findmatch(in_symbol, outs_symbol)
    in = string(in_symbol)
    outs = string.(outs_symbol)

    # try with full name
    idx = findall(isequal(in), outs)
    if isempty(idx)
        # try to strip comon pre and suffixes
        inshort = replace(in, r"(_meas|meas|_in|in)$" => "")
        outsshort = replace.(outs, r"(_meas|meas|_out|out)$" => "")

        idx = findall(isequal(inshort), outsshort)
    end

    if isempty(idx)
        error("Could not find a match for $in in $outs")
    elseif length(idx) == 1
        return outs_symbol[only(idx)]
    else
        error("Multiple matches found for $in in $outs")
    end
end
