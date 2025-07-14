"""
    Terminal

A ModelingToolkit connector for electrical terminals in power system components.

Represents an electrical connection point with complex voltage and current in dq coordinates.
The terminal defines the interface between power system components like buses, lines, and machines.

# Variables
- `u_r(t)`: d-axis voltage component
- `u_i(t)`: q-axis voltage component
- `i_r(t)`: d-axis current component (flow variable)
- `i_i(t)`: q-axis current component (flow variable)

# Notes
Current variables are defined as flow variables, meaning they sum to zero at connection points
according to Kirchhoff's current law.

See also: [`BusBar`](@ref), [`LineEnd`](@ref)
"""
@connector Terminal begin
    u_r(t), [description="d-voltage"]
    u_i(t), [description="q-voltage"]
    i_r(t), [guess=0, description="d-current", connect=Flow]
    i_i(t), [guess=0, description="q-current", connect=Flow]
end

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

"""
    @named busbar = BusBar()

A ModelingToolkit model representing the physical connection point within a bus in power systems.
It represents the physical busbar where all injectors and lines attach.

Within OpPoDyn.jl, it serves as an interface between the MTK world and the NetworkDynamics world:
A MTK model containing a `BusBar` the highest level is consdered a busmodel (see [`isbusmodel`](@ref))
and describes the dynamics of an entire bus. It can be transformed in a [`VertexModel`](@extref)
by calling [`Bus`](@ref).

See also: [`Terminal`](@ref), [`MTKBus`](@ref), [`Bus`](@ref)
"""
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

"""
    LineEnd

A ModelingToolkit model representing one end of a transmission line in power systems.
It represents the physical connection point at the end of a transmission line.

Within OpPoDyn.jl, it serves as an interface between the MTK world and the NetworkDynamics world:
A MTK model containing two `LineEnd`s (named `:src` and `:dst`) at the highest level is considered a linemodel (see [`islinemodel`](@ref))
and describes the dynamics of an entire line. It can be transformed in an [`EdgeModel`](@extref) by calling [`Line`](@ref).

See also: [`Terminal`](@ref), [`MTKLine`](@ref), [`Line`](@ref)
"""
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

"""
    MTKBus(injectors...; name=:bus)

Create a ModelingToolkit bus system by connecting multiple injector components.

Constructs a bus ODESystem by connecting all provided injector components to a
central [`BusBar`](@ref). Each injector component must satisfy the injector
model interface (see [`isinjectormodel`](@ref)).

# Arguments
- `injectors...`: Variable number of injector components (generators, loads, etc.)
- `name=:bus`: Name for the resulting bus system

# Returns
- An `ODESystem` representing the complete bus with all connected injectors

```
                                 ┌────────────────────┐
                                 │MTKBus   ┌─────────┐│
                                 │        ┌┤Generator││
        ┌─────────┐   ┌────┐     │┌──────┐│└─────────┘│
MTKBus(o┤Generator│, o┤Load│) => ││BusBar├o           │
        └─────────┘   └────┘     │└──────┘│┌────┐     │
                                 │        └┤Load│     │
                                 │         └────┘     │
                                 └────────────────────┘
```

See also: [`Bus`](@ref), [`BusBar`](@ref), [`isinjectormodel`](@ref)
"""
function MTKBus(injectors...; name=:bus)
    if !all(isinjectormodel.(injectors))
        throw(ArgumentError("All components must satisfy the bus component model interface!"))
    end
    @named busbar = BusBar()
    eqs = [connect(busbar.terminal, inj.terminal) for inj in injectors]
    ODESystem(eqs, t; systems=[busbar, injectors...], name)
end

"""
    MTKLine(branches...; name=:line)

Create a ModelingToolkit line system by connecting multiple branch components.

Constructs a line ODESystem by connecting all provided branch components between
source and destination line ends in parallel. Each branch component must satisfy
the branch model interface.

# Arguments
- `branches...`: Variable number of branch components (transmission lines, transformers, etc.)
- `name=:line`: Name for the resulting line system

# Returns
- An `ODESystem` representing the complete line with all connected branches

```
                                     ┌─────────────────────────────┐
                                     │MTKLine   ┌───────┐          │
                                     │         ┌┤BranchA├┐         │
         ┌───────┐    ┌───────┐      │┌───────┐│└───────┘│┌───────┐│
MTKLine(o┤BranchA├o, o┤BranchB├o) => ││LineEnd├o         o┤LineEnd││
         └───────┘    └───────┘      │└───────┘│┌───────┐│└───────┘│
                                     │  :src   └┤BranchB├┘  :dst   │
                                     │          └───────┘          │
                                     └─────────────────────────────┘
```

See also: [`Line`](@ref), [`LineEnd`](@ref), [`isbranchmodel`](@ref)
"""
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

For example, one could create a composite injector with three subsystems:
- a generator,
- a controller, and
- a load;
which is augmented with 2 equations
- one for the measurments, and
- one for the actuation.

The returned `CompositeInjector` system has the following structure: It will
automaticially create a new terminal `t` (thus satisfing the injector interface,
see [`isinjectormodel`](@ref)) which will be connected to the terminals of the
subsystems which have a terminal (machine and load in this case).

```
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
````
"""
function CompositeInjector(systems, eqs; name=Symbol(join(ModelingToolkit.getname.(systems), "_")))
    @named terminal = Terminal()
    ivs = ModelingToolkit.get_iv.(systems)
    @assert allequal(ivs) "Systems have different independent variables! $ivs"
    iv = first(ivs)
    termeqs = [connect(sys.terminal, terminal) for sys in systems if isinjectormodel(sys)]
    ODESystem(vcat(termeqs,eqs), iv; systems=vcat(terminal, systems), name)
end

function CompositeInjector(systems...; kwargs...)
    systems = collect(systems) # tuple -> vector
    with_terminal = filter(isinjectormodel, systems)
    without_terminal = filter(x->!isinjectormodel(x), systems)

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

Ibase(S, V) = S/V
Zbase(S, V) = V^2/S
Ybase(S, V) = S/V^2
