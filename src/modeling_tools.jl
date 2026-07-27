"""
    @named t = Terminal()

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

####
#### Global per-unit base references
####
const SBASE_DEFAULT = 100.0
const ωBASE_DEFAULT = 2π * 50
const VBASE_DEFAULT = 1.0
const SBASE = Ref(SBASE_DEFAULT)
const ωBASE = Ref(ωBASE_DEFAULT)
const VBASE = Ref(VBASE_DEFAULT)

"""
    get_Sbase()

The global system power base in MVA (`100.0` unless changed). Buses and lines default their
`Sbase` to this **at construction time**, so changing it later does not affect models that
already exist. Set it with [`set_Sbase!`](@ref).
"""
get_Sbase() = SBASE[]

"""
    set_Sbase!(val=100.0)

Set the global system power base to `val` MVA and return it. Called without an argument it
restores the default, so a script that temporarily changes the power base can undo that with
a bare `set_Sbase!()`. Read it back with [`get_Sbase`](@ref).
"""
set_Sbase!(val=SBASE_DEFAULT) = (SBASE[] = float(val))

"""
    get_ωbase()

The global system frequency base in rad/s (`2π*50` unless changed). Buses and lines default
their `ωbase` to this **at construction time**. Set it with [`set_ωbase!`](@ref) or
[`set_fbase!`](@ref); see [`get_fbase`](@ref) for the same value in Hz.
"""
get_ωbase() = ωBASE[]

"""
    set_ωbase!(val=2π*50)

Set the global system frequency base to `val` rad/s and return it. Called without an argument
it restores the default (50 Hz), so a script that temporarily changes the frequency base can
undo that with a bare `set_ωbase!()`.

See also [`set_fbase!`](@ref) to set it from a frequency in Hz, and [`get_ωbase`](@ref).
"""
set_ωbase!(val=ωBASE_DEFAULT) = (ωBASE[] = float(val))

"""
    get_fbase()

The global system frequency base in Hz, i.e. [`get_ωbase`](@ref)`()/2π` (`50.0` unless
changed). Set it with [`set_fbase!`](@ref).
"""
get_fbase() = get_ωbase()/(2π)

"""
    set_fbase!(f=50)

Set the global system frequency base from a frequency `f` in Hz (i.e. `ωbase = 2π·f`) and
return the resulting angular frequency in rad/s. Called without an argument it restores the
default of 50 Hz. Read it back with [`get_fbase`](@ref).
"""
set_fbase!(f=ωBASE_DEFAULT/(2π)) = set_ωbase!(2π * float(f))

"""
    get_Vbase()

The global *fallback* voltage base in kV (`1.0` unless changed). Unlike `Sbase`/`ωbase`,
voltage base is genuinely per-bus (`Vbase` on each busbar) — this global is only the fallback
used when nothing more specific is given. Set it with [`set_Vbase!`](@ref).
"""
get_Vbase() = VBASE[]

"""
    set_Vbase!(val=1.0)

Set the global fallback voltage base to `val` kV and return it. Called without an argument it
restores the default. Read it back with [`get_Vbase`](@ref).
"""
set_Vbase!(val=VBASE_DEFAULT) = (VBASE[] = float(val))

"""
    @named systembase = SystemBase()
    SystemBase(; name, Sbase=get_Sbase(), ωbase=get_ωbase(), ωframe=1.0)

The canonical home of the *global* per-unit bases and the global reference frame. A
parameter-only component, instantiated **exactly once per bus and once per line** by
[`MTKBus`](@ref) / [`MTKLine`](@ref) as a top-level sibling of the `busbar` (resp. of the
two [`LineEnd`](@ref)s), under the instance name `systembase`.

Every component that needs a global base declares a local shadow bound to it:

```julia
@parameters begin
    Sbase, [bound_to = :systembase₊Sbase]
    ωbase, [bound_to = :systembase₊ωbase]
end
```

Hand-rolled bus/line `System`s need not add a `SystemBase` themselves:
[`compile_bus`](@ref) / [`compile_line`](@ref) inject one if it is missing. Adding it
explicitly is still useful when you want to set the bases on the symbolic model
(`SystemBase(ωbase=2π*60)`) rather than on the compiled component.

# Parameters
- `Sbase`: system power base [MVA], defaults to the global [`get_Sbase`](@ref)`()`
- `ωbase`: system frequency base [rad/s], defaults to the global [`get_ωbase`](@ref)`()`
- `ωframe`: speed of the global dq reference frame [pu]. **Pinned to `1` in 5.0** — it is a
  gauge, present so that frame-dependent terms (`ωbase*(ω - ωframe)`, cross-coupling in EMT
  branches) name it explicitly instead of hiding an invisible `1`. Do not change it.

See also: [`BusBar`](@ref), [`LineEnd`](@ref), [`set_Sbase!`](@ref), [`set_ωbase!`](@ref)
"""
@component function SystemBase(; name, defaults...)
    params = @parameters begin
        Sbase = get_Sbase(), [description="System power base [MVA]"]
        ωbase = get_ωbase(), [description="System frequency base [rad/s]"]
        ωframe = 1.0,    [description="Global dq frame speed [pu] (pinned to 1 in PD@v5.0)"]
    end
    vars = @variables begin
        fbase(t), [description="System frequency [Hz] (ωbase/2π)"]
    end
    eqs = Equation[
        fbase ~ ωbase/(2π)
    ]
    set_mtk_defaults(System(eqs, t, vars, params; name), defaults)
end

"""
    add_systembase(sys::System)

Return `sys` with a [`SystemBase`](@ref) named `systembase` added at the top level, or `sys`
unchanged if it already has one.
"""
function add_systembase(sys::System)
    has_systembase(sys) && return sys
    @named systembase = SystemBase()
    @set sys.systems = vcat(ModelingToolkitBase.get_systems(sys), systembase)
end

"""
    has_systembase(sys::System)

Whether `sys` has a top-level subsystem named `systembase`. See `add_systembase`.
Internal; not exported.
"""
function has_systembase(sys::System)
    any(s -> getname(s) === :systembase, ModelingToolkitBase.get_systems(sys))
end

"""
    BusBase(; name, Vbase=get_Vbase())

The bare bus interface (the ND-facing input/output side of a bus): complex voltage
`u_r`/`u_i` (outputs), complex current `i_r`/`i_i` (inputs), pu power observables, and the
per-unit bases plus SI observables.
"""
@component function BusBase(; name, defaults...)
    params = @parameters begin
        # get_Vbase() via initf instead of default so we can distinguish between
        # soft default from global fallback vs hard default from user.
        # compile_bus special treats thes initf: directly applies for voltage
        # sources and drops for satellite buses
        Vbase, [initf_weak = get_Vbase(), description="Bus voltage base [kV]"]
        # global bases: structural aliases of the bus's `systembase` sibling
        Sbase, [bound_to = :systembase₊Sbase, description="System power base [MVA]"]
        ωbase, [bound_to = :systembase₊ωbase, description="System frequency base [rad/s]"]
        # auto demoted to observables (bound params)
        Ibase = Sbase/Vbase,   [description="Current base [kA] (Sbase/Vbase)"]
        Zbase = Vbase^2/Sbase, [description="Impedance base [Ω] (Vbase²/Sbase)"]
        Ybase = Sbase/Vbase^2, [description="Admittance base [S] (Sbase/Vbase²)"]
    end
    vars = @variables begin
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
        # SI observables (physical units, derived from pu quantities and the bus bases)
        u_kV(t), [description="bus voltage magnitude [kV]"]
        P_MW(t), [description="bus active power into network [MW]"]
        Q_MVAr(t), [description="bus reactive power into network [MVAr]"]
        i_kA(t), [description="bus current magnitude [kA]"]
        # ω(t), [description="bus angular frequency"]
    end
    eqs = [
        #observed equations
        # attension: flipped sign in P and Q, flow direction opposite to i
        P ~ u_r * (-i_r) + u_i * (-i_i)
        Q ~ u_i * (-i_r) - u_r * (-i_i)
        u_mag ~ sqrt(u_r^2 + u_i^2)
        u_arg ~ atan(u_i, u_r)
        i_mag ~ sqrt(i_r^2 + i_i^2)
        i_arg ~ atan(i_i, i_r)
        # SI observables: u_mag is pu on Vbase, P/Q pu on Sbase, i_mag pu on Ibase
        u_kV ~ u_mag * Vbase
        P_MW ~ P * Sbase
        Q_MVAr ~ Q * Sbase
        i_kA ~ i_mag * Ibase
        # ω ~ Dt(u_arg) # this can lead to Dt(i_r) and Dt(i_i) in the rhs of the equations
    ]
    set_mtk_defaults(System(eqs, t, vars, params; name), defaults)
end

"""
    @named busbar = BusBar()

A ModelingToolkit model representing the physical connection point within a bus in power systems.
It represents the physical busbar where all injectors and lines attach.

Within PowerDynamics.jl, it serves as an interface between the MTK world and the
NetworkDynamics world: A MTK model containing a `BusBar` the highest level is
consdered a busmodel (see [`isbusmodel`](@ref)) and describes the dynamics of an
entire bus. It can be transformed in a [`VertexModel`](@extref
NetworkDynamics.VertexModel-Tuple{}) by calling [`compile_bus`](@ref).

See also: [`Terminal`](@ref), [`MTKBus`](@ref), [`compile_bus`](@ref)
"""
@component function BusBar(; name, defaults...)
    @named base = BusBase()
    @named terminal = Terminal()
    @unpack u_r, u_i, i_r, i_i = base
    eqs = [
        u_r ~ terminal.u_r
        u_i ~ terminal.u_i
        i_r ~ terminal.i_r
        i_i ~ terminal.i_i
    ]
    sys = extend(System(eqs, t; name, systems=[terminal]), base)
    set_mtk_defaults(sys, defaults)
end


"""
    LineEnd(; name, side, defaults...)

A ModelingToolkit model representing one end of a transmission line in power systems.
It represents the physical connection point at the end of a transmission line.

Within PowerDynamics.jl, it serves as an interface between the MTK world and the
NetworkDynamics world: A MTK model containing two `LineEnd`s (named `:src` and
`:dst`) at the highest level is considered a linemodel (see
[`islinemodel`](@ref)) and describes the dynamics of an entire line. It can be
transformed in an [`EdgeModel`](@extref NetworkDynamics.EdgeModel-Tuple{}) by
calling [`compile_line`](@ref).

`side=:src/:dst` determines where the line end gets inherits it `Vbase` from.

See also: [`Terminal`](@ref), [`MTKLine`](@ref), [`compile_line`](@ref)
"""
@component function LineEnd(; name, side, defaults...)
    side in (:src, :dst) ||
        throw(ArgumentError("LineEnd: `side` must be :src or :dst, got $(repr(side)). \
                             It selects which incident bus this end inherits `Vbase` from."))
    @named terminal = Terminal()
    params = @parameters begin
        # inherit from Vbase from incident bus via default_from
        Vbase, [default_from = (side, :busbar₊Vbase), description="Line-end voltage base [kV]"]
        # global bases: structural aliases of the line's `systembase` sibling
        Sbase, [bound_to = :systembase₊Sbase, description="System power base [MVA]"]
        ωbase, [bound_to = :systembase₊ωbase, description="System frequency base [rad/s]"]
        # Derived bases (parameter bindings -> observables); see BusBase.
        Ibase = Sbase/Vbase,   [description="Current base [kA] (Sbase/Vbase)"]
        Zbase = Vbase^2/Sbase, [description="Impedance base [Ω] (Vbase²/Sbase)"]
        Ybase = Sbase/Vbase^2, [description="Admittance base [S] (Sbase/Vbase²)"]
    end
    vars = @variables begin
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
        # SI observables (physical units, derived from pu quantities and the end bases)
        u_kV(t), [description="line end voltage magnitude [kV]"]
        P_MW(t), [description="line end active power [MW]"]
        Q_MVAr(t), [description="line end reactive power [MVAr]"]
        i_kA(t), [description="line end current magnitude [kA]"]
    end
    eqs = [
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
        # SI observables
        u_kV ~ u_mag * Vbase
        P_MW ~ P * Sbase
        Q_MVAr ~ Q * Sbase
        i_kA ~ i_mag * Ibase
    ]
    set_mtk_defaults(System(eqs, t, vars, params; name, systems=[terminal]), defaults)
end

"""
    MTKBus(injectors...; name=:bus)

Create a ModelingToolkit bus system by connecting multiple injector components.

Constructs a bus `System` by connecting all provided injector components to a
central [`BusBar`](@ref). Each injector component must satisfy the injector
model interface (see [`isinjectormodel`](@ref)).

Additionally injects a [`SystemBase`](@ref) component to provide the global
per-unit bases and reference frame.

# Arguments
- `injectors...`: Variable number of injector components (generators, loads, etc.)
- `name=:bus`: Name for the resulting bus system

# Returns
- An `System` representing the complete bus with all connected injectors

```asciiart
                                 ┌────────────────────┐
                                 │MTKBus   ┌─────────┐│
                                 │┌──────┐┌┤Generator││
        ┌─────────┐   ┌────┐     ││BusBar├o└─────────┘│
MTKBus(o┤Generator│, o┤Load│) => │└──────┘│┌────┐     │
        └─────────┘   └────┘     │        └┤Load│     │
                                 │         └────┘     │
                                 │ + SystemBase       │
                                 └────────────────────┘
```

See also: [`compile_bus`](@ref), [`BusBar`](@ref), [`isinjectormodel`](@ref)
"""
function MTKBus(injectors...; name=:bus, Vbase=nothing)
    if !all(isinjectormodel.(injectors))
        throw(ArgumentError("All components must satisfy the injector model interface!"))
    end
    @named busbar = BusBar(; _busbar_defaults(Vbase)...)
    @named systembase = SystemBase()
    eqs = [connect(busbar.terminal, inj.terminal) for inj in injectors]
    System(eqs, t; systems=[busbar, systembase, injectors...], name)
end

function MTKBus(systems::Union{AbstractVector,Tuple,Set}, eqs=autoconnections(systems); name=:bus, Vbase=nothing)
    @named busbar = BusBar(; _busbar_defaults(Vbase)...)
    @named systembase = SystemBase()
    for sys in systems
        if isinjectormodel(sys)
            eq = connect(sys.terminal, busbar.terminal)
            push!(eqs, eq)
        end
    end
    System(eqs, t; systems=vcat(busbar, systembase, systems), name)
end

_busbar_defaults(Vbase) = isnothing(Vbase) ? (;) : (; Vbase)

"""
    MTKLine(branches...; name=:line)

Create a ModelingToolkit line system by connecting multiple branch components.

Constructs a line `System` by connecting all provided branch components between
source and destination line ends in parallel. Each branch component must satisfy
the branch model interface.

Additionally injects a [`SystemBase`](@ref) component to provide the global
per-unit bases and reference frame.

# Arguments
- `branches...`: Variable number of branch components (transmission lines, transformers, etc.)
- `name=:line`: Name for the resulting line system

# Returns
- An `System` representing the complete line with all connected branches

```asciiart
                                     ┌─────────────────────────────┐
                                     │MTKLine   ┌───────┐          │
                                     │         ┌┤BranchA├┐         │
         ┌───────┐    ┌───────┐      │┌───────┐│└───────┘│┌───────┐│
MTKLine(o┤BranchA├o, o┤BranchB├o) => ││LineEnd├o         o┤LineEnd││
         └───────┘    └───────┘      │└───────┘│┌───────┐│└───────┘│
                                     │  :src   └┤BranchB├┘  :dst   │
                                     │          └───────┘          │
                                     │ + SystemBase                │
                                     └─────────────────────────────┘
```

See also: [`compile_line`](@ref), [`LineEnd`](@ref), [`isbranchmodel`](@ref)
"""
function MTKLine(branches...; name=:line)
    if !all(isbranchmodel.(branches))
        throw(ArgumentError("All components must satisfy the branch model interface!"))
    end
    systems = @named begin
        src = LineEnd(; side=:src)
        dst = LineEnd(; side=:dst)
        systembase = SystemBase()
    end

    eqs = [[connect(src.terminal, branch.src) for branch in branches]...,
           [connect(dst.terminal, branch.dst) for branch in branches]...]

    System(eqs, t; systems=[ systems..., branches...], name)
end


# this is a hack to convice MTK that i_r and i_i do depend on u_r and u_i
# we tell MTK to not further resolve, which makes it accept
# the curent constraint as a valid constraint for ur/ui
@mtkmodel KirchoffBus begin
    @components begin
        busbar = BusBase()
        systembase = SystemBase()
    end
    @equations begin
        busbar.i_r ~ implicit_output(busbar.u_r)
        busbar.i_i ~ implicit_output(busbar.u_i)
    end
end

MTKBus(; name=:bus) = KirchoffBus(; name)

"""
    CompositeInjector(systems, eqs=autoconnections(systems); name=Symbol(join(getname.(systems), "_")))

Create an injector object which contains several subsystems. Every subsystem which has a `terminal` will be connected
to a newly created terminal of the composite injector. The subsystems are namespaced within the composite injector.

There are two options for additional connections between the subsystems:
- interconnections will be created automatically using some name-matching heuristics using `autoconnections(systems)`:
  It searches all `Blocks.RealOutput` and `Blocks.RealInput`, and tries to find a single matching output for each input.
- alternatively pass connecting equations of the form `[connect(sys1.output, sys2.input)]` explicitly

For example, one could create a composite injector with three subsystems:
- a generator,
- a controller, and
- a load;
which is augmented with 2 connection equations
- one for the measurements (generator -> controller), and
- one for the actuation (controller -> generator).

The returned model contains a new terminal `:terminal` at the toplevel, thus
satisfying the injector interface, see [`isinjectormodel`](@ref)). It can be used
as such in the [`MTKBus`](@ref) constructor.
```asciiart
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
```
"""
function CompositeInjector(systems, eqs=autoconnections(systems); name=Symbol(join(ModelingToolkitBase.getname.(systems), "_")))
    @named terminal = Terminal()
    ivs = ModelingToolkitBase.get_iv.(systems)
    @assert allequal(ivs) "Systems have different independent variables! $ivs"
    iv = first(ivs)
    termeqs = Equation[connect(sys.terminal, terminal) for sys in systems if isinjectormodel(sys)]
    System(vcat(termeqs, eqs), iv; systems=vcat(terminal, systems), name)
end

function autoconnections(systems)
    systems = collect(systems) # tuple -> vector

    outputs = mapreduce(vcat, systems) do sys
        subouts = filter(ModelingToolkitBase.get_systems(sys)) do subsys
            contains(repr(ModelingToolkitBase.get_gui_metadata(subsys).type), "Blocks.RealOutput")
        end
        subout_names = getname.(subouts)
        subouts_resolved = getproperty.(Ref(sys), subout_names)
        subout_names .=> subouts_resolved
    end
    inputs = mapreduce(vcat, systems) do sys
        subins = filter(ModelingToolkitBase.get_systems(sys)) do subsys
            contains(repr(ModelingToolkitBase.get_gui_metadata(subsys).type), "Blocks.RealInput")
        end
        subin_names = getname.(subins)
        subins_resolved = getproperty.(Ref(sys), subin_names)
        subin_names .=> subins_resolved
    end
    out_dict = Dict(outputs...)
    in_dict = Dict(inputs...)

    outnames = collect(keys(out_dict))
    eqs = Equation[]
    for (iname, isys) in in_dict
        out = _findmatch(iname, outnames)
        push!(eqs, connect(out_dict[out], isys))
    end
    return eqs
end
function _findmatch(in_symbol, outs_symbol)
    in = string(in_symbol)
    outs = string.(outs_symbol)

    # try with full name
    idx = findall(isequal(in), outs)
    if isempty(idx)
        # try to strip common pre and suffixes
        inshort = replace(in, r"(_meas|meas|_in|in)$" => "")
        outsshort = replace.(outs, r"(_meas|meas|_out|out)$" => "")

        idx = findall(isequal(inshort), outsshort)
    end

    if isempty(idx)
        error("Could not find a matching output for input :$in in $(Symbol.(outs))")
    elseif length(idx) == 1
        return outs_symbol[only(idx)]
    else
        error("Multiple possible matches found for input :$in. Candidates: $(Symbol.(outs[idx]))")
    end
end

"""
    Ibase(S, V)

Calculates current pu base based on Sbase and Vbase: Ibase = Sbase/Vbase.
"""
Ibase(S, V) = S/V

"""
    Zbase(S, V)

Calculates impedance pu base based on Sbase and Vbase: Zbase = Vbase²/Sbase.
"""
Zbase(S, V) = V^2/S

"""
    Ybase(S, V)

Calculates admittance pu base based on Sbase and Vbase: Ybase = Sbase/Vbase².
"""
Ybase(S, V) = S/V^2
