# Modeling Concepts

## Terminal

The `Terminal`─Connector is an important building block for every model.
It represents a connection point with constant voltage in dq─cordinates `u_r` and `u_i` and enforces the kirchoff constraints `sum(i_r)=0` and `sum(i_i)=0`.

## Modeling of Buses
### Model class `Injector`

An injector is a class of components with a single `Terminal()` (called `:terminal`).
Examples for injectors might be Generators, Shunts, Loads.
```
      ┌───────────┐
(t)   │           │
 o←───┤  Injector │
      │           │
      └───────────┘
```
The current for injectors is always in injector convention, i.e. positive currents flow *out* of the injector *towards* the terminal.

!!! note "Model classes"
    Model "classes" are nothing formalized. In this document, a model class is just a description for some `ODESystem` from `ModelingToolkit.jl`, which satisfies certain requirements.
    For example, any `ODESystem` is considered an "Injector" if it contains a connector `Terminal()` called `:terminal`.

!!! details "Code example: definition of PQ load as injector"
    ```@example concepts
    using PowerDynamics, PowerDynamics.Library, ModelingToolkit
    using ModelingToolkit: D_nounits as Dt, t_nounits as t
    @mtkmodel MyPQLoad begin
        @components begin
            terminal = Terminal()
        end
        @parameters begin
            Pset, [description="Active Power demand"]
            Qset, [description="Reactive Power demand"]
        end
        @variables begin
            P(t), [description="Active Power"]
            Q(t), [description="Reactive Power"]
        end
        @equations begin
            P ~ terminal.u_r*terminal.i_r + terminal.u_i*terminal.i_i
            Q ~ terminal.u_i*terminal.i_r - terminal.u_r*terminal.i_i
            # if possible, its better for the solver to explicitly provide algebraic equations for the current
            terminal.i_r ~ (Pset*terminal.u_r + Qset*terminal.u_i)/(terminal.u_r^2 + terminal.u_i^2)
            terminal.i_i ~ (Pset*terminal.u_i - Qset*terminal.u_r)/(terminal.u_r^2 + terminal.u_i^2)
        end
    end
    MyPQLoad(name=:pqload) #hide
    nothing #hide
    ```

### Model class `MTKBus`
A `MTKBus` isa class of models, which are used to describe the dynamic behavior of a full bus in a power grid.
Each `MTKBus` musst contain a predefined model of type `BusBar()` (named `:busbar`).
This busbar represents the connection point to the grid.
Optionally, it may contain various injectors.

```
 ┌───────────────────────────────────┐
 │ MTKBus             ┌───────────┐  │
 │  ┌──────────┐   ┌──┤ Generator │  │
 │  │          │   │  └───────────┘  │
 │  │  BusBar  ├───o                 │
 │  │          │   │  ┌───────────┐  │
 │  └──────────┘   └──┤ Load      │  │
 │                    └───────────┘  │
 └───────────────────────────────────┘
```
Sometimes it is not possible to connect all injectors directly but instead one needs or wants `Branches` between the busbar and injector terminal.
As long as the `:busbar` is present at the toplevel, there are few limitations on the overall model complexity.

For simple models (direct connections of a few injectors) it is possible to use the convenience method `MTKBus(injectors...)` to create the composite model based on provide injector models.

!!! details "Code example: definition of a Bus containing a swing equation and a load"
    ```@example concepts
    using PowerDynamics, PowerDynamics.Library, ModelingToolkit
    @mtkmodel MyMTKBus begin
        @components begin
            busbar = BusBar()
            swing = Swing()
            load = PQLoad()
        end
        @equations begin
            connect(busbar.terminal, swing.terminal)
            connect(busbar.terminal, load.terminal)
        end
    end
    MyMTKBus(name=:bus) #hide
    nothing #hide
    ```
    Alternativly, for that system you could have just called
    ```@example concepts
    mybus = MTKBus(Swing(;name=:swing), PQLoad(;name=:load))
    nothing #hide
    ```
    to get an instance of a model which is structually aquivalent to `MyMTKBus`.

## Line Modeling
### Model class `Branch`
A branch is the two-port equivalent to an injector.
I needs to have two `Terminal()`s, one is called `:src`, the other `:dst`.

Examples for branches are: PI─Model branches, dynamic RL branches or transformers.
```
      ┌───────────┐
(src) │           │ (dst)
  o←──┤  Branch   ├──→o
      │           │
      └───────────┘
```
*Both* ends follow the injector interface, i.e. current leaving the device towards the terminals is always positive.

!!! details "Code example: algebraic R-line" 
    ```@example conceps
    using PowerDynamics, PowerDynamics.Library, ModelingToolkit
    @mtkmodel MyRLine begin
        @components begin
            src = Terminal()
            dst = Terminal()
        end
        @parameters begin
            R=0, [description="Resistance"]
        end
        @equations begin
            dst.i_r ~ (dst.u_r - src.u_r)/R
            dst.i_i ~ (dst.u_i - src.u_i)/R
            src.i_r ~ -dst.i_r
            src.i_i ~ -dst.i_i
        end
    end
    MyRLine(name=:rline) #hide
    nothing #hide
    ```

### Model class: `MTKLine`
Similar to the `MTKBus`, a `MTKLine` is a model class which represents a transmission line in the network.

It musst contain two `LineEnd()` instances, one called `:src`, one called `:dst`.

```
 ┌────────────────────────────────────────────────┐
 │ MTKLine          ┌──────────┐                  │
 │  ┌─────────┐  ┌──┤ Branch A │──┐  ┌─────────┐  │
 │  │ LineEnd │  │  └──────────┘  │  │ LineEnd │  │
 │  │  :src   ├──o                o──┤  :dst   │  │
 │  │         │  │  ┌──────────┐  │  │         │  │
 │  └─────────┘  └──┤ Branch B │──┘  └─────────┘  │
 │                  └──────────┘                  │
 └────────────────────────────────────────────────┘
```

Simple line models, which consist only of valid `Branch` models can be instantiated using the `MTKLine(branches...)` constructor.

More complex models can be created manually.
For example if you want to chain multiple branches between the `LineEnds`, for example something like

```
LineEnd(:src) ──o── Transformer ──o── Pi─Line ──o── LineEnd(:dst)
```

!!! details "Code example: Transmission line with two pi-branches"
    ```@example concepts
    using PowerDynamics, PowerDynamics.Library, ModelingToolkit
    @mtkmodel MyMTKLine begin
        @components begin
            src = LineEnd()
            dst = LineEnd()
            branch1 = PiLine()
            branch2 = PiLine()
        end
        @equations begin
            connect(src.terminal, branch1.src)
            connect(src.terminal, branch2.src)
            connect(dst.terminal, branch1.dst)
            connect(dst.terminal, branch2.dst)
        end
    end
    MyMTKLine(name=:mtkline) #hide
    nothing #hide
    ```
    Alternatively, an equivalent model with multiple valid branch models in parallel could be created and instantiated with the convenience constructor
    ```@example concepts
    line = MTKLine(PiLine(;name=:branch1), PiLine(;name=:branch2))
    nothing #hide
    ```


## From MTK Models to NetworkDynamics
Valid `MTKLine` and `MTKBus` can be transformed into so called `Line` and `Bus` objects.

`Line` and `Bus` structs are no MTK models anymore, but rather containers.
Currently, they mainly contain a NetworkDynamics component function (`VertexModel`, `EdgeModel`).

Eventually, those models will contain more metadata. For example

 - static representation for powerflow,
 - possibly local information about PU system (for transforming parameters between SI/PU),
 - meta information for initialization, for example initialization model or the information which parameters are considered "tunable" in order to initialize the dynamical model

The exact structure here is not clear yet!

The result would look something like that:
```@example concepts
using PowerDynamics, PowerDynamics.Library, ModelingToolkit
using Graphs, NetworkDynamics
using OrdinaryDiffEqRosenbrock, OrdinaryDiffEqNonlinearSolve
using CairoMakie
nothing #hide
```

Define a swing bus with load
```@example concepts
# define injectors
@named swing = Swing(; Pm=1, V=1, D=0.1)
@named load = PQLoad(; Pset=-.5, Qset=0)
bus1mtk = MTKBus(swing, load; name=:swingbus)
vertex1f = Bus(bus1mtk) # extract component function
```
Define a second bus as a slack
```@example concepts
bus2mtk = SlackDifferential(; name=:slackbus)
vertex2f = Bus(bus2mtk) # extract component function
```
Define the powerline connecting both nodes
```@example concepts
@named branch1 = PiLine()
@named branch2 = PiLine()
linemtk = MTKLine(branch1, branch2; name=:powerline)
edgef = Line(linemtk) # extract component function
```
Define the graph, the network and extract initial conditions
```@example concepts
g = complete_graph(2)
nw = Network(g, [vertex1f, vertex2f], edgef)
u0 = NWState(nw) # extract parameters and state from models
u0.v[1, :swing₊θ] = 0 # set missing initial conditions
u0.v[1, :swing₊ω] = 1
```
Then we can solve the problem
```@example concepts
prob = ODEProblem(nw, uflat(u0), (0,1), pflat(u0))
sol = solve(prob, Rodas5P())
@assert OrdinaryDiffEqRosenbrock.SciMLBase.successful_retcode(sol) # hide
nothing #hide
```
And finally we can plot the solution:
```@example concepts
fig = Figure();
ax = Axis(fig[1,1])
lines!(ax, sol; idxs=VIndex(1,:busbar₊P), label="Power injection Bus", color=Cycled(1))
lines!(ax, sol; idxs=VIndex(1,:swing₊Pel), label="Power injection Swing", color=Cycled(2))
lines!(ax, sol; idxs=VIndex(1,:load₊P), label="Power injection load", color=Cycled(3))
axislegend(ax)

ax = Axis(fig[2,1])
lines!(ax, sol; idxs=VIndex(1,:busbar₊u_arg), label="swing bus voltage angle", color=Cycled(1))
lines!(ax, sol; idxs=VIndex(2,:busbar₊u_arg), label="slack bus voltage angle", color=Cycled(2))
axislegend(ax)
fig #hide
```

## Internals

Internally, we use different input/output conventions for bus and line models.
The predefined models `BusBar()` and `LineEnd()` are defined in the following way:

### Model: `BusBar()`
A busbar is a concrete model used in bus modeling.
It represents the physical connection within a bus, the thing where all injectors and lines attach.
```
           ┌──────────┐
i_lines ──→│          │  (t)
           │  Busbar  ├───o
  u_bus ←──│          │
           └──────────┘
```
It receives the sum of all line currents as an input and equals that to the currents flowing into the terminal.
As an output, it gives forwards the terminal voltage to the backend.

### Model: `LineEnd()`
A `LineEnd` model is very similar to the `BusBar` model.
It represents one end of a transmission line.
```
          ┌───────────┐
 u_bus ──→│           │  (t)
          │  LineEnd  ├───o
i_line ←──│           │
          └───────────┘
```

It has special input/output connectors which handle the network interconnection.
The main difference beeing the different input/output conventions for the network interface.
