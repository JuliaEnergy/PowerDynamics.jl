# Modeling Concepts

In general, PowerDynamics models power grids as a set of dynamical systems for both **nodes** and **edges** on a graph. Check out the [Mathematical Model](@extref NetworkDynamics) documentation of NetworkDynamics for the underlying concepts.

The simulation happens entirely in a synchronous dq-frame. Due to their conceptual similarity to
complex phasors, variables in this global dq-frame are referenced by subscripts `_r` and `_i` (for real and imaginary).
This helps distinguish the variables from local dq frames, e.g. a generator model might transform `u_r` and `u_i` into `u_d` and `u_q`.

Both edge and node models are so-called input-output-systems: the edges receive the voltage of adjacent nodes as an input, the nodes receive the currents on adjacent edges as an input.
In general, this leads to the following structure of a bus/node model:

```math
\begin{aligned}
M_{\mathrm v}\,\frac{\mathrm{d}}{\mathrm{d}t}x_{\mathrm v} &= f_{\mathrm v}\left(x_{\mathrm v}, \sum_k\begin{bmatrix}i^k_r\\ i^k_i\end{bmatrix}, p_{\mathrm v}, t\right)\\
\begin{bmatrix}u_r\\ u_i\end{bmatrix} &= g_{\mathrm v}(x_{\mathrm v},p_{\mathrm v}, t)
\end{aligned}
```
where $M_{\mathrm v}$ is the (possibly singular) mass-matrix, $x_{\mathrm v}$ are the internal states and $p_{\mathrm v}$ are parameters.
Function $f_{\mathrm v}$ describes the time evolution of the internal states while output equation $g_{\mathrm v}$ defines the output voltage. The input for the system is the sum of all inflowing currents from adjacent lines $k$.
Note how vertices are modeled as one-port systems; i.e., they receive the accumulated current from all connected lines, they can't distinguish which line provides which current.
For special cases, this limitation might be mitigated using [External Inputs](@extref NetworkDynamics).

!!! note "Nodal dynamics include injectors"
    An important distinction between our modeling and the modeling in many other libraries is that we include the **injector dynamics** inside the **node dynamics**. I.e. if you have a bus with a load and a generator, the overall node dynamics will include both machine and load dynamics within their equations. The modularity and model reuse on the bus level is provided by [ModelingToolkit.jl](https://github.com/SciML/ModelingToolkit.jl) integration.

The edge model on the other hand looks like this:
```math
\begin{aligned}
M_{\mathrm e}\,\frac{\mathrm{d}}{\mathrm{d}t}x_{\mathrm e} &= f_{\mathrm e}\left(x_{\mathrm e}, \begin{bmatrix} u_r^\mathrm{src}\\u_i^\mathrm{src}\end{bmatrix}, \begin{bmatrix} u_r^\mathrm{dst}\\u_i^\mathrm{dst}\end{bmatrix},p_\mathrm{e}, t\right)\\
\begin{bmatrix}i_r^\mathrm{src}\\i_i^\mathrm{src}\end{bmatrix} &= g^\mathrm{src}_{\mathrm e}\left(x_{\mathrm e}, \begin{bmatrix} u_r^\mathrm{src}\\u_i^\mathrm{src}\end{bmatrix}, \begin{bmatrix} u_r^\mathrm{dst}\\u_i^\mathrm{dst}\end{bmatrix}, p_\mathrm{e}, t\right)\\
\begin{bmatrix}i_r^\mathrm{dst}\\i_i^\mathrm{dst}\end{bmatrix} &= g^\mathrm{dst}_{\mathrm e}\left(x_{\mathrm e}, \begin{bmatrix} u_r^\mathrm{src}\\u_i^\mathrm{src}\end{bmatrix}, \begin{bmatrix} u_r^\mathrm{dst}\\u_i^\mathrm{dst}\end{bmatrix}, p_\mathrm{e}, t\right)\\
\end{aligned}
```
There are a few notable differences compared to bus models:
Edges are two-port systems, they have **two distinct** inputs and **two distinct** outputs. Namely, they receive the dq voltage from both source and destination end, and define the current for both ends separately. In very simple systems without losses, those output currents might be just antisymmetric, in general cases however the current on both ends can differ drastically.

!!! note
    Source and destination end of a line are purely conventional. It has nothing to do with the actual flow direction. Per convention from [Graphs.jl](https://github.com/JuliaGraphs/Graphs.jl), edges in undirected graphs always go from vertex with lower index to vertex with higher index, i.e. $15 \to 23$ never $23 \to 15$.


The above descriptions are important to understand what's happening inside the package.
However, since we use ModelingToolkit to define the individual models a lot of this complexity is hidden from the user.
In the following, we'll go through the most important concepts when designing models using ModelingToolkit.

## Relationship between ModelingToolkit and NetworkDynamics
A crucial part of using this library is understanding the relationship between ModelingToolkit models and NetworkDynamics.

In a nutshell, ModelingToolkit models are **symbolic models**; i.e., they consist of symbolic equations which are not yet "compiled" for use as a numeric model.
The modeling in MTK is very flexible and similar to the Modelica language.
What we need in the end is models in the structure defined in the equations above.
For that, we need the MTK models to have a specific structure. Then we can use the [`compile_bus`](@ref) and [`compile_line`](@ref) function to compile the MTK models and create [`EdgeModel`](@extref NetworkDynamics Component-Models-with-MTK) and [`VertexModel`](@extref NetworkDynamics Component-Models-with-MTK) objects from them.
Those objects are not symbolic anymore but compiled numeric versions of the symbolically created systems.

## ModelingToolkit Models
### Terminal

The `Terminal`-Connector is an important building block for every model.
It represents a connection point with constant voltage in dq-coordinates `u_r` and `u_i` and enforces the Kirchhoff constraints `sum(i_r)=0` and `sum(i_i)=0`.

### Modeling of Buses
#### [Model class `Injector`](@id Injector Interface)

An injector is a class of components with a single `Terminal()` (called `:terminal`).
Examples for injectors might be Generators, Shunts, Loads.
```asciiart
(t)   ┌──────────┐
 o←───┤ Injector │
      └──────────┘
```
The current for injectors is always in injector convention; i.e., positive currents flow *out* of the injector *towards* the terminal.

!!! note "Model classes"
    Model "classes" are nothing formalized. In this document, a model class is just a description for some `System` from `ModelingToolkit.jl`, which satisfies certain requirements.
    For example, any `System` is considered an "Injector" if it contains a connector `Terminal()` called `:terminal`.

You can check if a model satisfies the injector interface using the [`isinjectormodel`](@ref) function.

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
            # if possible, it's better for the solver to explicitly provide algebraic equations for the current
            terminal.i_r ~ (Pset*terminal.u_r + Qset*terminal.u_i)/(terminal.u_r^2 + terminal.u_i^2)
            terminal.i_i ~ (Pset*terminal.u_i - Qset*terminal.u_r)/(terminal.u_r^2 + terminal.u_i^2)
        end
    end
    MyPQLoad(name=:pqload) #hide
    nothing #hide
    ```

#### [Model class `MTKBus`](@id MTKBus Interface)
A `MTKBus` is a class of models, which are used to describe the dynamic behavior of a full bus in a power grid.
Each `MTKBus` must contain a predefined model of type `BusBar()` (named `:busbar`).
This busbar represents the connection point to the grid.
Optionally, it may contain various injectors.
If there are no injectors, the model just describes a junction bus; i.e., a bus that just satisfies the Kirchhoff constraint for the flows of connected lines.

```asciiart
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

For simple models (direct connections of a few injectors), it is possible to use the convenience method `MTKBus(injectors...)` to create the composite model based on provided injector models.

You can check if a model satisfies the bus interface using the [`isbusmodel`](@ref) function.

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
    Alternatively, for that system you could have just called
    ```@example concepts
    mybus = MTKBus(Swing(;name=:swing), PQLoad(;name=:load))
    nothing #hide
    ```
    to get an instance of a model which is structurally equivalent to `MyMTKBus`.

### Line Modeling
#### [Model class `Branch`](@id Branch Interface)
A branch is the two-port equivalent to an injector.
It needs to have two `Terminal()`s, one is called `:src`, the other `:dst`.

Examples for branches are: PI-model branches, dynamic RL branches or transformers.
```asciiart
(src) ┌──────────┐ (dst)
  o←──┤  Branch  ├──→o
      └──────────┘
```
*Both* ends follow the injector interface; i.e., current leaving the device towards the terminals is always positive.

You can check if a model satisfies the branch interface using the [`isbranchmodel`](@ref) function.

!!! details "Code example: algebraic R-line"
    ```@example concepts
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

#### [Model class: `MTKLine`](@id MTKLine Interface)
Similar to the `MTKBus`, a `MTKLine` is a model class which represents a transmission line in the network.

It must contain two `LineEnd()` instances, one called `:src`, one called `:dst`.

```asciiart
 ┌────────────────────────────────────────────────┐
 │ MTKLine          ┌──────────┐                  │
 │  ┌─────────┐  ┌──┤ Branch A ├──┐  ┌─────────┐  │
 │  │ LineEnd │  │  └──────────┘  │  │ LineEnd │  │
 │  │  :src   ├──o                o──┤  :dst   │  │
 │  │         │  │  ┌──────────┐  │  │         │  │
 │  └─────────┘  └──┤ Branch B ├──┘  └─────────┘  │
 │                  └──────────┘                  │
 └────────────────────────────────────────────────┘
```

Simple line models, which consist only of valid `Branch` models, can be instantiated using the `MTKLine(branches...)` constructor.

More complex models can be created manually.
For example, you could define a dynamic multi-branch DC line model by chaining
(possibly very complex and nested) source and destination inverter/rectifier models
with one or many dc branches.
```asciiart
┌───────────────────────────────────────────────────────────┐
│ MTKLine                ┌──────────┐                       │
│                       ┌┤ DC Br. A ├┐                      │
│┌─────────┐ ┌─────────┐│└──────────┘│┌─────────┐┌─────────┐│
││ LineEnd ├─┤ src inv ├o            o┤ dst inv ├┤ LineEnd ││
│└─────────┘ └─────────┘│┌──────────┐│└─────────┘└─────────┘│
│   :src                └┤ DC Br. B ├┘              :dst    │
│                        └──────────┘                       │
└───────────────────────────────────────────────────────────┘
```

You can check if a model satisfies the line interface using the [`islinemodel`](@ref) function.

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
Both `MTKLine` and `MTKBus` are still purely symbolic ModelingToolkit models.
However, they have an important property: they possess the correct input-output structure
and variable names to be compiled into [`VertexModel`](@extref NetworkDynamics Component-Models-with-MTK) and [`EdgeModel`](@extref NetworkDynamics Component-Models-with-MTK)
models.
To do so, PowerDynamics.jl provides the [`compile_line`](@ref) and [`compile_bus`](@ref) functions.

At their core, both `compile_*` functions use ModelingToolkit's [`mtkcompile`](@extref ModelingToolkit.mtkcompile) to perform **symbolic simplifications** on your models and reduce the number of states.
Most notably, this process can drastically reduce the number of equations, while all previously defined states remain "observable", i.e. inspectable after simulation.
For example, in the above code example of the PQ load we defined equations for active and reactive powers $P$ and $Q$. Those equations don't add anything to the actual behavior of the system,
however they will be kept around as so-called "observed" states, meaning we can reconstruct and plot them from dynamical simulations.

!!! tip "Tip: Introduce states of interest as observables"
    It is often useful to add derived quantities of interest explicitly to your models. For example, if you're interested in some internal
    voltage angle just define an equation `u_angle ~ atan(u_d, u_q)`. If nothing else depends on it, this equation will be symbolically reduced, 
    i.e. you don't add any overhead to your simulation but it will be accessible after the simulation. 
    This is often far more convenient than "reconstructing" such states manually!


### (MTK) Bus Model to VertexModel: `compile_bus`
```asciiart
                                      ╔═════════════════════════╗
                                      ║ VertexModel (compiled)  ║
┌────────────────────┐      Network   ║  ┌────────────────────┐ ║
│MTKBus   ┌─────────┐│     interface  ║  │MTKBus   ┌─────────┐│ ║
│        ┌┤Generator││                ║  │        ┌┤Generator││ ║
│┌──────┐│└─────────┘│      current ────→│┌──────┐│└─────────┘│ ║
││BusBar├o           │  =>            ║  ││BusBar├o           │ ║
│└──────┘│┌────┐     │      voltage ←────│└──────┘│┌────┐     │ ║
│        └┤Load│     │                ║  │        └┤Load│     │ ║
│         └────┘     │                ║  │         └────┘     │ ║
└────────────────────┘                ║  └────────────────────┘ ║
                                      ╚═════════════════════════╝
```

### (MTK) Line Model to EdgeModel: `compile_line`
```asciiart

                                        ╔══════════════════════════════╗
                                        ║ EdgeModel (compiled)         ║
┌─────────────────────────────┐     src ║ ┌──────────────────────────┐ ║ dst
│MTKLine   ┌───────┐          │  vertex ║ │MTKLine   ┌────┐          │ ║ vertex
│         ┌┤BranchA├┐         │         ║ │         ┌┤    ├┐         │ ║
│┌───────┐│└───────┘│┌───────┐│     u ───→│┌───────┐│└────┘│┌───────┐│←─── u
││LineEnd├o         o┤LineEnd││  =>     ║ ││LineEnd├o      o┤LineEnd││ ║
│└───────┘│┌───────┐│└───────┘│     i ←───│└───────┘│┌────┐│└───────┘│───→ i
│  :src   └┤BranchB├┘  :dst   │         ║ │         └┤    ├┘         │ ║
│          └───────┘          │         ║ │          └────┘          │ ║
└─────────────────────────────┘         ║ └──────────────────────────┘ ║
                                        ╚══════════════════════════════╝
```

### End to End Example
Putting the knowledge from this document together, we can start a short simulation of an example network:
```@example concepts
using PowerDynamics, PowerDynamics.Library, ModelingToolkit
using Graphs, NetworkDynamics
using OrdinaryDiffEqRosenbrock, OrdinaryDiffEqNonlinearSolve
using CairoMakie
nothing #hide
```

First, we define an `MTKBus` consisting of two predefined injector models from the
library: a `Swing` generator model and a `PQLoad`.
To do so, we use the [`MTKBus(injectors...)`](@ref MTKBus) constructor.
```@example concepts
@named swing = Swing(; Pm=1, V=1, D=0.1)
@named load = PQLoad(; Pset=-.5, Qset=0)
bus1mtk = MTKBus(swing, load; name=:swingbus)
show(stdout, MIME"text/plain"(), bus1mtk) #hide
```
This results in an MTK model, which fulfills the `MTKBus` interface and thus can be compiled into an actual `VertexModel` for simulation:

```@example concepts
vertex1f = compile_bus(bus1mtk) # extract component function
```

As a second bus in this example, we use a [`SlackDifferential`]() from the library.
This model is not an Injector but an MTKBus directly, as it does not make sense to connect anything else to a slack bus.
```@example concepts
bus2mtk = SlackDifferential(; name=:slackbus)
vertex2f = compile_bus(bus2mtk) # extract component function
```

For the connecting line, we instantiate two [`PiLine`]() from the library.
Each PiLine fulfills the Branch interface. Therefore we can define a `MTKLine` model
by putting both Branches in parallel:
```@example concepts
@named branch1 = PiLine()
@named branch2 = PiLine()
linemtk = MTKLine(branch1, branch2; name=:powerline)
show(stdout, MIME"text/plain"(), bus1mtk) #hide
```
Similar to before, we need to compile the MTKModel by calling [`compile_line`](@ref).
```@example concepts
edgef = compile_line(linemtk) # extract component function
```

To simulate the system, we place both components on a graph and define their network topology.
We define both graph topology as well as the models for the individual components.
```@example concepts
g = complete_graph(2)
nw = Network(g, [vertex1f, vertex2f], edgef)
u0 = NWState(nw) # extract parameters and state from models
u0.v[1, :swing₊θ] = 0 # set missing initial conditions
u0.v[1, :swing₊ω] = 1
nothing #hide
```
Then we can solve the problem
```@example concepts
prob = ODEProblem(nw, u0, (0,1))
sol = solve(prob, Rodas5P())
@assert OrdinaryDiffEqRosenbrock.SciMLBase.successful_retcode(sol) #hide
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
It represents the physical connection within a bus, the component where all injectors and lines attach.
```asciiart
           ┌──────────┐
i_lines ──→│          │  (t)
           │  Busbar  ├───o
  u_bus ←──│          │
           └──────────┘
```
It receives the sum of all line currents as an input and balances this with the currents flowing into the terminal.
As an output, it forwards the terminal voltage to the backend.

### Model: `LineEnd()`
A `LineEnd` model is very similar to the `BusBar` model.
It represents one end of a transmission line.
```asciiart
          ┌───────────┐
 u_bus ──→│           │  (t)
          │  LineEnd  ├───o
i_line ←──│           │
          └───────────┘
```

It has special input/output connectors which handle the network interconnection.
The main difference being the distinct input/output conventions for the network interface.
