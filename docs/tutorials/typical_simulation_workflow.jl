#=
# [Typical Simulation Workflow](@id typical-simulation-workflow)

The goal of this tutorial is to get you started with PowerDynamics.jl. We'll walk you through the different
"stages" of a typical simulation workflow while introducing key terminology along the way.

This tutorial can be downloaded as a normal Julia script [here](@__NAME__.jl). #md

The system to model is a simple 3 bus system:
- Bus 1: ideal droop inverter
- Bus 2: a constant Y load
- Bus 3: a second constant Y load

All buses are connected with standard pi-model power lines.
```asciiart
    ╭───────╮
2 ┯━┿       ┿━┯ 3
  ↓ │   ╭───╯ ↓
    ┷━┯━┷ 1
      │
     (~)
```
=#

using PowerDynamics
using PowerDynamics: Library
using OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqNonlinearSolve
using CairoMakie

#=
## Stage I: Defining Dynamical Models

The first phase of the workflow is typically to define your dynamical models.
Models may either come from a library (like the `IdealDroopInverter` and `ConstantYLoad` below)
or they can be defined by the user.

### Inverter Model
In this case, we want the first bus to be an ideal droop inverter.

Often, model definition will be a multi step process:

**First**: we define an ["injector model"](@ref Injector-Interface), in this case our
inverter:
```asciiart
(t) ┌────────────────┐
 o──┤ Droop Inverter │
    └────────────────┘
```
=#
inverter_model = Library.IdealDroopInverter(; name=:droop, Vset=1)
nothing #hide #md
#=
This model is an equation-based/symbolic model representing the dynamics.
It is based on the great [ModelingToolkit.jl Library](https://mtk.sciml.ai/stable/).

**Second**: we build a ["Bus Model"](@ref MTKBus-Interface), which connects
the injector to a Busbar.
This model still lives in the equation-based/symbolic domain.

```asciiart
┌───────────────────────────────┐
│BusModel                       │
│┌────────┐   ┌────────────────┐│
││ BusBar ├─o─┤ Droop inverter ││
│└────────┘   └────────────────┘│
└───────────────────────────────┘
```
=#
bus_model = MTKBus(inverter_model; name=:invbus)
nothing #hide #md
#=
**Third**: we compile the symbolic model into a julia function for numeric simulation.
Doing so, we get a [`VertexModel`](@extref NetworkDynamics.VertexModel-Tuple{}), which is an
object from our backend [NetworkDynamics.jl](https://juliadynamics.github.io/NetworkDynamics.jl/stable/)
```asciiart
           ╔════════════════════════════════╗
 Network   ║ VertexModel (compiled)         ║
interface  ║  ┌───────────────────────────┐ ║
           ║  │BusModel                   │ ║
 current ────→│┌──────┐ ┌────────────────┐│ ║
           ║  ││BusBar├o┤ Droop inverter ││ ║
 voltage ←────│└──────┘ └────────────────┘│ ║
           ║  └───────────────────────────┘ ║
           ╚════════════════════════════════╝
```
=#
bus1 = compile_bus(MTKBus(inverter_model); vidx=1)
#=
Note that this model is **no longer symbolic**. The equations have been reduced
and transformed into a nonlinear descriptor model.
For more information on the different model types, see the [Modeling Concepts](@ref) docs.
You can check out the NetworkDynamics.jl doc on the underlying [mathematical model](@extref Mathematical-Model).

In the printout above, you can see that we consider different types of variables in our models:
- the **input** is always the current flowing from the attached power lines into the bus,
- the **output** is always the voltage at the busbar,
- the **states** are dynamical or algebraic states in the sense of a Differential-Algebraic-Equation (DAE) model,
- and **parameters** are static values that stay mostly constant during simulation and define the system behavior.

There is a 5th class of states not shown above: **observables**. Observables are time dependent values, which are not
states in the sense of a DAE but can be reconstructed from the states, inputs, outputs and parameters.
Thus, they don't need to be "solved" for numerically, but they can be reconstructed in post-processing.
A simple example of an "observed" state would be a voltage angle or the active and reactive power at some bus.

### Load Models
For the two loads, we use the predefined `ConstantYLoad` model from the Library and compile them:
=#
load_model = Library.ConstantYLoad(; name=:load)
bus2 = compile_bus(MTKBus(load_model); name=:loadbus, vidx=2)
bus3 = compile_bus(MTKBus(load_model); name=:loadbus, vidx=3)
#=
### Power Line Models
Lastly, we need to define three power lines.
The workflow is similar to the bus models:
=#
l = MTKLine(Library.PiLine(; name=:piline))
line12 = compile_line(l; src=1, dst=2, piline₊R=0.01)
line13 = compile_line(l; src=1, dst=3, piline₊R=0.01)
line23 = compile_line(l; src=2, dst=3, piline₊R=0.01)
#=
Now we're all set for the next stage.

## Stage II: Initialization

When simulating power systems (or any large dynamical system for that matter),
it is quite typical to start from a steady state/equilibrium point.
In general, it is not trivial to find such a point for a large nonlinear system.
In power systems specifically, it is common to solve a simpler system first -- the so-called
"power flow" problem.

In the power flow problem, we neglect all the node dynamics and consider only 4 variables at
each bus:
- the active power $P$,
- the reactive power $Q$,
- the voltage magnitude $V$ and
- the voltage angle $\theta$.

In the simplest power flow, each bus can then be classified into one of three types:
- **Slack Bus**: The voltage magnitude and angle are fixed (typically used for one bus in the system), $P$ and $Q$ is considered free.
- **PV Bus** $P$ and $V$ are fixed, $Q$ and $\theta$ are free. Often used for generator buses or any buses with active voltage control.
- **PQ Bus** $P$ and $Q$ are fixed, $V$ and $\theta$ are free. Typically used for load buses.

So each component essentially introduces two algebraic equations and two free variables -- the system is
then solved for the free variables such that all equations are satisfied.

### Attaching Power Flow Models
In PowerDynamics.jl, we can attach the power flow models to the dynamic bus models using the [`set_pfmodel!`](@ref) function.
=#
set_pfmodel!(bus1, pfSlack(V=1))
set_pfmodel!(bus2, pfPQ(P=-0.4, Q=-0.3))
set_pfmodel!(bus3, pfPQ(P=-0.6, Q=-0.2))
nothing #hide #md

#=
### Building the Network
Now we can build the network using the [`Network`](@extref NetworkDynamics.Network) constructor from NetworkDynamics.jl
This constructor takes a list of VertexModels (the buses) and a list of EdgeModels (the powerlines) and connects them to a
network.
In general, we also need to define the topology of the undelying graph. In this case however, this is not necessary because
we told each component at the compile step where it is placed in the network (see the `vidx`, `src` and `dst` arguments above)
=#
nw = Network([bus1, bus2, bus3], [line12, line13, line23])
#=
The Network object tells us that we've just defined a system with 7 States and 39 parameters. We have 3 vertices of 2 unique types
(the inverter bus and the load bus) and 3 edges of a single unique type (all power lines are the same pi-line type).

The "states" and "parameters" already hint at a very important property of PowerDynamics/NetworkDynamics:
in the end, the whole network is just a big DAE system of the form
```math
\mathbf{M}\,\dot{\mathbf{x}} = f(\mathbf{x}, \mathbf{p}, t)
```
where $\mathbf{x}$ are the states and $\mathbf{p}$ the parameters.
This is very important to keep in mind, because it allows us to integrate seamlessly with the whole SciML ecosystem and,
most importantly, [DifferentialEquations.jl](https://diffeq.sciml.ai/stable/).

=#
@assert dim(nw) == 7 #hide
@assert pdim(nw) == 39 #hide
@assert dim(bus1) == 3 #hide
@assert pdim(bus1) == 8 #hide
@assert dim(bus2) == dim(bus2) == 2 #hide
@assert pdim(bus2) == pdim(bus3) == 2 #hide
@assert dim(line12) == dim(line13) == dim(line23) == 0 #hide
@assert pdim(line12) == pdim(line13) == pdim(line23) == 9 #hide
#=
The 7 states are essentially just the states of our models stacked on top of eachother.
Look at the representation of our Vertex and EdgeModels above to see their contribution:
- Bus 1: 3 states, 8 parameters
- Bus 2 & 3: 2 states, 2 parameters each
- Lines 1,2 and 3: 0 states, 9 parameters each
In sum, we get the 7 states and 39 parameters.

!!! tip "Advanced: State and Parameter Ordering"
    Even though the states and parameters are essential "just stacked" on top of eachother, the ordering is not
    trivial due to performance reasons. Never rely on the ordering of states or parameters in the full system!
    PowerDynamics and NetworkDynamics provides lots of helper functions for so-called "SymbolicIndexing" to circumvent this.

### Initializing the System via Power Flow
With the network constructed, we can finally find our equilibrium point.
We do so using [`initialize_from_pf`](@ref):
=#
s0 = initialize_from_pf(nw; verbose=true, subverbose=true)
nothing #hide #md
#=
This function actually does quite a lot.
1. It extracts the powerflow model for each component constructing the powerflow problem.
2. It solves the powerflow problem.
3. From the solution, we know the voltages and powers at each bus. Consequently, we also know the currents at each bus.
   This means we can directly "map" the powerflow solution to the inputs and outputs of each dynamic model.
4. For each dynamic model, we "fix" the inputs and outputs and try to find values for the states (and potentially parameters) such that
   the component model is in equilibrium. This is done using a nonlinear solver.
=#

@assert Set(free_p(nw[VIndex(1)])) == Set([:droop₊Pset, :droop₊Qset]) # hide
@assert Set(free_u(nw[VIndex(1)])) == Set([:droop₊δ, :droop₊Qfilt, :droop₊Pfilt]) # hide
#=
In the log statements, we see which variables/parameters were considered free during the initialization of each component.
This behavior can be fine-tuned in a lot of ways, which are beyond the scope of this tutorial.
However, here we see that, for example, the complex parameter $Y = B + j\,G$ of the constant Y-load was initially left free but
then initialized from the powerflow solution. This means that $Y$ is now set in a way that it draws the correct
amount of power at the given voltage.

Similarly, we see that the inverter bus had the parameters $P_{set}$ and $Q_{set}$ free, which were also initialized from the powerflow solution.
This is important, because we need to achieve power balance in the system, and due to the losses in the lines it's not possible to
know the exact power injections a priori.

The return value of the `initialize_from_pf` function is a so-called [`NWState`](@extref NetworkDynamics.NWState) object,
which wraps flat $x$ and $p$ vectors and provides a lot of helper functions to access and modify states (including observables) and parameters.

Lets inspect this object further:
=#
s0
#=
At the highest level, we see the values of the 7 states in the network and their symbolic indices.
Those indices can be used to access values directly:
=#
s0[VIndex(1, :droop₊Pfilt)]
#=
!!! tip
    In most julia dev environments you can type `\_+<TAB>` to autocomplete the MTK namespace separator `₊`.

Often, you want to access observables or parameters instead of states. There is a whole
filtering and access mechanism you can use for that. For example, in PD.jl, each bus has the states `:busbar₊P` and `:busbar₊Q`.
We can inspect them on all vertices using:
=#
s0.v(:, [:busbar₊P, :busbar₊Q])

#=
In the output, we clearly see how the load buses draw exactly the amount of power we specified in the power flow models.
On the inverter bus, however, we inject slightly more power than the loads demand to compensate for the line losses.

Similarly, we can access all node parameters at initial state using
=#
s0.v.p
#=
Here we see how $P_{set}$ and $Q_{set}$ of the inverter were initialized in a way, that they match the powerflow solution.

There is a lot more functionality in the `NWState` objects, see the [Symbolic
Indexing docs of ND.jl](@extref Symbolic-Indexing) and especially the
[`FilteringProxy`](@extref NetworkDynamics.FilteringProxy) for more details.

## Stage III: Time Domain Simulation
With the initialized state, we can finally simulate the system in time domain.

### Perturbing the System
Since we start from an equilibrium point, we expect the system to stay there if we don't
perturb it.
Therefore, to get interesting results, we need to perturb the system.

The simplest way to perturb the system is to change a parameter. For example, let's
increase the admittance at bus 2 by 10% after 0.1 seconds.

For that, we define a so-called "callback function", more specifically a preset time callback,
which is triggered at a specific simulation time and modifies the parameters. The simulation then
continues.
General information on callbacks in Differential Equations can be found in the
[DiffEq.jl docs](@extref DiffEq callbacks). Specific extensions for NetworkDynamics.jl can be found in the
[NetworkDynamics.jl callback docs](@extref NetworkDynamics Callbacks).

We define the callback and attach it to bus 2 (our first load) like this:
=#
affect = ComponentAffect([], [:load₊G, :load₊B]) do u, p, ctx
    @info "Increase load admittance Y by 10% at t=$(ctx.t)"
    p[:load₊G] = p[:load₊G] * 1.1
    p[:load₊B] = p[:load₊B] * 1.1
end
cb = PresetTimeComponentCallback(0.1, affect)
set_callback!(bus2, cb)
bus2 #hide


#=
With the callback defined, we can finally create and solve the [`ODEProblem`](@extref SciMLBase.ODEProblem(::NetworkDynamics.Network, ::NetworkDynamics.NWState, ::Any)):
=#
prob = ODEProblem(nw, s0, (0.0, 5.0))
sol = solve(prob, Rodas5P());
nothing #hide #md

#=
## Stage IV: Postprocessing and Visualization

Once we have the solution object, we can use it like any other solution from DifferentialEquations.jl.
Most importantly, we can use symbolic indices to access states, parameters and observables.

### Plotting Results
For example, we can quickly plot the frequency response of the droop using
=#
lines(sol, idxs=VIndex(1, :droop₊ω); axis=(;xlabel="Time [s]", ylabel="Frequency ω [pu]"))
#=
We clearly see how the increased active power leads to a drop in frequency. This drop is
then compensated by the droop control (i.e., we stabilize at a lower frequency).

Of course, we can also create more complex plots, such as this one showing the active and reactive power at each bus:
=#
let
    fig = Figure(size=(1000,600))
    ax = Axis(fig[1, 1]; xlabel="Time [s]", ylabel="Active Power Load [pu]")
    for i in 2:3
        lines!(ax, sol, idxs=VIndex(i, :busbar₊P), color=Cycled(i))
    end
    axislegend(ax)
    ax = Axis(fig[1,2]; xlabel="Time [s]", ylabel="Reactive Power Load [pu]")
    for i in 2:3
        lines!(ax, sol, idxs=VIndex(i, :busbar₊Q), color=Cycled(i))
    end
    axislegend(ax)
    ax = Axis(fig[2,1]; xlabel="Time [s]", ylabel="Active Power Injection [pu]")
    lines!(ax, sol, idxs=VIndex(1, :busbar₊P))
    axislegend(ax)
    ax = Axis(fig[2,2]; xlabel="Time [s]", ylabel="Reactive Power Injection [pu]")
    lines!(ax, sol, idxs=VIndex(1, :busbar₊Q))
    axislegend(ax)
    fig
end
#=
Here we see that the active and reactive power demand shoots up in the beginning
after we increase Y.
The power demand then slowly decreases again.
This is probably due to a drop in voltage, which leads to lower power demand
on constant Y loads.

Let's plot the voltage to verify this:
=#
let
    fig = Figure()
    ax = Axis(fig[1, 1]; xlabel="Time [s]", ylabel="Voltage [pu]")
    for i in 1:3
        lines!(ax, sol, idxs=VIndex(i, :busbar₊u_mag), color=Cycled(i))
    end
    axislegend(ax; position=:rc)
    fig
end
#=

### Programmatic Access to Variables
Instead of plotting, we can also always use the solution interpolation to access values programmatically.
For example
=#
sol(1.0, idxs=VIndex(1,:droop₊ω))
#=
gives us the frequency of the droop inverter at time t=1.0s.

Similarly, we can extract time series by passing a vector of time points rather than a single point:
=#
sol([0.1,0.2,0.3], idxs=VIndex(1:3, :busbar₊u_mag))
#=
This code gives us the voltage magnitude at all three buses at the time points 0.1s, 0.2s, and 0.3s.

To deeply inspect a single point, we can also construct a `NWState` object from the solution for a specific time:
=#
s095 = NWState(sol, 0.95)
#=
gives us the state at t=0.95s. We can use the state object for inspection as we did before.
For example, we can inspect the power at the destination end of all lines at this point in time:
=#
s095.e(:, :dst₊P)
