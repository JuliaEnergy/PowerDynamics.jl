#=
# Getting Started

The goal of this tutorial is to get you started with PowerDynamics.jl. Its main goal is to show you the different
"stages" of a typical simulation workflow and introducing some jargon along the way.

The system to model is a simple 3 bus system:
- Bus 1: ideal droop inverter
- Bus 2: a constant Y load
- Bus 3: a second constant Y load

```
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
Models may either come from a library (like the `IdealDroopInverter` and `ConstantYLoad` below
or can be defined by the user.

### Inverter Model
In this case, we want the first bus to be an ideal droop inverter.

Often, model definition will be a multi step process:
**First** we define an ["injector model"](@ref Injector Interface), in this case our
inverter:
=#
inverter_model = Library.IdealDroopInverter(; name=:droop, Vset=1)
show(stdout, MIME"text/plain"(), inverter_model) # hide
#=
This model is an equation-based/symbolic model representing the dynamics.
It is based on the great [ModelingToolkit.jl Library](https://mtk.sciml.ai/stable/).

**Second**: we build a ["Bus Model"](@ref MTKBus Interface), which connects
the injector to a Busbar.
This model still lives in the equation-based/symbolic domain.
=#
bus_model = MTKBus(inverter_model; name=:invbus)
show(stdout, MIME"text/plain"(), bus_model) # hide
#=
**Third** we compile the symbolic model into a julia function for numeric simulation.
Doing so, we get a [`VertexModel`](@extref NetworkDynamics.VertexModel-Tuple{}), which is an
object from our backend [NetworkDynamics.jl](https://juliadynamics.github.io/NetworkDynamics.jl/stable/)
=#
bus1 = compile_bus(MTKBus(inverter_model); vidx=1)
#=
Notably, this model **not** a symbolic model anymore. The equations have been reduced
and transformed into an nonlinear discriptor model.
For more information on the different modle types see the [Modeling Concepts](@ref) docs.
You can check out the NetworkDynamics.jl doc on the underlying [mathematical model](@extref Mathematical-Model).

In the printout above you can see that we consider different types of varaibles in our models:
- the **input** is allways the current comming flowing from the attatched powerlines into the bus,
- the **output** is always the voltage at the busbar,
- the **states** are dynamical or algebraic states in the sense of a Differential-Algebraic-Equation (DAE) model,
- and **parameters** are static values that stay mostly constant during simulation and define the system behavior.

There is a 5th class of states not shown above: **observables**. Observables are time values, which are not
states in the sense of a DAE but can be reconstructed from the states, inputs, outputs and parameters.
Thus, they don't need to be "solved" for numerically, but they can be reconstructed in post-processing.

### Load Models
For the two loads we use the predfined `ConstantYLoad` model from the Library, and compile them:
=#
load_model = Library.ConstantYLoad(; name=:load)
bus2 = compile_bus(MTKBus(load_model); name=:loadbus, vidx=2)
bus3 = compile_bus(MTKBus(load_model); name=:loadbus, vidx=3)
#=
### Powerline Models
Lastly, we need to define three powerlines.
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
In power systems specificially, it is common to solve a simpler system first -- the so called
"powerflow" problem.

In the powerflor problem, we neglect all the node dynamics and consider only 4 variables at
each bus:
- the active power $P$,
- the reactive power $Q$,
- the voltage magnitude $V$ and
- the voltage angle $\theta$.

In the simple-most powerflow, each bus can then be classified into one of three types:
- **Slack Bus**: The voltage magnitude and angle are fixed (typically used for one bus in the system), $P$ and $Q$ is considered free.
- **PV Bus** $P$ and $V$ are fixed, $Q$ and $\theta$ are free. Often used for generator buses or any buses with active voltage control.
- **PQ Bus** $P$ and $Q$ are fixed, $V$ and $\theta$ are free. Typically used for load buses.

So each component essentially introduces two algebraic equations and two free variables -- the system is
then solved for the free variables such that all equations are satisfied.

### Attaching Powerflow Models
In PowerDynamics.jl we can attach the powerflow models to the dynamic bus models using the [`set_pfmodel!`](@ref) function.
=#
set_pfmodel!(bus1, pfSlack(V=1))
set_pfmodel!(bus2, pfPQ(P=-0.4, Q=-0.3))
set_pfmodel!(bus3, pfPQ(P=-0.6, Q=-0.2))
nothing # hide

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
The Network object tells us, that we've just define a system with 7 States and 39 parameters. We have 3 vertices of 2 unique types
(the inverter bus and the load bus) and 3 edges of a single unique type (all powerlines are the same piline type).

The "states" and "parameters" allready hint at a very important property of PowerDynamics/NetworkDynamics:
in the end, the whole network is just a big DAE system of the form
```math
\mathbf{M}\,\dot{\mathbf{x}} = f(\mathbf{x}, \mathbf{p}, t)
```
where $\mathbf{x}$ are the states and $\mathbf{p}$ the parameters.
This is very important to keep in mind, because it allows us to integreate seamlessly with the whole SciMLEcosystem and
most importantly [DifferentailEquations.jl](https://diffeq.sciml.ai/stable/).

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
- Bus 2&3: 2 states, 2 parameters each
- Lines 1,2 and 3: 0 states, 9 parameters each
In sum, we get the 7 states and 39 parameters.

!!! tip "Advanced: State and Parameter Ordering"
    Even though the states and parameters are essential "just stacked" on top of eachother, the ordering is not
    trivial due to performance reasons. Never rely on the ordering of states or parameters in the full system!
    PowerDynamics and NetworkDynamics provides lots of helper functions for so-called "SymbolicIndexing" to circumvent this.

### Initializing the System via Powerflow
With the network constructed we can finally find our equilibrium point.
We do so using [`initialize_from_pf`](@ref):
=#
s0 = initialize_from_pf(nw; verbose=true, subverbose=true)
nothing #hide
#=
This function actually does quite a lot.
1. It extracts the powerflow model for each component constructing the powerflow problem.
2. It solves the powerflow problem.
3. From the solution, we know the voltages and powers at each bus. Consequently, we also know the currents at each bus.
   This means we can directly "map" the powerflow solution to the inputs and outputs of each dynamic model.
4. For each dynamic model, we "fix" the inputs and outputs and try to find values for the states (and potentially parameters) such that
   the component model is in equilibrium. This is done using a nonlinear solver.
=#

using Test #hide #src
@test_logs (:info, "Initialization problem is fully constrained. Created NonlinearLeastSquaresProblem for \
                    [:droop₊δ, :droop₊Qfilt, :droop₊Pfilt, :droop₊Pset, :droop₊Qset]") match_mode=:any initialize_from_pf(nw; verbose=false, subverbose=true); #hide #src
nothing #hide #src
#=
In the log statments we see, which variables/parameters where considered free during the initialization of each component.
This behavior can be finetuned in a lot of ways, which are beyond the scope of this tutorial.
However, here we see that, for example, the complex parameter $Y = B + j\,G$ ofthe constant Y-load was initialy left free but
was initialized from the powerflow solution. This means, that your $Y$ is no set in a way, that it draws the correct
amount of power at the given voltage.

Similarily, we see that the inverter bus had the parameters $P_{set}$ and $Q_{set}$ free, which were also initialized from the powerflow solution.
This is important, becaus we need to achive powerbalance in the system, and due to the losses in the lines its not possible to
known the exact power injections a priori.

The return value of the `initialize_from_pf` function is a so-called [`NWState`](@extref NetworkDynamics.NWState) object,
which wraps flat $x$ and $p$ vectors and provides a lot of helper functions to access and modify states (including observables) and parameters.

Lets inspect this object further:
=#
s0
#=
On the highes level, we see the values of the 7 states in the network and their symbolic indices.
Those indices, can be used to acces values direclty:
=#
s0[VIndex(1, :droop₊P_filt)]
#=
!!! tip
    In most julia dev environments you can use "\_+<TAB>" to autocomplete the MTK namespace separator "₊".

Often, you want to access observables or parameters instead of states. There is a whole
filtering and accessing mechanis you can use for that. For example, in PD.jl each bus has the states `:busbar₊P` and `:busbar₊Q`.
We can inspect them on all vertices using:
=#
s0.v(:, ["busbar₊P", "busbar₊Q"])
#=
In the output we clearly see, how the load busses draw exactly the amout of power we sepcified in the powerflow models.
On the inverter bus however, we inject slighly more power then then loads demand to compensate for the line losses.

Similarily, we can access all node parameters at initial state using
=#
s0.v.p
#=
Here we see, how $P_{set}$ and $Q_{set}$ of the inverter where initialized in a way, that they match the powerflow solution.

There is a lot more functionality in the `NWState` objects, see the [Symbolic
Indexing docs of ND.jl](@extref Symbolic-Indexing) and especially the
[`FilteringProxy`](@extref) for more details.

## Time Domain Simulation


=#

prob = ODEProblem(nw, uflat(s0), (0.0, 1.0), copy(pflat(s0)))
sol = solve(prob, Rodas5P())

affect = ComponentAffect([], [:load₊B]) do u, p, ctx
    @info "Increasing load B by 10%"
    p[:load₊B] = p[:load₊B] * 1.1
end
cb = PresetTimeComponentCallback(1.0, affect)
set_callback!(bus2, cb)

prob = ODEProblem(nw, uflat(s0), (0.0, 20.0), copy(pflat(s0)), callback=get_callbacks(nw))
sol = solve(prob, Rodas5P())

let
    fig = Figure()
    ts = range(sol.t[begin], sol.t[end]; length=1000)
    ax = Axis(fig[1, 1]; xlabel="Time [s]", ylabel="Voltage [pu]")
    for i in 1:3
        lines!(ax, sol, idxs=VIndex(i, :busbar₊u_mag), color=Cycled(i))
    end
    ax = Axis(fig[2,1]; xlabel="Time [s]", ylabel="Active Power Injection [rad]")
    lines!(ax, sol, idxs=VIndex(1, :busbar₊P))
    ax = Axis(fig[3,1]; xlabel="Time [s]", ylabel="Active Power Load [rad]")
    for i in 2:3
        lines!(ax, sol, idxs=VIndex(i, :busbar₊P), color=Cycled(i))
    end
    fig
end


nw[VIndex(2)]

bus2

uflat(s0)



dump_initial_state(bus1; obs=false)
dump_initial_state(bus2; obs=true)
dump_initial_state(bus3; obs=true)
