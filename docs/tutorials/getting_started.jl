#=
# [Getting Started with PowerDynamics.jl](@id getting-started)

This tutorial introduces the core ideas behind PowerDynamics.jl and its relationship to the SciML ecosystem.

This tutorial can be downloaded as a normal Julia script [here](@__NAME__.jl). #md

The most important distinction in contrast to other tools is that PowerDynamics.jl is a *modeling framework* rather than a *simulation tool*.
At its core, a dynamic powergrid model is just a set of differential-algebraic equations (DAEs) that describe the evolution of the system over time.
PowerDynamics.jl helps you to build these DAE models in a modular way, and then simulate them using the powerful solvers from the SciML ecosystem.

```asciiart
         ╭───────────────────────────────────────────────╮
         │ PowerGrid Model                               │
         │ Composite Model consisting of Buses and Lines │
         │                  ╭───────╮                    │
         │              2 ┯━┿       ┿━┯ 3                │
         │                ↓ │   ╭───╯ ↓                  │
         │                  ┷━┯━┷ 1                      │
         │                   (~)                         │
         ╰───────────────────────┬───────────────────────╯
    PowerDynamics.jl constructs  ▾
                       ╭─────────┴─────────╮
                       │ DAE System        │
                       │ M ̇x = f(x, p, t)  │
                       ╰─────────┬─────────╯
     RHS function + Mass Matrix  ▾
         ╭───────────────────────┴───────────────────────╮
         │ SciML-ODEProblem                              │
         │ Data structure for time-domain simulation     │
         ╰───────────────────────┬───────────────────────╯
   Any OrdinaryDiffEq.jl solver  ▾
         ╭───────────────────────┴───────────────────────╮
         │ SciML-ODESolution                             │
         │ Solution object containing the time series    │
         │ for all components                            │
         ╰───────────────────────┬───────────────────────╯
              Symbolic Indexing  ▾
         ╭───────────────────────┴───────────────────────╮
         │ Time-series Inspection                        │
         │ Symbolic indexing allows for easy access to   │
         │ all states of all subcomponents for detailed  │
         │ analysis.                                     │
         ╰───────────────────────────────────────────────╯
```

PowerDynamics.jl gives you direct access to the underlying DAE structure and purposely exposes you to the "raw" commands from the SciML Ecosystem.
While we could provide wrapper functions like `run_powergrid_simulation` for common tasks, this would limit flexibility for advanced use cases. By working directly with the DAE system, you have full control to customize initialization, event handling, and solver configuration for your specific needs -- and you start learning them from the very first usage.

This tight integration with the normal SciML workflow also means that you can use any of the powerful tools from the SciML ecosystem to analyze and manipulate your powergrid models, e.g. for sensitivity analysis, parameter estimation and plotting.

In this tutorial, we will first show you how you'd solve a simple example system which has nothing
to do with networks or powergrids and just uses "plain" SciML Packages.
We will then build a super simple power grid simulation to highlight the parallels in nomenclature and workflow.

For a quick overview, here are the key packages that appear throughout the documentation:

**Top-level Packages:**
- [PowerDynamics.jl](https://github.com/JuliaEnergy/PowerDynamics.jl): The main package for building powergrid models. It provides a library and modeling tools specific to power systems, such as powerflow models and component libraries.
- [NetworkDynamics.jl](https://github.com/JuliaDynamics/NetworkDynamics.jl): Our backend package that provides most of the core functionality. It is general-purpose and can model any kind of networked dynamical system.

**SciML Packages:**
- [ModelingToolkit.jl (MTK)](https://github.com/SciML/ModelingToolkit.jl): A symbolic modeling framework for defining and manipulating differential equations. The key word here is *symbolically* – you write equations, not numerical code. MTK automatically performs simplifications and generates efficient numerical code for simulation.
- [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl): Umbrella package for everything related to differential equations, including stochastic and delay differential equations. Since it's large, we typically import specific subpackages, i.e.:
- [OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl): Solvers for ordinary differential equations (ODEs and DAEs). You can reduce load time even further by only importing specific solver packages like OrdinaryDiffEqRosenbrock.jl or OrdinaryDiffEqTsit5.jl.
- [NonlinearSolve.jl](https://github.com/SciML/NonlinearSolve.jl): Solvers for nonlinear systems of equations, used for powerflow calculations and DAE initialization.

**Other Packages:**
- [Makie.jl](https://github.com/MakieOrg/Makie.jl): A powerful plotting package for visualizing results with its backends CairoMakie.jl for vector graphic output and GLMakie.jl/WGLMakie.jl for interactive visualizations.


## Simple ModelingToolkit System
To get things started, we'll first show you how to define and simulate a simple ODE system using just ModelingToolkit.jl.

Since we really like oscillators, we'll use a simple pendulum as an example.
```asciiart
     ╶┬╴
      ┆╲
      ┆ ╲ l
      ┆  ╲
      ┆θ,ω̇╲
           ╲
            ● m
            ↓ g
```
whose equations of motion can be written as
```math
\begin{aligned}
F &= -m\,g\,\sin{\theta}&&\text{(tangential force)} \\
m\,l^2\,\dot{\omega} &= F\,l&&\text{(Newton's second law)}\\
\dot{\theta} &= \omega
\end{aligned}
```

To simulate this system, we first need to import some packages...
=#

using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as Dt
using OrdinaryDiffEqRosenbrock
using CairoMakie
nothing #hide #md
#=
## MTK: Model Definition & Simulation
After importing the packages, we can define the *symbolic system* using the `@mtkmodel` macro:
=#
@mtkmodel Pendulum begin
    @parameters begin
        g = 9.81  # gravitational acceleration
        m = 1 # mass of the pendulum
        l = 1 # length of the pendulum
    end
    @variables begin
        θ(t) # angle of the pendulum
        ω(t) # angular velocity of the pendulum
        F(t) # accelerating force in tangential direction
    end
    @equations begin
        F ~ -m * g * sin(θ)
        m * l^2 * Dt(ω) ~ F * l
        Dt(θ) ~ ω
    end
end
nothing #hide #md
#=
The definition of the system is quite straightforward. Note how we defined 3
states, including one for the force $F$.
We can instantiate the system by calling its constructor `Pendulum()`:
=#
@named symbolic_system = Pendulum()
full_equations(symbolic_system) # show all equations
#=
In order to simulate the system, we need to call `mtkcompile`, which will essentially
perform a symbolic simplification of the system:
=#
compiled_system = mtkcompile(symbolic_system)
full_equations(compiled_system) # show all equations
#=
You can see that the "compiled" system only consists of two states, $\theta$ and $\omega$.
This is because $F$ is not really a state of the system, but rather an intermediate variable.
The resulting equations are much more closely aligned with the "canonical" form of the pendulum.

To simulate the system, we need to define initial conditions for the states $\theta$ and $\omega$.
Also, we need to define a time span for the simulation.
=#

u0 = [
    compiled_system.θ => pi/4,
    compiled_system.ω => 0.0,
]
tspan = (0.0, 5.0)
nothing #hide #md
#=
Combining the compiled system, initial conditions, and time span, we can define a so-called `ODEProblem`.
=#
prob = ODEProblem(compiled_system, u0, tspan)

#=
The ODEProblem contains all the information needed to simulate the system.
We can simulate the system using any of the solvers from OrdinaryDiffEq.jl.
In this case, we decided to use the `Rodas5P` solver from `OrdinaryDiffEqRosenbrock.jl`.
=#

sol = solve(prob, Rodas5P())
nothing #hide #md

#=
## MTK: Solution Handling
The solution object we get contains all the time series in the system. For low-level access, we can look at
=#
sol.t
#=
to get an array of all the points in time the solver stepped to. While
=#
sol.u
#=
gives the full state of the system for each of the time points.

However, this is far from all we can do with the solution object!
First off, since we use dense output by default, we can interpolate the solution at any point in time:
=#
sol(2.5) # interpolate at t=2.5s
#=
This interpolation is **far more accurate** than a simple linear interpolation between the time points, because the solution object also contains the derivatives of the states at each time point.

The output, however, is still not very user friendly, since we only get a vector of values.
This is where **symbolic indexing** comes to our help!
Using the syntax
=#
sol(1.0, idxs=compiled_system.θ) # get θ at t=1.0s
#=
we can extract a specific state at a specific time point.
This syntax has lots of variants; for example, we can efficiently interpolate multiple states at multiple time points:
=#
sol([0.0, 1.0], idxs=[compiled_system.θ, compiled_system.ω]) # get θ and ω at t=0.0s and t=1.0s
#=
Since ModelingToolkit keeps track of all of its simplifications, we can also extract so-called "observed" states, i.e., states that were part of the original symbolic system but got eliminated during compilation.
For this example, we can get the force $F$ at any time point even though it is not part of the solution itself:
=#
sol(0.0, idxs=compiled_system.F) # get F at t=0.0s
#=
Finally, we can use the same *symbolic indexing* syntax in plotting commands. The example below uses Makie.jl; however, the commands are very similar in Plots.jl.
=#
let
    fig = Figure()
    ax = Axis(fig[1,1], xlabel="Time (s)", ylabel="States")
    lines!(sol, idxs=compiled_system.θ, color=:darkred)
    lines!(sol, idxs=compiled_system.ω, color=:darkblue)
    lines!(sol, idxs=compiled_system.F, color=:darkgreen)
    axislegend(ax; position=:rt)
    fig
end

#=
Now that we've seen how to work with ModelingToolkit systems, let's apply these same concepts to build a simple power grid simulation.

## Simple PowerDynamics System

Now that we have seen how to define and simulate a simple ODE system using ModelingToolkit.jl, we can move on to a simple powergrid example using PowerDynamics.jl.

As an example, we'll build a simple 2-bus system with two swing equation models at two buses connected by a transmission line.
```asciiart
   ╭─────────╮
1 ━┿━       ━┿━ 2
   │         │
  (~)       (~)

```
First, we need to load the PowerDynamics.jl package:
=#
using PowerDynamics
using PowerDynamics: Library
#=
We start by loading a Swing model generator from the library.
=#
@named symbolic_swing_system = Library.Swing(V=1)
full_equations(symbolic_swing_system) # show all equations

#=
The equations represent a classic swing equation with no voltage dynamics.
We passed the keyword argument `V=1` to set the voltage magnitude to 1 p.u.
So far, this is a "pure" MTK model, similar to the `symbolic_system` from above.

We can then compile the model to get rid of intermediate variables and make it ready for simulation. We do so by calling `compile_bus`. The additional call to `MTKBus` can be ignored for now and is explained in further tutorials.
Additionally, we give the bus model an index using the `vidx` keyword (short for vertex index).
=#
bus1 = compile_bus(MTKBus(symbolic_swing_system); vidx=1)
#=
This object is a so-called `VertexModel`. VertexModels (and EdgeModels) are the building blocks
of systems in PowerDynamics.jl and NetworkDynamics.jl.
From the printout you can already see that it has different variables/parameters with some default values and so on.

Similarly, we can create a second bus:
=#
bus2 = compile_bus(MTKBus(symbolic_swing_system); vidx=2)
#=
And a powerline connecting the two:
=#
@named symbolic_piline = Library.PiLine()
line = compile_line(MTKLine(symbolic_piline); src=1, dst=2)
#=
The powerline got the `src` and `dst` keywords. This means our line is defined from bus 1 to bus 2.

Having defined all the components, we can now connect them to a network model.
=#
nw = Network([bus1, bus2], line)
#=
The `nw` object is somewhat similar to the `compiled_system` from above: it is a fully defined
DAE system (ODE system in this case) that can be simulated.
Similar to the `compiled_system`, it not only contains the right-hand-side function but also contains information necessary for **symbolic indexing**, i.e., which component has which states/parameters under which names and so on.

## PD: Symbolic Indexing

In contrast to the MTK example above, our symbolic indices are "hierarchical", i.e., we have to specify the component first and then the state/parameter name.

Symbolic network indices follow a very simple pattern consisting of a "component" specifier and a "state" specifier:
```asciiart
      ▷ VIndex(    1, :symbolic_swing_system₊ω)
      ▷ VIndex(:bus1, :symbolic_swing_system₊δ)
        ╶──┬─╴ ╶─┬─╴  ╶───────────┬──────────╴
Vertex╶────╯     ╰──────╮         │
Index or name of vertex╶╯         │
Name of parameter/state╶──────────╯

      ▷ EIndex(             1, :src₊P)
      ▷ EIndex(        :edge1, :src₊Q)
      ▷ EIndex(        1 => 2, :src₊Q)
      ▷ EIndex(:bus1 => :bus2, :symbolic_piline₊R)
        ╶──┬─╴ ╶──────┬─────╴  ╶───────┬────────╴
Edge╶──────╯          ╰────╮           │
Index/name or src-dst pair╶╯           │
Name of parameter/state╶───────────────╯
```

In order to simulate the system, we need to define initial conditions and a time span.

## PD: Manual Definition of Initial Conditions

For large systems, we can have thousands of states and thousands of parameters.
Finding a suitable initial state, i.e., initializing the system, is quite complex in general and is covered in depth in later tutorials.

For now, our system is quite simple and we can find a suitable initial state by hand.
We can create a "default" state by calling [`NWState`](@extref NetworkDynamics.NWState) on the network object:
=#
s0 = NWState(nw)
#=
This creates a state and parameter object, which is prefilled with all of the default values stored in the Network.
The undefined states/parameters are set to `NaN`.
Using the symbolic indexing syntax described above, we can now set the initial conditions for all states and parameters that are not already defined.

We initialize the frequency of both machines to 1 p.u. and the angles to small deviations around 0:
=#
s0[VIndex(1, :symbolic_swing_system₊θ)] = 0.01
s0[VIndex(2, :symbolic_swing_system₊θ)] = -0.01
s0[VIndex(1, :symbolic_swing_system₊ω)] = 1
s0[VIndex(2, :symbolic_swing_system₊ω)] = 1
nothing #hide #md
#=
For the parameters, we set the mechanical power input of bus 1 to 1, and bus 2 to -1:
=#
s0[VIndex(1, :symbolic_swing_system₊Pm)] =  1
s0[VIndex(2, :symbolic_swing_system₊Pm)] = -1
nothing #hide #md

#=
We can also do more complex operations, for example, we'll set the inertia and damping of both machines to 1 using broadcasting:
=#
s0.v[1:2][:symbolic_swing_system₊M] .= 1 # sets the parameters :…₊M of vertices 1 and 2 to 1
s0.v[1:2][:symbolic_swing_system₊D] .= 1 # sets the parameters :…₊D of vertices 1 and 2 to 1
nothing #hide #md

#=
It is important to understand that at its core, `NWState` objects are just "wrappers" around flat arrays.
Similar to the pure-MTK example above, where our state vector `u` was just a vector of 2 plain values, the `NWState` object contains a flat vector of all states and a flat vector of all parameters. The flat vectors can be accessed using the [`uflat`](@extref NetworkDynamics.uflat) and [`pflat`](@extref NetworkDynamics.pflat) functions:
=#
uflat(s0)
#-
pflat(s0)
#=
By wrapping those flat vectors in a `NWState` object we make them "human readable" by providing symbolic indexing.

There are lots of things you can do with `NWState` objects.
For example, once again it is possible to inspect "observed" states—states which are not actually part of the state vector but rather intermediate variables.
For example, we can inspect the active power at both src and destination end.
=#
s0.e[1=>2]([:src₊P, :dst₊P])
#=
The exact syntax will be covered in later tutorials.
However, you can already see that the pi line injects positive active power at the destination end and negative active power at the source end. This matches our expectation, since the line is defined from bus 1 to bus 2, and we gave bus 1 a slightly higher angle than bus 2—therefore we expected this net flow from 1 to 2.

## PD: Simulation of the System
Similar to before, we take our model and use it to define an `ODEProblem`.
We can then solve it using the `Rodas5P` solver again.
=#
prob = ODEProblem(nw, uflat(s0), (0.0, 10.0), pflat(s0))
sol = solve(prob, Rodas5P())
nothing #hide #md

#=
## PD: Solution Handling
The solution handling is analogous to the pure-MTK example above.
=#
sol(1.0, idxs=VIndex(1, :symbolic_swing_system₊θ)) # get θ of bus 1 at t=1.0s
#=
For generating lists of symbolic indices at once, NetworkDynamics.jl provides the
auxiliary functions [`vidxs`](@extref NetworkDynamics.vidxs) and [`eidxs`](@extref NetworkDynamics.eidxs):
=#
sol(1.0, idxs=vidxs(nw, 1:2, :symbolic_swing_system₊ω)) # get ω of bus 1&2 at t=1.0s
#=
Sometimes, you want to get the full `NWState` at a specific time point.
=#
s10 = NWState(sol, 1.0) # get full NWState at t=1.0s
#=
Lastly, we can do some plotting as before:
=#
let
    fig = Figure()
    ax = Axis(fig[1,1], xlabel="Time (s)", ylabel="Angles")
    lines!(sol, idxs=VIndex(1, :symbolic_swing_system₊θ), color=:darkred)
    lines!(sol, idxs=VIndex(2, :symbolic_swing_system₊θ), color=:darkblue)
    axislegend(ax; position=:rt)
    ax = Axis(fig[2,1], xlabel="Time (s)", ylabel="Frequencies")
    lines!(sol, idxs=VIndex(1, :symbolic_swing_system₊ω), color=:darkred)
    lines!(sol, idxs=VIndex(2, :symbolic_swing_system₊ω), color=:darkblue)
    axislegend(ax; position=:rt)
    ax = Axis(fig[3,1], xlabel="Time (s)", ylabel="Active Power")
    lines!(sol, idxs=EIndex(1=>2, :src₊P), color=:darkgreen, label="P injected towards bus 1")
    lines!(sol, idxs=EIndex(1=>2, :dst₊P), color=:lightgreen, label="P injected towards bus 2")
    axislegend(ax; position=:rt)
    fig
end
#=
We observe the expected behavior:
- both swing generators, well, swing around each other until they reach a steady state (due to the non-zero damping term)
- in steady state, the active power injected at bus 1 is equal to the active power extracted at bus 2
- there is a net power flow from bus 1 to bus 2
=#
