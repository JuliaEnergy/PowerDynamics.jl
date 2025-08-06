#=
# [IEEE39 Bus Tutorial - Part IV: Advanced Modeling & Parameter Optimization](@id ieee39-part4)

This is the fourth and final part of the IEEE 39-bus tutorial series:

- **Part I: Model Creation** - Build the network structure with buses, lines, and components
- **Part II: Initialization** - Perform power flow calculations and dynamic initialization
- **Part III: Dynamic Simulation** - Run time-domain simulations and analyze system behavior
- **Part IV: Advanced Modeling & Parameter Optimization** (this tutorial) - Create custom components and optimize system parameters

In this tutorial, we'll demonstrate advanced PowerDynamics.jl capabilities by:
1. Creating a custom droop-controlled inverter component
2. Integrating it into the IEEE 39-bus system
3. Optimizing its parameters to improve system performance

This tutorial showcases custom component creation and
the integration with Julia's optimization ecosystem for parameter tuning.

This script can be downloaded as a normal Julia script [here](@__NAME__.jl). #md

!!! note
    This tutorial is designed as a pedagogical example. It does not necessarily represent
    a realistic power system model and analysis, but rather serves to demonstrate the available tools
    while remaining relatively simple and concise.
=#

## Loading required packages and setup
using PowerDynamics
using PowerDynamics.Library
using ModelingToolkit
using ModelingToolkit: D_nounits as Dt, t_nounits as t
using NetworkDynamics
using NetworkDynamics: SII
using OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqNonlinearSolve
using SciMLSensitivity
using Optimization
using OptimizationOptimisers
using CairoMakie
using LinearAlgebra
using Graphs
using SparseConnectivityTracer
using Sparspak

## Load the network models from previous parts
EXAMPLEDIR = joinpath(pkgdir(PowerDynamics), "docs", "examples")
include(joinpath(EXAMPLEDIR, "ieee39_part1.jl"))  # Creates the basic network model
include(joinpath(EXAMPLEDIR, "ieee39_part3.jl"))  # Provides initialized network and reference solution
nothing # hide

#=
## Integration of a Droop-Controlled Inverter

In this section, we'll modify our network by adding a droop-controlled inverter.

### Mathematical Background

The droop-controlled inverter implements a decentralized control strategy commonly used in
microgrids and renewable energy integration. It establishes the following relationships:

**Power Measurement:**
```math
\begin{aligned}
P_{meas} &= u_r \cdot i_r + u_i \cdot i_i\\
Q_{meas} &= u_i \cdot i_r - u_r \cdot i_i
\end{aligned}
```

**Power Filtering (Low-pass filtering for measurement noise reduction):**
```math
\begin{aligned}
\tau \cdot \frac{dP_{filt}}{dt} &= P_{meas} - P_{filt} \\
\tau \cdot \frac{dQ_{filt}}{dt} &= Q_{meas} - Q_{filt}
\end{aligned}
```

**Droop Control:**
```math
\begin{aligned}
\omega &= \omega_0 - K_p \cdot (P_{filt} - P_{set}) \\
V &= V_{set} - K_q \cdot (Q_{filt} - Q_{set})
\end{aligned}
```

**Voltage Angle Dynamics:**
```math
\frac{d\delta}{dt} = \omega - \omega_0
```

**Output Voltage:**
```math
\begin{aligned}
u_r &= V \cdot \cos(\delta) \\
u_i &= V \cdot \sin(\delta)
\end{aligned}
```

These equations implement:
- **Frequency-Active Power Coupling (f-P)**: Frequency decreases when active power exceeds setpoint
- **Voltage-Reactive Power Coupling (V-Q)**: Voltage decreases when reactive power exceeds setpoint

This mimics the natural behavior of synchronous generators and enables stable power sharing
in islanded operation.

### Definition of the Droop Inverter Component

Network components in PowerDynamics must follow the [Injector Interface](@ref) - they connect to
the network through a single `Terminal`:

```
      ┌──────────────────────────┐
(t)   │                          │
 o←───┤ Droop Inverter Equations │
      │                          │
      └──────────────────────────┘
```

We can define the following MTKModel to represent the droop inverter:
=#

@mtkmodel DroopInverter begin
    @components begin
        terminal = Terminal()
    end

    @parameters begin
        Pset, [description="Active power setpoint", guess=1]
        Qset, [description="Reactive power setpoint", guess=0]
        Vset, [description="Voltage magnitude setpoint", guess=1]
        ω₀=1, [description="Nominal frequency"]
        Kp=1, [description="Active power droop coefficient"]
        Kq=0.1, [description="Reactive power droop coefficient"]
        τ_p = 1, [description="Active Power filter time constant"]
        τ_q = 1, [description="Reactive Power filter time constant"]
    end

    @variables begin
        Pmeas(t), [description="Measured active power", guess=1]
        Qmeas(t), [description="Measured reactive power", guess=0]
        Pfilt(t), [description="Filtered active power", guess=1]
        Qfilt(t), [description="Filtered reactive power", guess=1]
        ω(t), [description="Frequency"]
        δ(t), [description="Voltage angle", guess=0]
        V(t), [description="Voltage magnitude"]
    end

    @equations begin
        ## Power measurement from terminal quantities
        Pmeas ~ terminal.u_r*terminal.i_r + terminal.u_i*terminal.i_i
        Qmeas ~ terminal.u_i*terminal.i_r - terminal.u_r*terminal.i_i

        ## First-order low-pass filtering
        τ_p * Dt(Pfilt) ~ Pmeas - Pfilt
        τ_q * Dt(Qfilt) ~ Qmeas - Qfilt

        ## Droop control equations
        ω ~ ω₀ - Kp * (Pfilt - Pset)  # Frequency decreases with excess power
        V ~ Vset - Kq * (Qfilt - Qset)  # Voltage decreases with excess reactive power

        ## Voltage angle dynamics
        Dt(δ) ~ ω - ω₀

        ## Output voltage components
        terminal.u_r ~ V*cos(δ)
        terminal.u_i ~ V*sin(δ)
    end
end;

#=
### Creating a Bus with the Droop Inverter

Following the descriptions in [Modeling Concepts](@ref), we build an MTKBus using
the droop as the single injector and then compile the bus model, similar to
how we define the templates in [part I](@ieee39-part1):

```
           ╔═════════════════════════╗
           ║ Droop (compiled)        ║
 Network   ║  ┌────────────────────┐ ║
interface  ║  │ MTKBus             │ ║
 current ────→│┌──────┐ ┌────────┐ │ ║
           ║  ││BusBar├o┤Inverter│ │ ║
 voltage ←────│└──────┘ └────────┘ │ ║
           ║  └────────────────────┘ ║
           ╚═════════════════════════╝
```
=#

@named inverter = DroopInverter()
mtkbus = MTKBus(inverter)
droop_bus_template = Bus(mtkbus; name=:DroopInverter)

#=
We see that the droop inverter has 3 free parameters (you can check `free_p(droop_bus_template)` or `dump_initial_state(droop_bus_template)`).
Therefore, similar to what we did in [Part II](@ref ieee39-part2), we need to help the initialization
by attaching an additional initialization formula to the bus:
=#
set_initformula!(
    droop_bus_template,
    @initformula(:inverter₊Vset = sqrt(:busbar₊u_r^2 + :busbar₊u_i^2))
)
nothing # hide

#=
### Network Modification with Droop Inverter

We'll replace bus 32 (originally a controlled generator bus) with our new droop inverter bus.
=#
DROOP_BUS_IDX = 32
nothing #hide
#=
To do so, we first collect all the "old" vertex and edge models.

!!! tip
    We copy the components, to create individual instances for the new network.
    Since all the metadata (like the default values, initialize values and so on) are
    stored in the component models, we would otherwise share metadata between the
    old and the new network, which can lead to unexpected results.
=#

vertex_models = [copy(nw[VIndex(i)]) for i in 1:nv(nw)];
edge_models = [copy(nw[EIndex(i)]) for i in 1:ne(nw)];
nothing # hide

#=
Now we need to replace the original bus model at index `DROOP_BUS_IDX` with our new droop inverter bus.
However, we don't want to lose the original power flow model associated with this bus, so we
need to attach it to the droop bus model:
=#
original_pfmodel = get_pfmodel(vertex_models[DROOP_BUS_IDX])
nothing #hide

#=
We can then use the `Bus` constructor to essentially copy the droop_bus_template
and adjust some properties, like the powerflow model and the vertex index.
=#

droop_bus = Bus(droop_bus_template; pf=original_pfmodel, vidx=DROOP_BUS_IDX)

#=
We then replace the original bus model in the array with our droop bus and build a network again:
=#
vertex_models[DROOP_BUS_IDX] = droop_bus
nw_droop = Network(vertex_models, edge_models)
set_jac_prototype!(nw_droop; remove_conditions=true)
#=
Additionally, we've set the jacobian prototype for performance gains during
simulation and optimization, see NetworkDynamics docs on
[Sparsity Detection](@extref).
=#

#=
### Network Initialization with Droop Inverter

The modified network requires the same initialization steps as the original:
1. Power flow solution
2. Dynamic component initialization

This all happens within [`initialize_from_pf!`](@ref):
=#
s0_droop = initialize_from_pf!(nw_droop; verbose=false)
nothing #hide

#=
Let's examine the initialized state of our droop inverter:
=#
dump_initial_state(nw_droop[VIndex(DROOP_BUS_IDX)]; obs=false)
#=
We see that the filtered powers match the setpoints (steady state),
and both $P$ and $V_\mathrm{set}$ are initialized according to the parameters
of the PV powerflow model.
=#

#=
### Simulation with Droop Inverter

Now we'll simulate the modified network and compare it with the original system response:
=#
prob_droop = ODEProblem(nw_droop, uflat(s0_droop), (0.0, 15.0), copy(pflat(s0_droop));
                       callback=get_callbacks(nw_droop))
sol_droop = solve(prob_droop, Rodas5P())
@assert SciMLBase.successful_retcode(sol_droop)

#=
### Comparison of System Responses

Let's compare how the droop inverter affects the system's response to the short circuit disturbance:
=#

let fig = Figure(; size=(1000, 600))
    selected_buses = [3, 4, 25, DROOP_BUS_IDX]  # Representative buses including the droop bus
    ts = range(0, 10, length=1000)

    for (i, bus) in enumerate(selected_buses)
        row, col = divrem(i-1, 2) .+ (1, 1)
        ax = Axis(fig[row, col];
                  title="Voltage Magnitude at Bus $bus" * (bus==DROOP_BUS_IDX ? " (droop bus)" : ""),
                  xlabel="Time [s]",
                  ylabel="Voltage [pu]")

        ## Original system response
        lines!(ax, ts, sol(ts; idxs=VIndex(bus, :busbar₊u_mag)).u;
               label="Original System", color=:blue, linewidth=2)

        ## Droop inverter system response
        lines!(ax, ts, sol_droop(ts; idxs=VIndex(bus, :busbar₊u_mag)).u;
               label="With Droop Inverter", color=:red, linewidth=2)

        ylims!(ax, 0.85, 1.15)
        i == 1 && axislegend(ax; position=:rb)
    end

    fig
end

#=
We see that the overall system reacts similarly but distinctly differently to the identical disturbance.

## Parameter Optimization

To showcase advanced capabilities of the SciML-ecosystem and the integration with PowerDynamics.jl,
we now want to try to *tune* the droop inverter parameters so that the overall system behavior
more closely resembles the original behavior, i.e., to reduce the difference between the system with generator
and the system with droop inverter.

### Optimization Problem Formulation

We define a loss function that measures the deviation between the original system response
and the modified system response:

```math
L(p) = \sum_{i,t} |x_{ref}(t)_i - x(t;p)_i|^2
```

Where we have
- Parameters $p = [K_p, K_q, \tau]$ to be optimized
- the reference solution $x_{ref}(t)$ (original system)
- the solution of the modified system $x(t;p)$ with updated parameters $p$

**Goal:** Find parameters $p$ that minimize this loss function, making the droop inverter system
behave as closely as possible to the original system.

### Setting Up the Optimization

First, we define the reference solution and identify the tunable parameters.

We probe the original solution at fixed timepoints, exporting `u_r` and `u_i` for every
bus:
=#
opt_ref = sol(0.3:0.1:10, idxs=[VIndex(1:39, :busbar₊u_r), VIndex(1:39, :busbar₊u_i)])
nothing #hide

#=
Next, we need to identify the "tunable" parameters. This is a bit tricky, because
the overall `nw_droop` has 1271 parameters, so we need to find the indices of the
parameters we want to tune in the flat parameter array.
We can do so, by leveraging NetworkDynamics implementation of the [SymbolicIndexingInterface](https://github.com/SciML/SymbolicIndexingInterface.jl):
=#
# tunable_parameters = [:inverter₊Kp, :inverter₊Kq, :inverter₊τ]
tunable_parameters = [:inverter₊Kp, :inverter₊Kq, :inverter₊τ_p, :inverter₊τ_q]
tp_idx = SII.parameter_index(sol_droop, VIndex(DROOP_BUS_IDX, tunable_parameters))

#=
We also get their initial values, which we use as the starting point for the optimization.
=#
p0 = sol_droop(sol_droop.t[begin], idxs=collect(VIndex(DROOP_BUS_IDX, tunable_parameters)))

#=
### Loss Function Implementation

The loss function simulates the system with given parameters and compares the result
to the reference:
=#
function loss(p)
    ## Create parameter vector for the full system
    ## allp = similar(p, pdim(nw_droop)) # create a vector of the full length
    allp = similar(p, length(s0_droop.p)) # create a vector of the full length
    allp .= pflat(s0_droop.p) # copy all "initial" parameters to that vector
    allp[tp_idx] .= p  # Update only the tunable parameters with the parameters for the given optimization iteration

    ## Solve the system with new parameters
    _sol = solve(prob_droop, Rodas5P(autodiff=true);
        p = allp,
        saveat = opt_ref.t,
        tspan=(0.0, opt_ref.t[end]),
        initializealg = SciMLBase.NoInit(),
        abstol=0.01,
        reltol=0.01
    )

    ## Return infinite loss if simulation failed
    if !SciMLBase.successful_retcode(_sol)
        @warn "Retcode $(_sol.retcode) indicates a failed simulation, returning Inf loss"
        return Inf
    end

    ## Extract solution at reference time points
    x = _sol(opt_ref.t; idxs=[VIndex(1:39, :busbar₊u_r), VIndex(1:39, :busbar₊u_i)])

    ## Compute L2 norm of the difference
    res = opt_ref.u - x.u
    l2loss = sum(abs2, reduce(vcat, res))
end
nothing #hide

#=
### Optimization Execution

We use the Optimization.jl ecosystem with the Adam optimizer:
=#

## Create optimization function with automatic differentiation
optf = Optimization.OptimizationFunction((x, p) -> loss(x), Optimization.AutoForwardDiff())
nothing # hide

#=
To better monitor the optimization progress, we want to store the optimized parameters
at every iteration of the optimizer. We can do so by defining a callback function
for the optimizer:
=#
optimization_states = Any[] # global variable to store the optimization parameters at each step
callback = function (state, l)
    push!(optimization_states, state)
    println("Iteration $(state.iter): loss = $l\t p = $(state.u)")
    return false  # Continue optimization
end
nothing #hide
#=
That callback will snapshot the current parameter values at every step of the gradient descent.

With that, we can run the optimization:
=#
optprob = Optimization.OptimizationProblem(optf, p0; callback)
VERBOSE_CALLBACK = false #hide

@time optsol = Optimization.solve(optprob, Optimisers.Adam(0.06), maxiters = 50)

println("\nOptimization completed!")
println("Initial parameters: ", p0)
println("Optimized parameters: ", optsol.u)
println("Initial loss: ", loss(p0))
println("Final loss: ", loss(optsol.u))
nothing #hide

#=
### Optimization Results Analysis

Let's visualize how the optimization improved the system behavior:
=#

function plot_optimization_comparison(p_initial, p_current)
    fig = Figure(; size=(1200, 800))
    selected_buses = [3, 4, 25, DROOP_BUS_IDX]
    ts = range(0, 10, length=1000)

    ## Simulate with optimized parameters
    allp_opt = @lift let
        _p = copy(pflat(s0_droop))
        _p[tp_idx] .= $p_current
        _p
    end
    sol_opt = @lift solve(prob_droop, Rodas5P(); p=$allp_opt)


    for (i, bus) in enumerate(selected_buses)
        row, col = divrem(i-1, 2) .+ (1, 1)
        ax = Axis(fig[row, col];
                  title="Voltage Magnitude at Bus $bus" * (bus==DROOP_BUS_IDX ? " (droop bus)" : ""),
                  xlabel="Time [s]",
                  ylabel="Voltage [pu]")

        ## Reference (original system)
        lines!(ax, ts, sol(ts; idxs=VIndex(bus, :busbar₊u_mag)).u;
               label="Reference", linestyle=:solid, color=:blue, linewidth=2)

        ## Initial droop parameters
        lines!(ax, ts, sol_droop(ts; idxs=VIndex(bus, :busbar₊u_mag)).u;
               label="Initial Droop", linestyle=:dash, color=:red, linewidth=2)

        ## Optimized droop parameters
        dat = @lift $(sol_opt)(ts; idxs=VIndex(bus, :busbar₊u_mag)).u
        lines!(ax, ts, dat; label="Optimized Droop", color=:green, linewidth=2)

        ylims!(ax, 0.85, 1.15)
        i == 1 && axislegend(ax; position=:rb)
    end

    fig
end

pobs = Observable(p0)
comparison_fig = plot_optimization_comparison(p0, pobs)
record(comparison_fig, "parameter_evolution.mp4", optimization_states; framerate=10) do s
    pobs[] = s.u
end
nothing #hide

#=
![timeseries evolution animation](parameter_evolution.mp4)
=#

#=
### Parameter Evolution During Optimization

Let's see how each parameter changed during the optimization process:
=#
let fig = Figure(; size=(1000, 800))
    param_names = ["Kp", "Kq", "τ_p", "τ_q"]
    for (i, param_name) in enumerate(param_names)
        row = (i-1) ÷ 2 + 1
        col = (i-1) % 2 + 1
        ax = Axis(fig[row, col];
                  title="Parameter Evolution: $param_name",
                  xlabel="Iteration",
                  ylabel="Parameter Value")

        ## Extract parameter values over iterations
        param_values = [state.u[i] for state in optimization_states[1:end-1]]
        iterations = [state.iter for state in optimization_states[1:end-1]]

        scatterlines!(ax, iterations, param_values; linewidth=3, markersize=6, color=:blue)

        ## Mark initial and final values
        hlines!(ax, [p0[i]]; linestyle=:dash, color=:gray, alpha=0.7)
        text!(ax, 1, p0[i]; text="Initial: $(round(p0[i], digits=3))",
              fontsize=10, color=:gray)

        lossax = Axis(fig[row, col],
            yticklabelcolor=:black,
            yaxisposition = :right,
            ylabel="loss", yscale=log10,
            xgridvisible=false, ygridvisible=false
        )
        scatterlines!(
            lossax,
            iterations,
            [loss(s.u) for s in optimization_states[1:end-1]],
            color=:black, linewidth=1, markersize=3,
        )
    end
    fig
end
