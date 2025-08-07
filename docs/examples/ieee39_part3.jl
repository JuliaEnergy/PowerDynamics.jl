#=
# [IEEE39 Bus Tutorial - Part III: Dynamic Simulation](@id ieee39-part3)

This is the third part of the IEEE 39-bus tutorial series:

- **Part I: Model Creation** - Build the network structure with buses, lines, and components
- **Part II: Initialization** - Perform power flow calculations and dynamic initialization
- **Part III: Dynamic Simulation** (this tutorial) - Run time-domain simulations and analyze system behavior
- **Part IV: Advanced Modeling & Parameter Optimization** - Create custom components and optimize system parameters

In this tutorial, we'll demonstrate how to perform dynamic simulations of power system disturbances
using PowerDynamics.jl. We'll simulate a short circuit fault followed by line disconnection and
analyze the system's dynamic response.

## Short Circuit Disturbance Scenario

We will simulate a three-phase short circuit fault on a transmission line, which is a common
and severe disturbance in power systems. The disturbance scenario consists of:

1. **t = 0.1s**: A short circuit occurs on line 11, drastically reducing its impedance
2. **t = 0.2s**: The protective relay trips the line, completely disconnecting it from the system
3. **t = 0.2s onwards**: The system operates with the line permanently out of service

This scenario tests the system's ability to:
- Survive the initial fault (maintain synchronism)
- Stabilize after the line disconnection
- Operate reliably with reduced transmission capacity

This script can be downloaded as a normal Julia script [here](@__NAME__.jl). #md
=#

## Loading required packages and setup
using PowerDynamics
using PowerDynamics.Library
using ModelingToolkit
using NetworkDynamics
using OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqNonlinearSolve
using CairoMakie

## Load the network model from Part I
EXAMPLEDIR = joinpath(pkgdir(PowerDynamics), "docs", "examples")
include(joinpath(EXAMPLEDIR, "ieee39_part1.jl"))

#=
## Network Initialization

Before we can run dynamic simulations, we need to initialize the network as described in Part II.
This involves solving the power flow and initializing all dynamic components.

For buses with both generators and loads (buses 31 and 39), we need to add initialization
formulas to resolve structural underconstraints by setting the load voltage setpoint
equal to the bus voltage magnitude.
=#

## Add initialization formulas as described in Part II
formula = @initformula :ZIPLoad₊Vset = sqrt(:busbar₊u_r^2 + :busbar₊u_i^2)
set_initformula!(nw[VIndex(31)], formula)
set_initformula!(nw[VIndex(39)], formula)

## Initialize the complete network from power flow solution
s0 = initialize_from_pf!(nw; verbose=false)
nothing #hide

#=
## Short Circuit Disturbance Definition

To simulate realistic power system dynamics, we need to introduce a disturbance that will
excite the system's dynamic behavior. We'll simulate a short circuit fault on transmission
line 11, which connects buses 5 and 8.

### Understanding Line Models with Fault Capability

The transmission line models in our network include built-in parameters for fault simulation:

- `pibranch₊shortcircuit`: When set to 1, this simulates a three-phase to ground short circuit along the line. The position in percentage can be given as a parameter too.
- `pibranch₊active`: When set to 0, this completely disconnects the line from the network (no current flowing into line or out of line, i.e., the line is disconnected at both ends)

### Callback Functions for Disturbance Events

In order to simulate discrete perturbations, such as enabling a short circuit or disabling
a line, we need to use **callbacks**.
Callbacks are a neat [feature of DifferentialEquations.jl](https://docs.sciml.ai/DiffEqDocs/stable/features/callback_functions/),
which allow you to stop the solver under certain conditions and trigger a user-defined affect function to change the state of the system.

NetworkDynamics [inherits this functionality as well](@extref Callbacks).
In addition, ND.jl provides a new type of callback: component callbacks.
Those are callbacks which are attached to a single component rather than the full network,
bringing the effect handling to the component level.

Here we'll define two component callbacks: one to enable a short circuit at a given time and
one to disable the line at a given time.
=#

## Select the line to be affected by the short circuit
AFFECTED_LINE = 11
nothing #hide

#=
Let's examine the transmission line that will experience the short circuit. This line connects
two important buses in the network and its outage will test the system's stability.

Now we define callback functions to model the disturbance sequence:
=#

## Define callback to enable short circuit at t=0.1s
VERBOSE_CALLBACK = true #hide
_enable_short = ComponentAffect([], [:piline₊shortcircuit]) do u, p, ctx
    if VERBOSE_CALLBACK #hide
    @info "Short circuit activated on line $(ctx.src)→$(ctx.dst) at t = $(ctx.t)s"
    end #hide
    p[:piline₊shortcircuit] = 1
end
shortcircuit_cb = PresetTimeComponentCallback(0.1, _enable_short)

## Define callback to disconnect line at t=0.2s (fault clearing)
_disable_line = ComponentAffect([], [:piline₊active]) do u, p, ctx
    if VERBOSE_CALLBACK #hide
    @info "Line $(ctx.src)→$(ctx.dst) disconnected at t = $(ctx.t)s"
    end #hide
    p[:piline₊active] = 0
end
deactivate_cb = PresetTimeComponentCallback(0.2, _disable_line)

## Attach both callbacks to the selected line
set_callback!(nw, EIndex(AFFECTED_LINE), (shortcircuit_cb, deactivate_cb))
nothing #hide

#=
The callbacks are now attached to line 11. During simulation:
1. At t=0.1s, the short circuit callback activates, simulating the fault
2. At t=0.2s, the line disconnection callback activates, simulating relay action

Let's verify the callbacks are properly attached:
=#
nw[EIndex(AFFECTED_LINE)]

#=
## Dynamic Simulation

Now we're ready to perform the dynamic simulation. We'll set up and solve an ordinary
differential equation (ODE) problem that represents the network's dynamic behavior.

### Simulation Setup

The simulation process involves:
1. Creating the initial state vector from our initialized network
2. Setting up the ODE problem with appropriate time span and callbacks
3. Solving the ODE using a suitable numerical method
4. Analyzing the results

Note that we use [`get_callbacks`](@extref NetworkDynamics.get_callbacks-Tuple{Network}) to collect the
component callbacks, transform them into a `CallbackSet` compatible with the DifferentialEquations.jl ecosystem and pass them to the ODEProblem.
=#

u0 = NWState(nw) # state is stored in metadata because of mutating init function!
prob = ODEProblem(nw, uflat(u0), (0.0, 15.0), copy(pflat(u0)); callback=get_callbacks(nw))
## Solve the ODE using Rodas5P (suitable for stiff differential-algebraic systems)
sol = solve(prob, Rodas5P());
@assert SciMLBase.successful_retcode(sol) # ensure the simulation was successful

#=
The simulation is complete! The `sol` object contains the time-domain solution of all
network variables. We can now analyze how the system responded to the short circuit disturbance.

!!! tip
    The ODEProblem contains a reference to exactly one copy of the *flat parameter array*.
    If you use callbacks to change those parameters (as we do), it is advised to
    `copy` the parameter array before passing it to the ODEProblem! Otherwise the callback
    will change our `u0` object.
    Also, this means you need to be careful when using the same `prob` for multiple subsequent
    `solve` calls, as the initial state of the `prob` object might have changed!

### Simulation Results Overview

The solution object contains:
- Time-domain trajectories of all state variables (generator angles, voltages, etc.)
- Network interface variables (bus voltages, line currents and power flows)
- Full system response from t=0 to t=15 seconds

Let's examine some key aspects of the system response through various plots.
=#

#=
## Results Analysis

### Power Flow in the Affected Line

First, let's examine how the power flow through the faulted line changes during the disturbance.
This plot shows the most direct impact of our short circuit and line disconnection events.
=#

let fig = Figure(; size=(800, 400))
    ax = Axis(fig[1, 1];
        title="Active Power Flow in Line $AFFECTED_LINE During Short Circuit",
        xlabel="Time [s]",
        ylabel="Active Power [pu]")

    ## Focus on the disturbance period to see the fault clearly
    ts = range(0, 0.35, length=1000)

    ## Plot power flow at the destination end of the line
    lines!(ax, ts, sol(ts; idxs=EIndex(AFFECTED_LINE, :src₊P)).u;
           label="Active Power towards src", linewidth=2)
    lines!(ax, ts, sol(ts; idxs=EIndex(AFFECTED_LINE, :dst₊P)).u;
           label="Active Power towards dst", linewidth=2)

    axislegend(ax; position=:rt)

    fig
end

#=
**Observations:**
- **Normal operation (t < 0.1s)**: Steady power flow through the line. Power towards destination is positive while power towards source is negative. This means we have net transmission from bus 5 to bus 8.
- **Short circuit (0.1s < t < 0.2s)**: Dramatic power flow change due to the short circuit. Both source and destination show negative power, which means we have active power flowing from both sides towards the short circuit.
- **Line disconnection (t > 0.2s)**: Zero power flow as the line is permanently out of service

The power absorption during the fault demonstrates the severe electrical stress that short circuits
place on the system. The protective relay action at t=0.2s successfully isolates the fault.
=#

#=
### Voltage Response at Adjacent Buses

Next, let's examine how the buses directly connected to the faulted line respond to the disturbance.
These buses experience the most severe voltage impacts during the fault.
=#

let fig = Figure(; size=(800, 500))
    ax = Axis(fig[1, 1];
        title="Voltage Magnitudes at Buses Adjacent to Faulted Line",
        xlabel="Time [s]",
        ylabel="Voltage Magnitude [pu]")

    ## Full simulation time to see both disturbance and recovery
    ts = range(0, 15, length=1000)

    ## Get the source and destination buses of the affected line
    src_bus, dst_bus = get_graphelement(nw[EIndex(AFFECTED_LINE)])

    ## Plot voltage magnitudes at both ends of the faulted line
    lines!(ax, ts, sol(ts; idxs=VIndex(src_bus, :busbar₊u_mag)).u;
           label="Bus $src_bus (source)", linewidth=2)
    lines!(ax, ts, sol(ts; idxs=VIndex(dst_bus, :busbar₊u_mag)).u;
           label="Bus $dst_bus (destination)", linewidth=2)

    axislegend(ax; position=:rb)
    ylims!(ax, 0.85, 1.15)

    fig
end

#=
**Observations:** Both buses experience voltage depression during the short circuit but are able to recover after the short circuit is cleared by disconnection of the line.

The voltage recovery demonstrates the system's ability to adapt to the new network topology
after the line outage. The generators' automatic voltage regulators help maintain voltage stability.
=#

#=
### System-Wide Voltage Response

To get a complete picture of the system's response, let's examine the voltage profiles
across all 39 buses. This "spaghetti plot" shows how the disturbance propagates through
the entire network.
=#

let fig = Figure(; size=(800, 600))
    ax = Axis(fig[1, 1];
        title="Voltage Magnitudes Across All 39 Buses",
        xlabel="Time [s]",
        ylabel="Voltage Magnitude [pu]")

    ## Full simulation time range
    ts = range(0, 15, length=1000)

    ## Plot voltage magnitude for all buses
    for i in 1:39
        voltage_data = sol(ts; idxs=VIndex(i, :busbar₊u_mag)).u
        lines!(ax, ts, voltage_data; linewidth=2)
    end
    ylims!(ax, 0.85, 1.15)
    fig
end

#=
Once again, we see how all bus voltages are affected by the short circuit and the overall voltage drops. However, after the fault is cleared, the system achieves a synchronous steady-state again.
=#
