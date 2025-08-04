#=
# Define a Custom Bus

In this Tutorial, we will define a custom bus model that can be used in PowerDynamics.jl.

The model we set out to recreate is the classical machine from Chapter 15.1 from Milanos book

> F. Milano, Power System Modelling and Scripting,  Berlin, Heidelberg: Springer Berlin Heidelberg, 2010. doi: 10.1007/978-3-642-13669-6.

## Defining the Machine as Injector

In order to use this model in a Bus, we need to define it in a way that it
specifies the [Injector Interface](@ref).

```
            ┌───────────────────┐
terminal    │                   │
   o←───────┤ Machine Equations │
u_r, u_i    │                   │
i_r, i_i    └───────────────────┘

```

The received values for $u_r$, $u_i$, $i_r$, and $i_i$ at the terminal are in the global
synchronous dq frame. The internal state $δ$ describs the rotor angle of the machin in this frame.
In order to obtain the **local** dq-frame voltages and currents, we need to apply a Park transformation.

```@raw html
<picture>
<source srcset="../../assets/dqgrafic-dark.svg" media="(prefers-color-scheme: dark)">
<img src="../../assets/dqgrafic.svg" width="70%" height="70%"/>
</picture>
```

Additional to the transformation, the mode is defined by the following equations:
```math
\begin{aligned}
\frac{d\delta}{dt} &= \omega_b(\omega - 1) &\text{(Milano 15.5)} \\
2H \frac{d\omega}{dt} &= \frac{P_m}{\omega} - \tau_e &\text{(Power form of Milano 15.5)} \\
\psi_d &= V_q + R_s I_q &\text{(Milano 15.11)} \\
\psi_q &= -V_d - R_s I_d &\text{(Milano 15.11)} \\
\tau_e &= \psi_d I_q - \psi_q I_d &\text{(Milano 15.6)} \\
0 &= V_q + R_s I_q + X'_d I_d - v_{f,\text{set}} &\text{(Milano 15.36)} \\
0 &= V_d + R_s I_d - X'_d I_q &\text{(Milano 15.36)}
\end{aligned}
```

We can use the ModelingToolkit DSL to define the full injector model:
=#
using PowerDynamics, NetworkDynamics, ModelingToolkit
using PowerDynamics.Library
using ModelingToolkit: t_nounits as t, D_nounits as Dt
using ModelingToolkitStandardLibrary.Blocks
using OrdinaryDiffEqRosenbrock, OrdinaryDiffEqNonlinearSolve
using CairoMakie

@mtkmodel MilanoClassicalMachine begin
    @components begin
        terminal=Terminal()
    end
    @parameters begin
        R_s=0.000124, [description="stator resistance"]
        X′_d=0.0608, [description="d-axis transient reactance"]
        H=23.64, [description="inertia constant"]
        ω_b=2π*50, [description="System base frequency in rad/s"]
        vf_set, [guess=1, description="field voltage"]
        P_m, [guess=1, description="mechanical power"]
    end
    @variables begin
        δ(t), [guess=0, description="rotor angle"]
        ω(t), [guess=1, description="rotor speed"]
        τ_e(t), [description="electrical torque"]
        I_d(t), [description="d-axis current"]
        I_q(t), [description="q-axis current"]
        V_d(t), [description="d-axis voltage"]
        V_q(t), [description="q-axis voltage"]
        ψ_d(t), [description="d-axis flux linkage"]
        ψ_q(t), [description="q-axis flux linkage"]
    end
    begin
        T_to_loc(α)  = [ sin(α) -cos(α);
                         cos(α)  sin(α)]
        T_to_glob(α) = [ sin(α)  cos(α);
                        -cos(α)  sin(α)]
    end
    @equations begin
        ## Park's transformations
        [terminal.u_r, terminal.u_i] .~ T_to_glob(δ)*[V_d, V_q]
        [I_d, I_q] .~ T_to_loc(δ)*[terminal.i_r, terminal.i_i]

        ## mechanical swing equation Milano 15.5
        Dt(δ) ~ ω_b*(ω - 1)
        2*H * Dt(ω) ~ P_m/ω - τ_e

        ## static flux linkage equations Milano 15.11
        ψ_d ~  V_q + R_s*I_q
        ψ_q ~ -V_d - R_s*I_d

        ## electrical torque Milano 15.6
        τ_e ~ ψ_d*I_q - ψ_q*I_d

        ## magnetic equations from static model Milano 15.36
        0 ~ V_q + R_s*I_q + X′_d*I_d - vf_set
        0 ~ V_d + R_s*I_d - X′_d*I_q
    end
end


@named machine = MilanoClassicalMachine();

#=
We can assure, that the model satisfies the [Injector Interface](@ref) by checking
=#
isinjectormodel(machine)

#=
## Attaching the Machine to a Busbar

In order to use the machine model, we need to attach it to a
busbar, thus forming a system which satisfies the [MTKBus Interface](@ref).
There are two ways of doing so: manually and using the `MTKBus` constructor.

**Manual Construction**

We need to define a new ODESystem, which has 2 components: a busbar and the machine.
Both components have a `terminal` as a subcomponent, we can use the `connect` function
to hook the machine on the busbar.
=#
@mtkmodel MyMTKBus begin
    @components begin
        busbar = BusBar()
        machine = MilanoClassicalMachine()
    end
    @equations begin
        connect(busbar.terminal, machine.terminal)
    end
end
mtkbus = MyMTKBus(name=:bus)
isbusmodel(mtkbus) # assert that the created model satisfies the interface

#=
**Automatic Construction**

We can also use the [`MTKBus`](@ref) constructor to create a busbar with a machine attached.
This constructor takes a list of *injector* models and hooks them all to the same busbar.
=#
mtkbus = MTKBus(machine; name=:bus)
isbusmodel(mtkbus) # assert that the created model satisfies the interface

#=
## Compiling bus to `VertexModel`
To actual simulate the system, we need to *compile* the model, i.e. transformig it from a
purly symbolic representation to a numerical one.
=#
Bus(mtkbus)

#=
## Defining a Simulation Scenario

To simulate the model, we need to define some kind of scenario. We'll create a simple
two-bus system where our custom Milano machine is connected to a slack bus through
a transmission line. This will allow us to observe the machine's dynamic behavior
in response to a frequency disturbance.

First, we create a slack bus that provides the voltage and frequency reference for the system.
We'll also add a frequency disturbance event to observe the machine's response.
=#
slackbus = Bus(
    PowerDynamics.VariableFrequencySlack(name=:variable_slack),
    vidx=1,
    pf=pfSlack(V=1)
)

#=
We define a frequency event that increases the system frequency at t=1 second.
This disturbance will cause our machine to respond dynamically as it tries to
maintain synchronism with the network.
=#
freq_event = PresetTimeComponentCallback(
    1,
    ComponentAffect([],[:V,:ω]) do u, p, ctx
        @info "Adapt frequency"
        p[:ω] = 1.01 # set frequency to 1.01 pu
        # p[:V] = 1.001
    end
)
set_callback!(slackbus, freq_event)

#=
Next, we create the generator bus using our custom Milano machine model.
We specify it as a PV bus in the power flow with 1 pu voltage and 1 pu active power.
=#
genbus = Bus(
    mtkbus,
    vidx=2,
    pf=pfPV(V=1, P=1)
)

#=
We connect the two buses with a simple PI transmission line model.
=#
line = Line(MTKLine(PiLine(; name=:piline)); src=1,dst=2)

#=
Now we can build the complete network with our two buses and the connecting line.
=#
nw = Network([slackbus, genbus], line)

#=
Before running dynamic simulation, we initialize the system from power flow.
This ensures that all dynamic states start from a steady-state condition.
=#
initialize_from_pf!(nw) # initialize the network and store the steady state in the components

#=
Let's examine the initial state of our generator bus to verify proper initialization.
=#
dump_initial_state(nw[VIndex(2)])

#=
## Dynamic Simulation

With the system properly initialized, we can set up and solve the dynamic simulation.
We simulate for 10 seconds to capture the machine's response to the frequency disturbance.
=#
s0 = NWState(nw)
prob = ODEProblem(nw, uflat(s0), (0,100), pflat(s0), callback=get_callbacks(nw))
sol = solve(prob, Rodas5P())
nothing #hide

#=
## Visualizing the Results

Now let's create comprehensive plots to visualize how our custom Milano machine
responds to the frequency disturbance. We'll plot several key variables that
demonstrate the machine's electromechanical dynamics.
=#

let
    fig = Figure(size=(800, 1000));

    ax1 = Axis(fig[1, 1];
        title="Rotor Angle",
        xlabel="Time [s]",
        ylabel="Rotor Angle δ [rad]")
    lines!(ax1, sol; idxs=VIndex(2, :machine₊δ), color=:blue, linewidth=2)
    vlines!(ax1, [1.0], color=:red, linestyle=:dash, alpha=0.7, label="Frequency disturbance")
    axislegend(ax1)

    ax2 = Axis(fig[2, 1];
        title="Rotor Speed",
        xlabel="Time [s]",
        ylabel="Rotor Speed ω [pu]")
    lines!(ax2, sol; idxs=VIndex(2, :machine₊ω), color=:green, linewidth=2)
    vlines!(ax2, [1.0], color=:red, linestyle=:dash, alpha=0.7)
    hlines!(ax2, [1.0], color=:black, linestyle=:dot, alpha=0.5, label="Synchronous speed")
    axislegend(ax2)

    ax3 = Axis(fig[3, 1];
        title="Machine Currents",
        xlabel="Time [s]",
        ylabel="Current [pu]")
    lines!(ax3, sol; idxs=VIndex(2, :machine₊I_d), label="I_d", color=:orange, linewidth=2)
    lines!(ax3, sol; idxs=VIndex(2, :machine₊I_q), label="I_q", color=:purple, linewidth=2)
    vlines!(ax3, [1.0], color=:red, linestyle=:dash, alpha=0.7)
    axislegend(ax3)

    ax4 = Axis(fig[4, 1];
        title="Machine Voltages",
        xlabel="Time [s]",
        ylabel="Voltage [pu]")
    lines!(ax4, sol; idxs=VIndex(2, :machine₊V_d), label="V_d", color=:cyan, linewidth=2)
    lines!(ax4, sol; idxs=VIndex(2, :machine₊V_q), label="V_q", color=:magenta, linewidth=2)
    vlines!(ax4, [1.0], color=:red, linestyle=:dash, alpha=0.7)
    axislegend(ax4)

    ax5 = Axis(fig[5, 1];
        title="Electrical Torque",
        xlabel="Time [s]",
        ylabel="Torque [pu]")
    lines!(ax5, sol; idxs=VIndex(2, :machine₊τ_e), color=:brown, linewidth=2)
    vlines!(ax5, [1.0], color=:red, linestyle=:dash, alpha=0.7)

    ax6 = Axis(fig[6, 1];
        title="Bus Power",
        xlabel="Time [s]",
        ylabel="Power [pu]")
    lines!(ax6, sol; idxs=VIndex(2, :busbar₊P), label="Active Power", color=:darkblue, linewidth=2)
    lines!(ax6, sol; idxs=VIndex(2, :busbar₊Q), label="Reactive Power", color=:darkred, linewidth=2)
    vlines!(ax6, [1.0], color=:red, linestyle=:dash, alpha=0.7)
    axislegend(ax6)
    save("machine_response.png", fig)
    fig
end |> display

#=
## Observing the Poor Damping Problem

From the plots above, we can see that the Milano classical machine exhibits very
lightly damped oscillations that persist for a very long time. The rotor angle
and speed oscillate for hundreds of seconds without settling to a steady state.

This poor damping behavior occurs because:
1. **No damper windings**: The model lacks electromagnetic damping mechanisms
2. **Constant field voltage**: No dynamic response to help stabilize the machine
3. **No mechanical damping**: The swing equation has no friction losses

To solve this problem, real power systems use control systems, particularly
**Power System Stabilizers (PSS)** that are specifically designed to damp
electromechanical oscillations.

## Adding a Power System Stabilizer

Let's create an improved machine model with controllable field voltage and
add the simplest possible PSS to demonstrate the damping improvement.
=#

#=
### Controllable Machine Model

First, we create a modified Milano machine with control inputs/outputs:
=#

@mtkmodel MilanoControllableMachine begin
    @components begin
        terminal=Terminal()
        # Control interface
        vf_in = RealInput(guess=1)  # field voltage input
        ω_out = RealOutput()        # frequency output for PSS
    end
    @parameters begin
        R_s=0.000124, [description="stator resistance"]
        X′_d=0.0608, [description="d-axis transient reactance"]
        H=23.64, [description="inertia constant"]
        ω_b=2π*50, [description="System base frequency in rad/s"]
        P_m, [guess=1, description="mechanical power"]
    end
    @variables begin
        δ(t), [guess=0, description="rotor angle"]
        ω(t), [guess=1, description="rotor speed"]
        τ_e(t), [description="electrical torque"]
        I_d(t), [description="d-axis current"]
        I_q(t), [description="q-axis current"]
        V_d(t), [description="d-axis voltage"]
        V_q(t), [description="q-axis voltage"]
        ψ_d(t), [description="d-axis flux linkage"]
        ψ_q(t), [description="q-axis flux linkage"]
    end
    begin
        T_to_loc(α)  = [ sin(α) -cos(α);
                         cos(α)  sin(α)]
        T_to_glob(α) = [ sin(α)  cos(α);
                        -cos(α)  sin(α)]
    end
    @equations begin
        ## Park's transformations
        [terminal.u_r, terminal.u_i] .~ T_to_glob(δ)*[V_d, V_q]
        [I_d, I_q] .~ T_to_loc(δ)*[terminal.i_r, terminal.i_i]

        ## mechanical swing equation Milano 15.5
        Dt(δ) ~ ω_b*(ω - 1)
        2*H * Dt(ω) ~ P_m/ω - τ_e

        ## static flux linkage equations Milano 15.11
        ψ_d ~  V_q + R_s*I_q
        ψ_q ~ -V_d - R_s*I_d

        ## electrical torque Milano 15.6
        τ_e ~ ψ_d*I_q - ψ_q*I_d

        ## magnetic equations from static model Milano 15.36
        0 ~ V_q + R_s*I_q + X′_d*I_d - vf_in.u  # Use controllable input
        0 ~ V_d + R_s*I_d - X′_d*I_q

        ## Control interface - output frequency for PSS
        ω_out.u ~ ω
    end
end

#=
### Simple Power System Stabilizer

The simplest PSS consists of a washout filter with gain. The washout filter
ensures the PSS only responds to frequency changes, not steady-state errors.
=#

@mtkmodel SimplePSS begin
    @components begin
        ω_in = RealInput()         # frequency input from machine
        vst = RealOutput()         # stabilizer output signal
    end
    @parameters begin
        Tw=10, [description="washout time constant"]
        Ks=10, [description="stabilizer gain"]
    end
    @variables begin
        x(t), [guess=0, description="washout filter state"]
    end
    @equations begin
        # Washout filter: removes DC, responds only to frequency changes
        Dt(x) ~ (ω_in.u - 1 - x) / Tw  # (ω_in.u - 1) gives frequency deviation
        vst.u ~ Ks * x
    end
end

#=
### Complete Generator with PSS

Now we create a composite model that combines the controllable machine with
the PSS, properly connecting the control loop.
=#

@mtkmodel GeneratorWithPSS begin
    @components begin
        terminal = Terminal()
        machine = MilanoControllableMachine()
        pss = SimplePSS()
    end
    @parameters begin
        vf_base, [guess=1.0, description="base field voltage"]
    end
    @equations begin
        # Connect terminals
        connect(terminal, machine.terminal)
        # Connect control loop: machine frequency → PSS → back to machine field voltage
        connect(machine.ω_out, pss.ω_in)
        # Sum base field voltage with PSS output
        machine.vf_in.u ~ vf_base + pss.vst.u
    end
end

@named gen_with_pss = GeneratorWithPSS()
isinjectormodel(gen_with_pss) # Verify it's still an injector


#=
## Simulation with PSS

Now let's run the same simulation scenario with the PSS-equipped generator
to observe the damping improvement.
=#

# Create the improved generator bus with simple PSS
genbus_pss = Bus(
    MTKBus(gen_with_pss; name=:bus_pss),
    vidx=2,
    pf=pfPV(V=1, P=1)
)

# Create network with PSS-equipped generator
nw_pss = Network([slackbus, genbus_pss], line)
initialize_from_pf!(nw_pss)

# Run simulation with simple PSS
s0_pss = NWState(nw_pss)
prob_pss = ODEProblem(nw_pss, uflat(s0_pss), (0,100), pflat(s0_pss), callback=get_callbacks(nw_pss))
sol_pss = solve(prob_pss, Rodas5P())

nothing #hide

#=
## Comparing Results: With and Without PSS

Let's create comparison plots to clearly see the damping improvement:
=#

let
    fig = Figure(size=(1000, 800));

    # Compare rotor speeds
    ax1 = Axis(fig[1, 1];
        title="Rotor Speed Comparison: Effect of PSS on Damping",
        xlabel="Time [s]",
        ylabel="Rotor Speed ω [pu]")
    lines!(ax1, sol; idxs=VIndex(2, :machine₊ω), label="No PSS (Unstable)", color=:red, linewidth=3)
    lines!(ax1, sol_pss; idxs=VIndex(2, :gen_with_pss₊machine₊ω), label="Simple PSS", color=:blue, linewidth=2)
    vlines!(ax1, [1.0], color=:black, linestyle=:dash, alpha=0.7, label="Frequency Disturbance")
    hlines!(ax1, [1.0], color=:gray, linestyle=:dot, alpha=0.5, label="Synchronous Speed")
    axislegend(ax1, position=:rt)
    xlims!(ax1, 0, 30)  # Focus on first 30 seconds

    # Compare rotor angles
    ax2 = Axis(fig[2, 1];
        title="Rotor Angle Comparison",
        xlabel="Time [s]",
        ylabel="Rotor Angle δ [rad]")
    lines!(ax2, sol; idxs=VIndex(2, :machine₊δ), label="No PSS", color=:red, linewidth=3)
    lines!(ax2, sol_pss; idxs=VIndex(2, :gen_with_pss₊machine₊δ), label="Simple PSS", color=:blue, linewidth=2)
    vlines!(ax2, [1.0], color=:black, linestyle=:dash, alpha=0.7)
    axislegend(ax2, position=:rt)
    xlims!(ax2, 0, 30)

    # Detailed view of first 10 seconds
    ax3 = Axis(fig[3, 1];
        title="Detailed View: First 10 Seconds After Disturbance",
        xlabel="Time [s]",
        ylabel="Rotor Speed ω [pu]")
    lines!(ax3, sol; idxs=VIndex(2, :machine₊ω), label="No PSS", color=:red, linewidth=3)
    lines!(ax3, sol_pss; idxs=VIndex(2, :gen_with_pss₊machine₊ω), label="Simple PSS", color=:blue, linewidth=2)
    vlines!(ax3, [1.0], color=:black, linestyle=:dash, alpha=0.7)
    hlines!(ax3, [1.0], color=:gray, linestyle=:dot, alpha=0.5)
    axislegend(ax3, position=:rt)
    xlims!(ax3, 0, 10)
    ylims!(ax3, 0.995, 1.015)

    save("pss_action.png", fig)
    fig
end |> display

#=
## Summary

This tutorial demonstrated the complete process of implementing and improving
a custom generator model in PowerDynamics.jl:

**Part 1: Basic Milano Classical Machine**
- Implemented Milano's classical machine equations from Chapter 15.1
- Showed how to use ModelingToolkit DSL for power system components
- Demonstrated proper Park transformations and interface compliance

**Part 2: Identifying the Stability Problem**
- Revealed poor damping characteristics typical of classical machine models
- Explained the physical reasons: no damper windings, constant field voltage, no mechanical damping
- Showed that oscillations persist for hundreds of seconds without control

**Part 3: Control System Solution**
- Created a controllable machine model with input/output interfaces
- Implemented the simplest possible PSS: washout filter + gain
- Used composite models to connect machine and PSS properly

**Part 4: Demonstrating the Improvement**
- Compared simulations with no PSS and simple PSS
- Showed dramatic damping improvement with properly tuned control systems
- Illustrated the importance of PSS parameter selection for stability

**Key Learning Points:**
- Classical machine models need control systems for realistic behavior
- PSS parameter tuning is critical: too aggressive gains can cause instability
- Simple PSS (washout + gain): Use moderate gains (Ks ≈ 5-15) and appropriate time constants (Tw ≈ 2s)
- ModelingToolkit allows easy creation of both component models and control systems
- Proper input/output interfaces enable flexible system composition

**PSS Design Guidelines:**
- Start with conservative parameters and increase damping gradually
- Washout time constant: 1-10 seconds (removes DC component)
- Simple PSS gain: 5-20 (POSITIVE for damping - increases field voltage when frequency rises)
- Key insight: Higher field voltage → higher electrical torque → slows down accelerating rotor

This example provides a foundation for implementing other custom generator models
and control systems in PowerDynamics.jl.
=#
