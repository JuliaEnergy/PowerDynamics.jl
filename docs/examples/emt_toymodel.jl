#=
# [EMT Toy Model Example](@id emt-toymodel)

This script can be downloaded as a normal Julia script [here](@__NAME__.jl). #md

This example demonstrates an electromagnetic transient (EMT) simulation of a simple
two-bus system using PowerDynamics.jl. The system consists of a slack bus connected
to a load bus through an RL transmission line, with the load bus having both a
dynamic PQ load and a capacitive shunt element.

We compare our simulation results with PowerFactory reference data to validate
the EMT modeling approach.

!!! note "Pedagogical Example"
    This is a **pedagogical example** that demonstrates the modeling concepts in
    PowerDynamics.jl are generally compatible with EMT simulations. However, this
    is far from being an actual interesting simulation study. The way we want to
    handle EMT simulations in PowerDynamics.jl is not yet fully clear and remains
    an active area of development.

    The example serves to illustrate the flexibility of the modeling framework
    rather than provide a production-ready EMT simulation tool.

## System Description

The test system includes:
- Bus 1: Slack bus (infinite bus with constant voltage)
- Bus 2: Load bus with PQ load and shunt capacitor
- Transmission line: RL model with distributed capacitance
- Disturbance: Load disconnection at t=0.1s
=#

using PowerDynamics
using PowerDynamics.Library
using NetworkDynamics
using ModelingToolkit
using ModelingToolkit: D_nounits as Dt, t_nounits as t
using CSV
using SteadyStateDiffEq
using OrdinaryDiffEqRosenbrock
using DataFrames
using CairoMakie

#=
## System Parameters

First, we define the base system parameters and component values.
=#

ω0    = 2π*50    ## Nominal frequency [rad/s]
Sbase = 300      ## Base power [MW]
Vbase = 110      ## Base voltage [kV]

Rline = 1        ## Line resistance [Ω]
Lline = (1/100π) ## Line inductance [H]
Cline = (2e-6)   ## Line capacitance [F]
Pload = -300     ## Load power (negative = consumption) [MW]

## Convert to per-unit values
Rline_pu = Rline / Zbase(Sbase, Vbase)
Lline_pu = Lline / Zbase(Sbase, Vbase)
Cline_pu = Cline / Ybase(Sbase, Vbase)
Pload_pu = Pload / Sbase
nothing # hide


#=
## Bus Definitions

### Slack Bus
The slack bus maintains constant voltage magnitude and angle,
representing an infinite bus or strong grid connection.
=#

slackbus = compile_bus(pfSlack(; V=1), vidx=1)

#=
### Dynamic Shunt Capacitor Model

The shunt capacitor is modeled as a dynamic component in the dq-frame.
This allows us to observe the three-phase voltages (`u_a`, `u_b`, `u_c`) by
transforming from the dq-frame back to abc coordinates.

The capacitor dynamics are given by:
```math
\begin{aligned}
\frac{du_r}{dt} &= \phantom{-}\omega_0 u_i - \frac{1}{C} i_r \\
\frac{du_i}{dt} &= -\omega_0 u_r - \frac{1}{C} i_i
\end{aligned}
```
where the $\omega_0 u_i$ and $-\omega_0 u_r$ terms account for the rotating dq-frame. 
The terminal currents $i_r$ and $i_i$ follow the injector current sign convention: a positive current is defined to flow out of the injector and into the terminal.

=#

@mtkmodel DynamicShunt begin
    @components begin
        terminal = Terminal()
    end
    @variables begin
        u_r(t), [guess=1, description="Real part of voltage"]
        u_i(t), [guess=0, description="Imaginary part of voltage"]
        ## Three-phase voltages as observables
        u_a(t), [description="Voltage in a phase"]
        u_b(t), [description="Voltage in b phase"]
        u_c(t), [description="Voltage in c phase"]
    end
    @parameters begin
        C, [description="Capacitance"]
        ω0, [description="Angular frequency of dq Frame"]
    end
    begin
        ## Transformation matrix from dq to abc coordinates
        Tdqinv(δ) = [cos(δ)       -sin(δ)
                     cos(δ-2pi/3) -sin(δ-2pi/3)
                     cos(δ+2pi/3) -sin(δ+2pi/3)]
    end
    @equations begin
        ## Capacitor dynamics in rotating dq-frame
        Dt(u_r) ~  ω0*u_i - 1/C * terminal.i_r
        Dt(u_i) ~ -ω0*u_r - 1/C * terminal.i_i
        ## Terminal connections
        terminal.u_r ~ u_r
        terminal.u_i ~ u_i
        ## Transform to three-phase voltages
        [u_a, u_b, u_c] ~ Tdqinv(ω0*t) * [u_r, u_i]
    end
end
nothing #hide #md


#=
### Load Bus Components

The load bus combines two components:
1. A PQ load consuming constant active power (injector model from Library)
2. A dynamic shunt capacitor representing line charging
=#

@named load = PQLoad(Pset=Pload_pu, Qset=0)
@named shunt = DynamicShunt(C=Cline_pu, ω0=ω0)
loadbus = compile_bus(
    MTKBus(load, shunt);
    vidx=2
)

#=
## Transmission Line Model

### Dynamic RL Branch

The transmission line is modeled as a dynamic RL branch in the dq-frame.
The line current dynamics are given by:
```math
\begin{aligned}
\frac{di_r}{dt} &= \phantom{-}\omega_0 i_i - \frac{R}{L} i_r + \frac{1}{L}(u_{\text{src},r} - u_{\text{dst},r}) \\
\frac{di_i}{dt} &= -\omega_0 i_r - \frac{R}{L} i_i + \frac{1}{L}(u_{\text{src},i} - u_{\text{dst},i})
\end{aligned}
```
where the voltage difference drives the current through the line impedance.
Similarly, the line model follows the injector interface at both ends: 
a positive current is defined as leaving the device and flowing toward the terminals.
=#

@mtkmodel DynamicRLBranch begin
    @components begin
        src = Terminal()
        dst = Terminal()
    end
    @variables begin
        i_r(t)=0, [description="Current in real part"]
        i_i(t)=-1, [description="Current in imaginary part"]
    end
    @parameters begin
        R, [description="Resistance"]
        L, [description="Inductance"]
        ω0, [description="Angular frequency of dq Frame"]
    end
    @equations begin
        ## RL line dynamics in rotating dq-frame
        Dt(i_r) ~  ω0 * i_i  - R/L * i_r + 1/L*(src.u_r - dst.u_r)
        Dt(i_i) ~ -ω0 * i_r  - R/L * i_i + 1/L*(src.u_i - dst.u_i)
        ## Terminal current connections (KCL enforcement)
        src.i_r ~ -i_r  ## Current flows out of source
        src.i_i ~ -i_i
        dst.i_r ~ i_r   ## Current flows into destination
        dst.i_i ~ i_i
    end
end

@named branch = DynamicRLBranch(; R=Rline_pu, L=Lline_pu, ω0=ω0)
line_model = compile_line(
    MTKLine(branch);
    src=1, dst=2
)

#=
## Network Assembly and Initialization

We assemble the complete network and attempt initialization.
=#

nw = Network([slackbus, loadbus], line_model)
try #hide #md
s0 = find_fixpoint(nw; alg=DynamicSS(Rodas5P()))
catch e #hide #md
    @error e #hide #md
end #hide #md

#=
### Initialization Challenge

The direct initialization fails due to the stiffness of the PQ load model.
When the load current is computed algebraically from $i = P \frac{u}{|u|^2}$,
the system becomes numerically challenging to initialize.

To overcome this, we use a "less stiff" load model with dynamics that
smooth out the algebraic singularity during initialization.

We create a "less stiff" version of the PQ load that introduces first-order
dynamics with a fast time constant (1/1000 s). This smooths the algebraic
relation and makes initialization more robust:

```math
\begin{aligned}
\frac{di_r}{dt} &= 1000 \left( P \frac{u_r}{u_r^2 + u_i^2} - i_r \right) \\
\frac{di_i}{dt} &= 1000 \left( P \frac{u_i}{u_r^2 + u_i^2} - i_i \right)
\end{aligned}
```

This approaches the algebraic PQ load behavior but avoids initialization issues.
=#

@mtkmodel LessStiffPQLoad begin
    @components begin
        terminal = Terminal()
    end
    @variables begin
        i_r(t)=0, [description="Current in real part"]
        i_i(t)=0, [description="Current in imaginary part"]
    end
    @parameters begin
        Pset, [description="Active Power demand"]
    end
    @equations begin
        ## First-order dynamics with fast time constant
        Dt(i_r) ~ 1e3*(-Pset * terminal.u_r/(terminal.u_r^2 + terminal.u_i^2) - i_r)
        Dt(i_i) ~ 1e3*(-Pset * terminal.u_i/(terminal.u_r^2 + terminal.u_i^2) - i_i)
        terminal.i_r ~ i_r
        terminal.i_i ~ i_i
    end
end

@named less_stiff_load = LessStiffPQLoad(Pset=-Pload_pu)
less_stiff_loadbus = compile_bus(
    MTKBus(less_stiff_load, shunt);
    vidx=2
)
less_stiff_nw = Network([slackbus, less_stiff_loadbus], line_model)
less_stiff_s0 = find_fixpoint(less_stiff_nw; alg=DynamicSS(Rodas5P()))

#=
### Initialize Target System

Perfect! The less stiff load initialization worked. Now we use this solution
as an initial guess for our target system with the algebraic PQ load.
=#

s0guess = NWState(nw)
## Transfer key state variables from less stiff solution
s0guess[VIndex(2, :busbar₊u_i)] = less_stiff_s0[VIndex(2, :busbar₊u_i)]
s0guess[VIndex(2, :busbar₊u_r)] = less_stiff_s0[VIndex(2, :busbar₊u_r)]
s0guess[EIndex(1, :branch₊i_i)] = less_stiff_s0[EIndex(1, :branch₊i_i)]
s0guess[EIndex(1, :branch₊i_r)] = less_stiff_s0[EIndex(1, :branch₊i_r)]
s0 = find_fixpoint(nw, s0guess; alg=DynamicSS(Rodas5P()))

#=
## Disturbance Setup

Excellent! The initialization succeeded. Now we set up a disturbance to
observe the system's transient response. We'll disable the load at t=0.1s
to simulate a sudden load disconnection.
=#

disable_load_affect = ComponentAffect([], [:load₊Pset]) do u, p, ctx
    println("Disabling load at time $(ctx.t)")
    p[:load₊Pset] = 0  ## Set load power to zero
end
set_callback!(loadbus, PresetTimeComponentCallback(0.1, disable_load_affect))
nothing #hide #md

#=
## Dynamic Simulation

With the system properly initialized and the disturbance configured,
we can now run the electromagnetic transient simulation.
=#

prob = ODEProblem(nw, s0, (0.0, 0.15))
sol = solve(prob, Rodas5P());
nothing #hide #md

#=
## Results and Validation

We compare our EMT simulation results with PowerFactory reference data
to validate the modeling approach. The comparison focuses on the three-phase
voltages at bus 2 during the load disconnection transient.

The thick gray lines show the PowerFactory reference, while our PowerDynamics.jl
results are overlaid in color. The close agreement validates our EMT modeling approach.
=#

fig = let
    fig = Figure()
    ax = Axis(fig[1,1];
        title="Three-Phase Voltage at Bus 2",
        xlabel="Time [s]",
        ylabel="Voltage [pu]")
    ts = range(0.09, 0.13; length=2000)

    ## Load PowerFactory reference data
    df = CSV.read(
        joinpath(pkgdir(PowerDynamics),"docs","examples", "emt_data_minimal.csv.gz"),
        DataFrame
    )
    ## Plot PowerFactory results (thick gray lines)
    lines!(ax, df.t, df.u_2_a; label="PowerFactory A", color=:lightgray, linewidth=5)
    lines!(ax, df.t, df.u_2_b; label="PowerFactory B", color=:lightgray, linewidth=5)
    lines!(ax, df.t, df.u_2_c; label="PowerFactory C", color=:lightgray, linewidth=5)

    ## Extract and plot our simulation results
    a = sol(ts, idxs=VIndex(2, :shunt₊u_a)).u
    b = sol(ts, idxs=VIndex(2, :shunt₊u_b)).u
    c = sol(ts, idxs=VIndex(2, :shunt₊u_c)).u
    lines!(ax, ts, a, label="PowerDynamics A", color=Cycled(1))
    lines!(ax, ts, b, label="PowerDynamics B", color=Cycled(2))
    lines!(ax, ts, c, label="PowerDynamics C", color=Cycled(3))

    axislegend(ax, position=:rt)
    xlims!(ax, ts[begin], ts[end])
    fig
end

#=
### Detailed View of Transient Response

Let's zoom in on the critical period around the load disconnection
to better observe the transient behavior and compare with the reference.
=#

xlims!(0.0995, 0.105)
fig
