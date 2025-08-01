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
2H \frac{d\omega}{dt} &= \frac{\tau_m}{\omega} - \tau_e &\text{(Milano 15.5)} \\
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
using ModelingToolkit: t_nounits as t, D_nounits as Dt
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
        τ_m, [guess=1, description="mechanical torque"]
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
        2*H * Dt(ω) ~ τ_m / ω - τ_e

        ## static flux linkage equations Milaon 15.11
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

To simulate the model, we need to define some kind of scenario.
=#
slackbus = Bus(
    VariableFrequencySlack(name=:variable_slack),
    vidx=1,
    pf=pfSlack(V=1)
)
freq_event = PresetTimeComponentCallback(
    1,
    ComponentAffect([],[:ω]) do u, p, ctx
        @info "Adapt frequency"
        p[:ω] = 1.01 # set frequency to 1.01 pu
    end
)
set_callback!(slackbus, freq_event)

#
genbus = Bus(
    mtkbus,
    vidx=2,
    pf=pfPV(V=1, P=1)
)
#
line = Line(MTKLine(PiLine(; name=:piline)); src=1,dst=2)
nw = Network([slackbus, genbus], line)
initialize_from_pf!(nw) # initialize the network and store the steady state in the components

dump_initial_state(nw[VIndex(2)])

s0 = NWState(nw)
prob = ODEProblem(nw, uflat(s0), (0,10), pflat(s0), callback=get_callbacks(nw))
sol = solve(prob, Rodas5P())



#=
To see the model in action, we need to place it in some kind of scenario.
=#
