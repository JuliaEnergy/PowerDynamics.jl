# # Define a Custom Bus Model
#
# In this Tutorial, we will define a custom bus model that can be used in PowerDynamics.jl.
#
# The model we set out to recreate is the classical machine from Chapter 15.1 from Milano's book
#
# > F. Milano, Power System Modelling and Scripting,  Berlin, Heidelberg: Springer Berlin Heidelberg, 2010. doi: 10.1007/978-3-642-13669-6.
#
# ## Defining the Machine as Injector
#
# In order to use this model in a Bus, we need to define it in a way that it
# specifies the Injector Interface.
#
# ```
#             ┌───────────────────┐
# terminal    │                   │
#    o←───────┤ Machine Equations │
# u_r, u_i    │                   │
# i_r, i_i    └───────────────────┘
#
# ```
#
# The received values for $u_r$, $u_i$, $i_r$, and $i_i$ at the terminal are in the global
# synchronous dq frame. The internal state $δ$ describes the rotor angle of the machine in this frame.
# In order to obtain the **local** dq-frame voltages and currents, we need to apply a Park transformation.
#
# ```@raw html
# <picture>
# <source srcset="../../assets/dqgrafic-dark.svg" media="(prefers-color-scheme: dark)">
# <img src="../../assets/dqgrafic.svg" width="70%" height="70%"/>
# </picture>
# ```
#
# In addition to the transformation, the model is defined by the following equations:
# ```math
# \begin{aligned}
# \frac{d\delta}{dt} &= \omega_b(\omega - 1) &\text{(Milano 15.5)} \\
# 2H \frac{d\omega}{dt} &= \frac{P_m}{\omega} - \tau_e &\text{(Power form of Milano 15.5)} \\
# \psi_d &= V_q + R_s I_q &\text{(Milano 15.11)} \\
# \psi_q &= -V_d - R_s I_d &\text{(Milano 15.11)} \\
# \tau_e &= \psi_d I_q - \psi_q I_d &\text{(Milano 15.6)} \\
# 0 &= V_q + R_s I_q + X'_d I_d - v_{f,\text{set}} &\text{(Milano 15.36)} \\
# 0 &= V_d + R_s I_d - X'_d I_q &\text{(Milano 15.36)}
# \end{aligned}
# ```
#
# We can use the ModelingToolkit DSL to define the full injector model:

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
        # Park's transformations
        [terminal.u_r, terminal.u_i] .~ T_to_glob(δ)*[V_d, V_q]
        [I_d, I_q] .~ T_to_loc(δ)*[terminal.i_r, terminal.i_i]

        # mechanical swing equation Milano 15.5
        Dt(δ) ~ ω_b*(ω - 1)
        2*H * Dt(ω) ~ P_m/ω - τ_e

        # static flux linkage equations Milano 15.11
        ψ_d ~  V_q + R_s*I_q
        ψ_q ~ -V_d - R_s*I_d

        # electrical torque Milano 15.6
        τ_e ~ ψ_d*I_q - ψ_q*I_d

        # magnetic equations from static model Milano 15.36
        0 ~ V_q + R_s*I_q + X′_d*I_d - vf_set
        0 ~ V_d + R_s*I_d - X′_d*I_q
    end
end


@named machine = MilanoClassicalMachine();

# We can verify that the model satisfies the Injector Interface by checking

isinjectormodel(machine)

# ## Attaching the Machine to a Busbar
#
# In order to use the machine model, we need to attach it to a
# busbar, thus forming a system which satisfies the MTKBus Interface.
# There are two ways of doing so: manually and using the `MTKBus` constructor.
#
# **Manual Construction**
#
# We need to define a new MTK model, which has 2 components: a busbar and the machine.
# Both components have a `terminal` as a subcomponent, we can use the `connect` function
# to hook the machine on the busbar.

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

# **Automatic Construction**
#
# We can also use the `MTKBus` constructor to create a busbar with a machine attached.
# This constructor takes a list of *injector* models and hooks them all to the same busbar.

mtkbus = MTKBus(machine; name=:bus)
isbusmodel(mtkbus) # assert that the created model satisfies the interface

# ## Compiling bus to `VertexModel`
# To actually simulate the system, we need to *compile* the model, i.e. transforming it from a
# purely symbolic representation to a numerical one.

Bus(mtkbus)

# ## Defining a Simulation Scenario
#
# To simulate the model, we need to define some kind of scenario. We'll create a simple
# two-bus system where our custom Milano machine is connected to a slack bus through
# a transmission line. This will allow us to observe the machine's dynamic behavior
# in response to a frequency disturbance.
#
# First, we create a slack bus that provides the voltage and frequency reference for the system.

slackbus = Bus(
    PowerDynamics.VariableFrequencySlack(name=:variable_slack),
    vidx=1,
    pf=pfSlack(V=1)
)

# We define a frequency event that increases the system frequency at t=1 second
# (see ND docs on Callbacks for details).
# This disturbance will cause our machine to respond dynamically as it tries to
# maintain synchronism with the network.

freq_event = PresetTimeComponentCallback(
    1, # trigger at time 1
    ComponentAffect([],[:V,:ω]) do u, p, ctx
        p[:ω] = 1.01 # set frequency to 1.01 pu
    end
)
set_callback!(slackbus, freq_event)

# Next, we create the generator bus using our custom Milano machine model.
# We specify it as a PV bus for the power flow with 1 pu voltage and 1 pu active power.

genbus = Bus(
    mtkbus,
    vidx=2,
    pf=pfPV(V=1, P=1)
)

# We connect the two buses with a simple PI transmission line model.

line = Line(MTKLine(PiLine(; name=:piline)); src=1,dst=2)

# Now we can build the complete network with our two buses and the connecting line.

nw = Network([slackbus, genbus], line)

# Before running dynamic simulation, we initialize the system from power flow.
# This ensures that all dynamic states start from a steady-state condition.
#
# To do so, we use the function `initialize_from_pf!`, which does several steps:
# 1. Calculate the powerflow according to the powerflow models.
# 2. Initialize the "free" states and parameters of the dynamical components, such that
#    the system is in a steady state.
#
# More information on initialization can be found in the docs on Powergrid Initialization.

initialize_from_pf!(nw)

# Let's examine the initial state of our generator bus to verify proper initialization.

dump_initial_state(nw[VIndex(2)])

# The printout shows us several important aspects:
# The free internal states $\delta$, $\omega$ and the free internal parameters
# $P_{\mathrm m}$ and $vf_{\mathrm{set}}$ have been initialized.
# We see, that both power and excitation voltage are slightly above the given (1,1) for the powerflow,
# which is expected since there are some losses in the model.
# However the initialized state matches the powerflow solution at the **network interface**,
# i.e. `busbar₊P` and `busbar₊u_mag` are both 1 pu.

# ## Dynamic Simulation
#
# With the system properly initialized, we can set up and solve the dynamic simulation.
# We simulate for 100 seconds to capture the machine's response to the frequency disturbance.

s0 = NWState(nw)
prob = ODEProblem(nw, uflat(s0), (0,100), pflat(s0), callback=get_callbacks(nw))
sol = solve(prob, Rodas5P())

# ## Visualizing the Results
#
# Now let's create comprehensive plots to visualize how our custom Milano machine
# responds to the frequency disturbance. We'll plot several key variables that
# demonstrate the machine's electromechanical dynamics.

let
    fig = Figure(size=(800, 600));

    ax1 = Axis(fig[1, 1];
        title="Rotor Angle",
        xlabel="Time [s]",
        ylabel="Rotor Angle δ [rad]")
    lines!(ax1, sol; idxs=VIndex(2, :machine₊δ), linewidth=2)
    axislegend(ax1)

    ax2 = Axis(fig[2, 1];
        title="Rotor Speed",
        xlabel="Time [s]",
        ylabel="Rotor Speed ω [pu]")
    lines!(ax2, sol; idxs=VIndex(2, :machine₊ω), linewidth=2)
    axislegend(ax2)

    ax3 = Axis(fig[3, 1];
        title="Machine Voltages",
        xlabel="Time [s]",
        ylabel="Voltage [pu]")
    lines!(ax3, sol; idxs=VIndex(2, :machine₊V_d), color=Cycled(1), linewidth=2)
    lines!(ax3, sol; idxs=VIndex(2, :machine₊V_q), color=Cycled(2), linewidth=2)
    axislegend(ax3)
    fig
end

# ## Observing the Poor Damping Problem
#
# From the plots above, we can see that the Milano classical machine exhibits very
# lightly damped oscillations that persist for a very long time. The rotor angle
# and speed oscillate for hundreds of seconds without settling to a steady state.
#
# This poor damping behavior occurs because:
# 1. **No damper windings**: The model lacks electromagnetic damping mechanisms
# 2. **Constant field voltage**: No dynamic response to help stabilize the machine
# 3. **No mechanical damping**: The swing equation has no friction losses
#
# The only source of damping here is, that we have specified a *constant
# mechanical power* rather than a constant mechanical torque.
#
# To solve this problem, real power systems use control systems, particularly
# **Power System Stabilizers (PSS)** that are specifically designed to damp
# electromechanical oscillations.
#
# ## Adding a Power System Stabilizer (PSS)
#
# Let's create an improved machine model with controllable field voltage and
# add the simplest possible PSS to demonstrate the damping improvement.
#
# The implemented PSS is a simple device, which adjusts the excitation voltage
# based on frequency deviation. It consists of a washout filter to remove steady-state errors
# and only react to frequency changes, and a gain to amplify the response.
#
# To achieve this goal we will:
# 1. Modify the Milano machine model to include a controllable field voltage input and a rotor frequency measurement output.
# 2. Create a simple PSS model that takes the frequency input and outputs a stabilizing signal to the field voltage.
# 3. Combine the machine and PSS into a new composite model that forms an injector.
# 4. Repeat the simulation above with our new controlled-generator model and compare the results.

# ### Controllable Machine Model
#
# First, we create a modified Milano machine with control inputs/outputs:
# `vf_in` for field voltage and `ω_out` for frequency output.

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
        # Park's transformations
        [terminal.u_r, terminal.u_i] .~ T_to_glob(δ)*[V_d, V_q]
        [I_d, I_q] .~ T_to_loc(δ)*[terminal.i_r, terminal.i_i]

        # mechanical swing equation Milano 15.5
        Dt(δ) ~ ω_b*(ω - 1)
        2*H * Dt(ω) ~ P_m/ω - τ_e

        # static flux linkage equations Milano 15.11
        ψ_d ~  V_q + R_s*I_q
        ψ_q ~ -V_d - R_s*I_d

        # electrical torque Milano 15.6
        τ_e ~ ψ_d*I_q - ψ_q*I_d

        # magnetic equations from static model Milano 15.36
        0 ~ V_q + R_s*I_q + X′_d*I_d - vf_in.u  # Use controllable input
        0 ~ V_d + R_s*I_d - X′_d*I_q

        # Control interface - output frequency for PSS
        ω_out.u ~ ω
    end
end

# ### Simple Power System Stabilizer
#
# The simplest PSS consists of a washout filter with gain. The washout filter
# ensures the PSS only responds to frequency changes, not steady-state errors.

@mtkmodel SimplePSS begin
    @components begin
        ω_in = RealInput() # frequency input from machine
        vst = RealOutput() # stabilizer output signal
    end
    @parameters begin
        Tw=10, [description="washout time constant"]
        Ks=20, [description="stabilizer gain"]
    end
    @variables begin
        y(t), [guess=0, description="washout filter output"]
    end
    @equations begin
        # Washout filter: dy/dt = (ω - y)/Tw
        Dt(y) ~ (ω_in.u - y) / Tw
        # output gain
        vst.u ~ Ks * (ω_in.u - y)
    end
end

# ### Complete Generator with PSS
#
# The PSS only adds an **offset** to the field voltage based on the frequency input.
# Therefore, our combined injector model needs to look something like this:
# ```
#     ┌───────────────────────────┐
#     │GeneratorWithPss           │
#     │         ╭─────→─────╮     │
# (t) │ ┌───────┴─┐ ω_out ┌─┴───┐ │
#  o──┼─┤ Machine │       │ PSS │ │
#     │ └───────┬─┘       └─┬───┘ │
#     │   vf_in ╰──←─(+)──←─╯ vst │
#     │               ↑           │
#     │            vf_base        │
#     └───────────────────────────┘
# ```
# Notably, similar to how we left `vf_set` free for initialization in the previous example,
# now we need to leave `vf_base` free.
#
# We define a new mtkmodel which combines machine with controller and forms a new injector:

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

# Since this is an injector, we can use `MTKBus(gen_with_pss)` to build the symbolic bus model.
# However, this leads to another level of namespacing, as the overall bus will have variable names like
# `gen_with_pss₊machine₊δ` due to the encapsulation.
#
# Alternatively, we could define a model which directly implements the `MTKBus` interface:
# ```
# ┌─────────────────────────────────────┐
# │MyMTKBus                             │
# │                   ╭─────→─────╮     │
# │┌──────┐   ┌───────┴─┐ ω_out ┌─┴───┐ │
# ││BusBar├─o─┤ Machine │       │ PSS │ │
# │└──────┘   └───────┬─┘       └─┬───┘ │
# │             vf_in ╰──←─(+)──←─╯ vst │
# │                         ↑           │
# │                      vf_base        │
# └─────────────────────────────────────┘
# ```

@mtkmodel CustomMTKBus begin
    @components begin
        busbar = BusBar()
        machine = MilanoControllableMachine()
        pss = SimplePSS()
    end
    @parameters begin
        vf_base, [guess=1.0, description="base field voltage"]
    end
    @equations begin
        connect(busbar.terminal, machine.terminal)
        connect(machine.ω_out, pss.ω_in)
        machine.vf_in.u ~ vf_base + pss.vst.u
    end
end
@named genbus_custom = CustomMTKBus()
@assert isbusmodel(genbus_custom)

# In practice, it doesn't really matter which approach you choose, as both will work.
# However this highlights the flexibility of the MTK modeling framework **before**
# you go to the compiled-model domain by calling `Bus` on the model fulfilling the `MTKBus` interface.

# ## Simulation with PSS
#
# Now let's run the same simulation scenario with the PSS-equipped generator
# to observe the damping improvement.

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

# ## Comparing Results: With and Without PSS
#
# Let's create comparison plots to clearly see the damping improvement:

let
    fig = Figure(size=(800, 600));

    # Compare rotor speeds
    ax1 = Axis(fig[1, 1];
        title="Rotor Speed Comparison: Effect of PSS on Damping",
        xlabel="Time [s]",
        ylabel="Rotor Speed ω [pu]")
    lines!(ax1, sol; idxs=VIndex(2, :machine₊ω), label="No PSS", color=Cycled(2))
    lines!(ax1, sol_pss; idxs=VIndex(2, :gen_with_pss₊machine₊ω), label="Simple PSS", color=Cycled(1), linewidth=2)
    axislegend(ax1, position=:rt)
    xlims!(ax1, 0, 30)  # Focus on first 30 seconds

    # PSS Output - shows the actual stabilizer signal
    ax2 = Axis(fig[2, 1];
        title="PSS Output Signal",
        xlabel="Time [s]",
        ylabel="PSS Output [pu]")
    lines!(ax2, sol_pss; idxs=VIndex(2, :gen_with_pss₊pss₊vst₊u), label="PSS Output", linewidth=2)
    axislegend(ax2, position=:rt)
    xlims!(ax2, 0, 30)

    fig
end

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
