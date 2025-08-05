#=
# [IEEE 9Bus Example](@id ieee9bus)

In this example, we're going to model the IEEE 9 bus system.

The parameters are mainly adopted from the RTDS data.
=#
using PowerDynamics
using PowerDynamics.Library
using ModelingToolkit
using NetworkDynamics
using Graphs
using OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqNonlinearSolve
using CairoMakie

#=
## Generator Busses

The 3 generator buses are modeld using a SauerPai 6th order machine model with
variable field voltage and mechanical torque input.
The field voltage is provided by an `AVRTypeI`, the torque is provide by a `TGOV1` model.
=#
function GeneratorBus(; machine_p=(;), avr_p=(;), gov_p=(;))
    @named machine = SauerPaiMachine(;
        vf_input=true,
        τ_m_input=true,
        S_b=100,
        V_b=1,
        ω_b=2π*60,
        R_s=0.000125,
        T″_d0=0.01,
        T″_q0=0.01,
        machine_p... # unpack machine parameters
    )
    @named avr = AVRTypeI(; vr_min=-5, vr_max=5,
        Ka=20, Ta=0.2,
        Kf=0.063, Tf=0.35,
        Ke=1, Te=0.314,
        E1=3.3, Se1=0.6602, E2=4.5, Se2=4.2662,
        tmeas_lag=false,
        avr_p... # unpack AVR parameters
    )
    @named gov = TGOV1(; R=0.05, T1=0.05, T2=2.1, T3=7.0, DT=0, V_max=5, V_min=-5,
        gov_p... # unpack governor parameters
    )
    ## generate the "injector" as combination of multiple components
    injector = CompositeInjector([machine, avr, gov]; name=:generator)

    ## generate the MTKBus (i.e. the MTK model containg the busbar and the injector)
    mtkbus = MTKBus(injector)
end
nothing # hide


#=
## Load Busses

The dynamic loads are modeld as static Y-loads. Those have 3 parameters: `Pset`, `Qset` and `Vset`.
For now, those parameters will be left free. We'll initialize them later on from the powerflow results.

The `Vset` parameter is left free for now. Later on it is automaticially determined to match the
behavior of the static power flow load model.
=#
function ConstantYLoadBus()
    @named load = ConstantYLoad()
    MTKBus(load; name=:loadbus)
end
nothing #hide


#=
## Generate Dynamical models

The parameters of the machines are obtaind from the data table from the RTDS datasheet.
=#
gen1p = (;X_ls=0.01460, X_d=0.1460, X′_d=0.0608, X″_d=0.06, X_q=0.1000, X′_q=0.0969, X″_q=0.06, T′_d0=8.96, T′_q0=0.310, H=23.64)
gen2p = (;X_ls=0.08958, X_d=0.8958, X′_d=0.1198, X″_d=0.11, X_q=0.8645, X′_q=0.1969, X″_q=0.11, T′_d0=6.00, T′_q0=0.535, H= 6.40)
gen3p = (;X_ls=0.13125, X_d=1.3125, X′_d=0.1813, X″_d=0.18, X_q=1.2578, X′_q=0.2500, X″_q=0.18, T′_d0=5.89, T′_q0=0.600, H= 3.01)
nothing #hide

#=
We instantiate all models as modeling toolkit models.
=#
mtkbus1 = GeneratorBus(; machine_p=gen1p)
mtkbus2 = GeneratorBus(; machine_p=gen2p)
mtkbus3 = GeneratorBus(; machine_p=gen3p)
mtkbus4 = MTKBus() # <- bus with no injectors, essentially
mtkbus5 = ConstantYLoadBus()
mtkbus6 = ConstantYLoadBus()
mtkbus7 = MTKBus()
mtkbus8 = ConstantYLoadBus()
mtkbus9 = MTKBus()
nothing #hide

#=
After this, we can build the `NetworkDynamics` components using the `Bus`-constructor.

The `Bus` constructor is essentially a thin wrapper around the `VertexModel` constructor which,
per default, adds some metadata. For example the `vidx` property which later on allows for
"graph free" network dynamics instantiation.

We use the `pf` keyword to specify the models which should be used in the powerflow calculation.
Here, generator 1 is modeld as a slack bus while the other two generators are modeled as a PV bus.
The loads are modeled as PQ buses.
=#
@named bus1 = Bus(mtkbus1; vidx=1, pf=pfSlack(V=1.04))
@named bus2 = Bus(mtkbus2; vidx=2, pf=pfPV(V=1.025, P=1.63))
@named bus3 = Bus(mtkbus3; vidx=3, pf=pfPV(V=1.025, P=0.85))
@named bus4 = Bus(mtkbus4; vidx=4)
@named bus5 = Bus(mtkbus5; vidx=5, pf=pfPQ(P=-1.25, Q=-0.5))
@named bus6 = Bus(mtkbus6; vidx=6, pf=pfPQ(P=-0.9, Q=-0.3))
@named bus7 = Bus(mtkbus7; vidx=7)
@named bus8 = Bus(mtkbus8; vidx=8, pf=pfPQ(P=-1.0, Q=-0.35))
@named bus9 = Bus(mtkbus9; vidx=9)
nothing #hide

#=
Later on in initialization, we want to ensure that the internal parameters
of the loads are initialized to match the powerflow results.
We have 3 free parameters: `Vset`, `Pset` and `Qset`.
To help the initialization, we can provide a [`InitFormula`](@extref NetworkDynamics.InitFormula)
to set the `Vset` parameter to the voltage magnitude of the bus.
=#
vset_formula = @initformula :load₊Vset = sqrt(:busbar₊u_r^2 + :busbar₊u_i^2)
add_initformula!(bus5, vset_formula)
add_initformula!(bus6, vset_formula)
add_initformula!(bus8, vset_formula)
nothing #hide

#=
## Branches

Branches and Transformers are build from the same `PILine` model with optional
transformer on both ends. However, the data is provided in a way that the actual
transformer values are 1.0. Apparently, the transforming action has been absorbed
into the line parameters according to the base voltage on both ends.

For the lines we again make use of the `src` and `dst` metadata of the
`EdgeModel` objects for automatic graph construction.
=#
function piline(; R, X, B)
    @named pibranch = PiLine(;R, X, B_src=B/2, B_dst=B/2, G_src=0, G_dst=0)
    MTKLine(pibranch)
end
function transformer(; R, X)
    @named transformer = PiLine(;R, X, B_src=0, B_dst=0, G_src=0, G_dst=0)
    MTKLine(transformer)
end

@named l45 = Line(piline(; R=0.0100, X=0.0850, B=0.1760), src=4, dst=5)
@named l46 = Line(piline(; R=0.0170, X=0.0920, B=0.1580), src=4, dst=6)
@named l57 = Line(piline(; R=0.0320, X=0.1610, B=0.3060), src=5, dst=7)
@named l69 = Line(piline(; R=0.0390, X=0.1700, B=0.3580), src=6, dst=9)
@named l78 = Line(piline(; R=0.0085, X=0.0720, B=0.1490), src=7, dst=8)
@named l89 = Line(piline(; R=0.0119, X=0.1008, B=0.2090), src=8, dst=9)
@named t14 = Line(transformer(; R=0, X=0.0576), src=1, dst=4)
@named t27 = Line(transformer(; R=0, X=0.0625), src=2, dst=7)
@named t39 = Line(transformer(; R=0, X=0.0586), src=3, dst=9)
nothing #hide

#=
## Build Network

Finally, we can build the network by providing the vertices and edges.
=#
vertexfs = [bus1, bus2, bus3, bus4, bus5, bus6, bus7, bus8, bus9];
edgefs = [l45, l46, l57, l69, l78, l89, t14, t27, t39];
nw = Network(vertexfs, edgefs; warn_order=false)

#=
## System Initialization

To initialize the system for dynamic simulation, we use `initialize_from_pf!` which 
performs a unified powerflow solving and component initialization process.

Internally, this function:
1. Builds an equivalent static powerflow network from the dynamic models
2. Solves the static powerflow equations using the specified powerflow models (pfSlack, pfPV, pfPQ)
3. Uses the powerflow solution to initialize all dynamic component states and parameters

This ensures that the dynamic model starts from a steady-state condition that matches
the powerflow solution. Specifically, it determines:
- Bus voltages and currents from the powerflow solution
- Unknown `Vset` parameters for load buses
- Internal machine states (flux linkages, rotor angles, etc.)
- Controller states and references for AVRs and governors
=#
u0 = initialize_from_pf(nw);
nothing #hide

#=
We could check the initial state of some of the variables, we expect the model to be initialized in
a way, that the setpoint of the constant Y-loads matches the powerfow constraints.
=#
u0[VIndex(5, :load₊Pset)] ≈ -1.25 && u0[VIndex(5, :load₊Qset)] ≈ -0.5

#=
## Disturbance

To see some dynamics, we need to introduce some disturbance.
For that we use a [`PresetTimeComponentCallback`](@extref NetworkDynamics.PresetTimeComponentCallback)
to deactivate a line at a certain time.
=#
deactivate_line = ComponentAffect([], [:pibranch₊active]) do u, p, ctx
    @info "Deactivate line $(ctx.src)=>$(ctx.dst) at t=$(ctx.t)"
    p[:pibranch₊active] = 0
end
cb = PresetTimeComponentCallback([1.0], deactivate_line)
set_callback!(l46, cb)
l46 # printout shows that the callback is set

#=
## Dynamic Simulation

With the system properly initialized, we can now set up and run the dynamic simulation.
We create an ODE problem using the initialized state and simulate the system response
to the line outage disturbance.
=#
prob = ODEProblem(nw, uflat(u0), (0,15), pflat(u0); callback=get_callbacks(nw))
sol = solve(prob, Rodas5P())
nothing #hide

#=
## Plotting the Solution

Finally, we visualize the simulation results showing the system response to the 
line outage at t=1.0 seconds. The plots show active power, voltage magnitudes,
and generator frequencies across the simulation time.
=#
fig = Figure(size=(600,800));

## Active power at selected buses
ax = Axis(fig[1, 1]; title="Active Power", xlabel="Time [s]", ylabel="Power [pu]")
for i in [1,2,3,5,6,8]
    lines!(ax, sol; idxs=VIndex(i,:busbar₊P), label="Bus $i", color=Cycled(i))
end
axislegend(ax)

## Voltage magnitude at all buses
ax = Axis(fig[2, 1]; title="Voltage Magnitude", xlabel="Time [s]", ylabel="Voltage [pu]")
for i in 1:9
    lines!(ax, sol; idxs=VIndex(i,:busbar₊u_mag), label="Bus $i", color=Cycled(i))
end

## Generator frequencies
ax = Axis(fig[3, 1]; title="Generator Frequency", xlabel="Time [s]", ylabel="Frequency [pu]")
for i in 1:3
    lines!(ax, sol; idxs=VIndex(i,:generator₊machine₊ω), label="Gen $i", color=Cycled(i))
end
axislegend(ax)

fig #hide
