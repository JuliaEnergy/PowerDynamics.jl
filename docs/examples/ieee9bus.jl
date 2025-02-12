#=
# IEEE 9Bus Example

In this example, we're going to model the IEEE 9 bus system.

The parameters are mainly adopted from the RTDS data.
=#
using OpPoDyn
using OpPoDyn.Library
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
@mtkmodel GenBus begin
    @components begin
        machine = SauerPaiMachine(;
            vf_input=true,
            τ_m_input=true,
            S_b=100,
            V_b=1,
            ω_b=2π*60,
            R_s=0.000125,
            T″_d0=0.01,
            T″_q0=0.01,
            X_d, X′_d, X″_d, X_q, X′_q, X″_q, X_ls, T′_d0, T′_q0, H # free per machine parameter
        )
        avr = AVRTypeI(vr_min=-5, vr_max=5,
            Ka=20, Ta=0.2,
            Kf=0.063, Tf=0.35,
            Ke=1, Te=0.314,
            E1=3.3, Se1=0.6602, E2=4.5, Se2=4.2662,
            tmeas_lag=false)
        gov = TGOV1(R=0.05, T1=0.05, T2=2.1, T3=7.0, DT=0, V_max=5, V_min=-5)
        busbar = BusBar()
    end
    @equations begin
        connect(machine.terminal, busbar.terminal)
        connect(machine.v_mag_out, avr.vh)
        connect(avr.vf, machine.vf_in)
        connect(gov.τ_m, machine.τ_m_in)
        connect(machine.ωout, gov.ω_meas)
    end
end
nothing # hide


#=
## Load Busses

The dynamic loads are modeld as static Y-loads. Those have 3 parameters: `Pset`, `Qset` and `Vset`.
The `Vset` parameter is left free for now. Later on it is automaticially determined to match the
behavior of the static power flow load model.
=#
@mtkmodel LoadBus begin
    @components begin
        busbar = BusBar()
        load = ConstantYLoad(Pset, Qset)
    end
    @equations begin
        connect(load.terminal, busbar.terminal)
    end
end
nothing #hide


#=
## Generate Dynamical models

The parameters of the machines are obtaind from the data table from the RTDS datasheet.
=#
gen1p = (;machine__X_ls=0.01460, machine__X_d=0.1460, machine__X′_d=0.0608, machine__X″_d=0.06, machine__X_q=0.1000, machine__X′_q=0.0969, machine__X″_q=0.06, machine__T′_d0=8.96, machine__T′_q0=0.310, machine__H=23.64)
gen2p = (;machine__X_ls=0.08958, machine__X_d=0.8958, machine__X′_d=0.1198, machine__X″_d=0.11, machine__X_q=0.8645, machine__X′_q=0.1969, machine__X″_q=0.11, machine__T′_d0=6.00, machine__T′_q0=0.535, machine__H= 6.40)
gen3p = (;machine__X_ls=0.13125, machine__X_d=1.3125, machine__X′_d=0.1813, machine__X″_d=0.18, machine__X_q=1.2578, machine__X′_q=0.2500, machine__X″_q=0.18, machine__T′_d0=5.89, machine__T′_q0=0.600, machine__H= 3.01)
nothing #hide

#=
We instantiate all models as modeling toolkit models.
=#
@named mtkbus1 = GenBus(; gen1p...)
@named mtkbus2 = GenBus(; gen2p...)
@named mtkbus3 = GenBus(; gen3p...)
@named mtkbus4 = MTKBus()
@named mtkbus5 = LoadBus(;load__Pset=-1.25, load__Qset=-0.5)
@named mtkbus6 = LoadBus(;load__Pset=-0.90, load__Qset=-0.3)
@named mtkbus7 = MTKBus()
@named mtkbus8 = LoadBus(;load__Pset=-1.0, load__Qset=-0.35)
@named mtkbus9 = MTKBus()
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
nw = Network(vertexfs, edgefs)
break

#=
## Powerflow

To initialize the system, we first solve the static powerflow problem.
Internally, `OpPoDyn.jl` builds an equivalent network but replaces each dynamic model
with the given static power flow model.
The static powerflow problem is solved, after which the power flow solution (bus voltages
and currents) are stored as `default` values in the vertex models.
=#
solve_powerflow!(nw)
nothing #hide

#=
## Component initialization

The power flow solution provided all the "interface states" (i.e. voltages and currents).
With that information, we can initialize the free states and parameters of the dynamic models,
such that the dynamic steady state matches the static power flow solution.

When calling [`initialize!`](@ref),  `OpPoDyn.jl` will loop through all the dynamic models
in the system, automaticially creating and solving a nonlinear initialization problem for each of them.

Concretly, here were solving for the following things:
- unknown `Vset` for load busses,
- unknown internal machine and controller states as well as the free govenor and
  avr references (parameters) of the generator busses.
=#
initialize!(nw)
nothing #hide

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
noting # hide

#=
## Build Network

Finally, we can build the network by providing the vertices and edges.
=#
u0 = NWState(nw)
prob = ODEProblem(nw, uflat(u0), (0,15), pflat(u0); callback=get_callbacks(nw))
sol = solve(prob, Rodas5P())
nothing

#=
## Plotting the Solution
=#
fig = Figure(size=(600,800));
ax = Axis(fig[1, 1]; title="Active power")
for i in [1,2,3,5,6,8]
    lines!(ax, sol; idxs=VIndex(i,:busbar₊P), label="Bus $i", color=Cycled(i))
end
axislegend(ax)
ax = Axis(fig[2, 1]; title="Voltage magnitude")
for i in 1:9
    lines!(ax, sol; idxs=VIndex(i,:busbar₊u_mag), label="Bus $i", color=Cycled(i))
end
axislegend(ax)
ax = Axis(fig[3, 1]; title="Frequency")
for i in 1:3
    lines!(ax, sol; idxs=VIndex(i,:machine₊ω), label="Bus $i", color=Cycled(i))
end
axislegend(ax)
fig #hide
