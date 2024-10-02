using OpPoDyn
using OpPoDyn.Library
using ModelingToolkit
using NetworkDynamics
using Graphs
using OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqNonlinearSolve
using DiffEqCallbacks
using CairoMakie

@mtkmodel LoadBus begin
    @components begin
        busbar = BusBar()
        load = ConstantYLoad(Pset, Qset, Vset=nothing)
    end
    @equations begin
        connect(load.terminal, busbar.terminal)
    end
end

@mtkmodel ClassicBus begin
    @components begin
        machine = Library.ClassicalMachine(;
            τ_m_input=false,
            S_b=100,
            V_b=18,
            ω_b=2π*60,
            X′_d,
            R_s=0.0026,
            vf_set=nothing,
            τ_m_set=nothing,
            H,
        )
        busbar = BusBar()
    end
    @equations begin
        connect(machine.terminal, busbar.terminal)
    end
end

# generate all MTK bus models
@named mtkbus1 = ClassicBus(; machine__H=23.64, machine__X′_d=0.0608)
@named mtkbus2 = ClassicBus(; machine__H= 6.40, machine__X′_d=0.1198)
@named mtkbus3 = ClassicBus(; machine__H= 3.01, machine__X′_d=0.1813)
@named mtkbus4 = MTKBus()
@named mtkbus5 = LoadBus(;load__Pset=-1.25, load__Qset=-0.5)
@named mtkbus6 = LoadBus(;load__Pset=-0.90, load__Qset=-0.3)
@named mtkbus7 = MTKBus()
@named mtkbus8 = LoadBus(;load__Pset=-1.0, load__Qset=-0.35)
@named mtkbus9 = MTKBus()


# generate the dynamic component functions
@named bus1 = Bus(mtkbus1; vidx=1)
@named bus2 = Bus(mtkbus2; vidx=2)
@named bus3 = Bus(mtkbus3; vidx=3)
@named bus4 = Bus(mtkbus4; vidx=4)
@named bus5 = Bus(mtkbus5; vidx=5)
@named bus6 = Bus(mtkbus6; vidx=6)
@named bus7 = Bus(mtkbus7; vidx=7)
@named bus8 = Bus(mtkbus8; vidx=8)
@named bus9 = Bus(mtkbus9; vidx=9)

# Branches
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

# set powerflow for initialization
set_voltage!(bus1; mag=1.040, arg=deg2rad( 0.0))
set_voltage!(bus2; mag=1.025, arg=deg2rad( 9.3))
set_voltage!(bus3; mag=1.025, arg=deg2rad( 4.7))
set_voltage!(bus4; mag=1.026, arg=deg2rad(-2.2))
set_voltage!(bus5; mag=0.996, arg=deg2rad(-4.0))
set_voltage!(bus6; mag=1.013, arg=deg2rad(-3.7))
set_voltage!(bus7; mag=1.026, arg=deg2rad( 3.7))
set_voltage!(bus8; mag=1.016, arg=deg2rad( 0.7))
set_voltage!(bus9; mag=1.032, arg=deg2rad( 2.0))

# set currents (via PQ) for buses which need initialization
set_current!(bus1; P=0.716, Q= 0.270)
set_current!(bus2; P=1.630, Q= 0.067)
set_current!(bus3; P=0.850, Q=-0.109)
set_current!(bus5; P=-1.25, Q=-0.5)
set_current!(bus6; P=-0.90, Q=-0.3)
set_current!(bus8; P=-1.00, Q=-0.35)

# initialize generators
initialize_component!(bus1)
initialize_component!(bus2)
initialize_component!(bus3)

# we also initialize the loads to get the correct nominal voltage
initialize_component!(bus5)
initialize_component!(bus6)
initialize_component!(bus8)

# build network
vertexfs = [bus1, bus2, bus3, bus4, bus5, bus6, bus7, bus8, bus9];
edgefs = [l45, l46, l57, l69, l78, l89, t14, t27, t39];
nw = Network(vertexfs, edgefs)
u0 = NWState(nw)

# create fault
affect! = (integrator) -> begin
    if integrator.t == 0.0833
        @info "Deactivate line 6 at t = $(integrator.t)"
        p = NWParameter(integrator)
        p.e[6, :pibranch₊active] = 0
    else
        error("Schould not be reached.")
    end
end
cb = PresetTimeCallback([0.0833], affect!)

prob = ODEProblem(nw, uflat(u0), (0,2), pflat(u0); callback=cb)
sol = solve(prob, Rodas5P());
nothing

break # stop execution of script here

fig = Figure();
ax = Axis(fig[1, 1]; title="Läuferwinkel")
ts = range(0,2,length=1000)
δ2 = sol(ts, idxs=VIndex(2, :machine₊δ)).u - sol(ts, idxs=VIndex(1, :machine₊δ)).u
δ3 = sol(ts, idxs=VIndex(3, :machine₊δ)).u - sol(ts, idxs=VIndex(1, :machine₊δ)).u
lines!(ax, ts, rad2deg.(δ2))
lines!(ax, ts, rad2deg.(δ3))
fig

# get first and maximum value
rad2deg(first(δ2)) => rad2deg(maximum(δ2))
rad2deg(first(δ3)) => rad2deg(maximum(δ3))


fig = Figure();
ax = Axis(fig[1, 1]; title="Bus voltage magnitude")
ts = range(0,2,length=1000)
umag5 = sqrt.(sol(ts; idxs=VIndex(5, :busbar₊u_r)).^2 + sol(ts; idxs=VIndex(5, :busbar₊u_i)).^2)
umag7 = sqrt.(sol(ts; idxs=VIndex(7, :busbar₊u_r)).^2 + sol(ts; idxs=VIndex(7, :busbar₊u_i)).^2)
lines!(ax, ts, umag5.u; label="Bus5")
lines!(ax, ts, umag7.u; label="Bus7")
axislegend(ax)
fig

i = findfirst(>(0.0833), ts)
umag5[1] => umag5[i]
umag7[1] => umag7[i]


# Plotting the Solution
fig = Figure(size=(1000,2000));
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
ax = Axis(fig[3, 1]; title="Voltag angel")
for i in 1:9
    lines!(ax, sol; idxs=VIndex(i,:busbar₊u_arg), label="Bus $i", color=Cycled(i))
end
axislegend(ax)
ax = Axis(fig[4, 1]; title="Frequency")
for i in 1:3
    lines!(ax, sol; idxs=VIndex(i,:machine₊ω), label="Bus $i", color=Cycled(i))
end
axislegend(ax)
fig
