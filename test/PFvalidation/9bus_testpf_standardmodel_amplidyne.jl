using OpPoDyn
using OpPoDyn.Library
using ModelingToolkit
using NetworkDynamics
using NetworkDynamicsInspector
using Graphs
using OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqNonlinearSolve
using DiffEqCallbacks
using CairoMakie
using CSV
using DataFrames

@mtkmodel LoadBus begin
    @components begin
        busbar = BusBar()
        load = ConstantYLoad(Pset, Qset, Vset=nothing)
    end
    @equations begin
        connect(load.terminal, busbar.terminal)
    end
end

@mtkmodel StandardBus begin
    @components begin
        machine = Library.StandardModel_pf_testneu(;
            S_b,
            V_b,
            Sn,
            Vn,
            ω_b=2π*60,
            vf_input=false,
            τ_m_input=false,
            R_s,
            X_rld,
            X_rlq,
            X″_d,
            X″_q,
            X_ls,
            X_ad,
            X_aq,
            X_1q,
            X_det_d,
            X_det_q,
            X_fd_loop,
            X_1d_loop,
            X_1q_loop,
            X_2q_loop,
            k_fd,
            k_1d,
            k_1q,
            k_2q,
            R_fd,
            R_1d,
            R_1q,
            R_2q,
            H,
            D,
            dpe,
            dkd,
            cosn,
            salientpole,
            dpu,
            addmt,
            xmdm
        )
        busbar = BusBar()
    end
    @equations begin
        connect(machine.terminal, busbar.terminal)
    end
end

# Branches
function piline(; R, X, B)
    @named pibranch = PiLine(;R, X, B_src=B/2, B_dst=B/2, G_src=0, G_dst=0)
    MTKLine(pibranch)
end

function piline_shortcircuit(; R, X, B, pos, G_fault=0, B_fault=0)
    #faultimp = if (G_fault + B_fault) == 0
       # 0
    #else
        #1
    #end
    @named pibranch = PiLine_fault(;R, X, B_src=B/2, B_dst=B/2, G_src=0, G_dst=0, G_fault, B_fault, pos)#, faultimp)
    MTKLine(pibranch)
end

function transformer(; R, X)
    @named transformer = PiLine(;R, X, B_src=0, B_dst=0, G_src=0, G_dst=0)
    MTKLine(transformer)
end


#calculate generator parameters externally
function secondary_from_primary(; ω_b, X_rld, X_rlq, X_d, X_q, X′_d, X′_q, X″_d, X″_q, X_ls, T′_d0, T″_d0, T′_q0, T″_q0, salientpole, kwargs...)
    #conversion of time constants (not exact conversion) (45)
    T″_d = T″_d0 * X″_d/X′_d
    T″_q = T″_q0 * X″_q/(X′_q * (1-salientpole) + salientpole * X_q)
    T′_d = T′_d0 * X′_d/X_d
    T′_q = T′_q0 * X′_q/X_q

    #(47)-(53): calculation of equivalent model parameters
    X_ad = X_d - X_ls
    X_aq = X_q - X_ls

    X_1 = X_d - X_ls + X_rld
    X_2 = X_1 - (X_d - X_ls)^2 / X_d
    X_3 = (X_2 - X_1 * X″_d/X_d) / (1 - X″_d/X_d)
    T_1 = X_d / X′_d * T′_d + (1 - X_d/X′_d + X_d/X″_d) * T″_d
    T_2 = T″_d + T′_d
    a = (X_2 * T_1 - X_1 * T_2) / (X_1 - X_2)
    b = X_3 * T″_d * T′_d / (X_3 - X_2)
    T_σfd = -a/2 + sqrt(a^2/4 - b)
    T_σ1d = -a/2 - sqrt(a^2/4 - b)
    X_fd = (T_σfd - T_σ1d) / ((T_1 - T_2)/(X_1 - X_2) + T_σ1d / X_3)
    X_1d = (T_σ1d - T_σfd) / ((T_1 - T_2)/(X_1 - X_2) + T_σfd / X_3)
    R_fd = X_fd / (ω_b * T_σfd)
    R_1d = X_1d / (ω_b * T_σ1d)

    #selbst berechnet ggf Fehler: The q-axis model parameters can be calculated analogously to the d-axis parameters in case of a round-rotor machine
    X_4 = X_q - X_ls + X_rlq
    X_5 = X_4 - (X_q - X_ls)^2 / X_q
    X_6 = (X_5 - X_4 * X″_q/X_q) / (1 - X″_q/X_q)
    T_3 = X_q / X′_q * T′_q + (1 - X_q/X′_q + X_q/X″_q) * T″_q
    T_4 = T″_q + T′_q
    c = (X_5 * T_3 - X_4 * T_4) / (X_4 - X_5)
    d = X_6 * T″_q * T′_q / (X_6 - X_5)
    T_σ2q = -c/2 + sqrt(c^2/4 - d)
    T_σ1q = -c/2 - sqrt(c^2/4 - d)
    X_2qr = (T_σ2q - T_σ1q) / ((T_3 - T_4)/(X_4 - X_5) + T_σ1q / X_6) #round-rotor
    X_1qr = (T_σ1q - T_σ2q) / ((T_3 - T_4)/(X_4 - X_5) + T_σ2q / X_6) #round-rotor
    R_2qr = X_2qr / (ω_b * T_σ2q) #round-rotor
    R_1qr = X_1qr / (ω_b * T_σ1q) #round-rotor
    X_1qs = (X_q - X_ls) * (X″_q - X_ls) / (X_q - X″_q) #salient pole
    R_1qs = X″_q / X_q * (X_q - X_ls + X_1qs) / (ω_b * T″_q) #salient pole
    X_1q = salientpole * X_1qs + (1-salientpole) * X_1qr
    R_1q = salientpole * R_1qs + (1-salientpole) * R_1qr
    X_2q = (1-salientpole) * X_2qr
    R_2q = (1-salientpole) * R_2qr

    #(63)-(65)
    k_fd = (X_ad * X_1d) / ((X_ad + X_rld) * (X_1d + X_fd) + X_fd * X_1d)
    k_1d = (X_ad * X_fd) / ((X_ad + X_rld) * (X_1d + X_fd) + X_fd * X_1d)
    # XXX: X″_d wird hier definiert aber schon oben benutzt?
    X″_d = X_ad + X_ls - (k_1d + k_fd) * X_ad #??
    k_1qs = X_aq / (X_aq + X_rlq + X_1q) #salient pole
        #k_2qs = 0 #salient pole
    X″_qs = X_aq + X_ls - k_1qs * X_aq #salient pole
    k_1qr = (X_aq * X_2q)/((X_aq + X_rlq) * (X_2q + X_1q) + X_2q * X_1q) #round rotor
    k_2qr = (X_aq * X_1q)/((X_aq + X_rlq) * (X_2q + X_1q) + X_2q * X_1q) #round rotor
    X″_qr = X_aq + X_ls - (k_2qr + k_1qr) * X_aq #round rotor
    k_1q = salientpole * k_1qs + (1-salientpole) * k_1qr
    k_2q = (1-salientpole) * k_2qr
    # XXX: X″_q wird hier definiert aber schon oben benutzt?
    X″_q = salientpole * X″_qs + (1-salientpole)* X″_qr

    #(69), (71)
    X_det_d = (X_ad + X_rld) * (X_1d + X_fd) + X_fd * X_1d
    X_det_q = (X_aq + X_rlq) * (X_2q + X_1q) + X_2q * X_1q
    X_fd_loop = X_ad + X_rld + X_fd
    X_1d_loop = X_ad + X_rld + X_1d
    X_1q_loop = X_aq + X_rlq + X_1q
    X_2q_loop = X_aq + X_rlq + X_2q
    return (;X″_d, X″_q, X_ad, X_aq, X_det_d, X_det_q, X_fd_loop, X_1d_loop, X_1q_loop, X_2q_loop, k_fd, k_1d, k_1q, k_2q, X_1q, R_1q, X_2q,
            R_2q, X_fd, R_1d, R_fd, X_1d)
end


# generate all MTK bus models
#aus PF Implementierung -> S_b, Sn in MVA, V_b, Vn in kV
primary_parameters_gen1 = (;
    S_b=100e6,
    Sn=247.5e6,
    V_b=16.5,
    Vn=16.5,
    ω_b=2π*60,
    H=9.551516,
    D=0,
    R_s=0,
    X_rld=0,
    X_rlq=0,
    X_d=0.36135,
    X_q=0.2398275,
    X′_d=0.15048,
    X′_q=0.0001,
    X″_d=0.1,
    X″_q=0.1,
    X_ls=0.08316,
    T′_d0=8.96,
    T″_d0=0.075,
    T′_q0=0.0001,
    T″_q0=0.15,
    cosn=1,
    dkd=0,
    dpe=0,
    salientpole=1,
    dpu=0,
    addmt=0,
    xmdm=0,
)
secondary_parameters_gen1 = secondary_from_primary(; primary_parameters_gen1...)
allp_gen1 = Dict(pairs(primary_parameters_gen1)..., pairs(secondary_parameters_gen1)...)
# get rid of parametes which are not needed
for p in [:X_d, :X′_q, :T′_d0, :T′_q0, :X′_d, :T″_q0, :T″_d0, :X_q, :X_fd, :X_1d, :X_2q]
    delete!(allp_gen1, p)
end
renamedp_gen1 = Dict(map(s->Symbol("machine__", s), collect(keys(allp_gen1))) .=> Base.values(allp_gen1))


primary_parameters_gen2 = (;
    S_b=100e6,
    Sn=192e6,
    V_b=18,
    Vn=18,
    ω_b=2π*60,
    H=3.921568,
    D=0,
    R_s=0.005,
    X_rld=0,
    X_rlq=0,
    X_d=1.719936,
    X_q=1.65984,
    X′_d=0.230016,
    X′_q=0.378048,
    X″_d=0.2,
    X″_q=0.2,
    X_ls=0.100032,
    T′_d0=6,
    T″_d0=0.0575,
    T′_q0=0.535,
    T″_q0=0.0945,
    cosn=0.85,
    dkd=0,
    dpe=0,
    salientpole=0,
    dpu=0,
    addmt=0,
    xmdm=0,
)
secondary_parameters_gen2 = secondary_from_primary(; primary_parameters_gen2...)
allp_gen2 = Dict(pairs(primary_parameters_gen2)..., pairs(secondary_parameters_gen2)...)
# get rid of parametes which are not needed
for p in [:X_d, :X′_q, :T′_d0, :T′_q0, :X′_d, :T″_q0, :T″_d0, :X_q, :X_fd, :X_1d, :X_2q]
    delete!(allp_gen2, p)
end
#renamedp_gen2 = Dict(map(s->Symbol("machine__", s), collect(keys(allp_gen2))) .=> Base.values(allp_gen2))

primary_parameters_gen3 = (;
    S_b=100e6,
    Sn=128e6,
    V_b=13.8,
    Vn=13.8,
    ω_b=2π*60,
    H=2.766544,
    D=0,
    R_s=0.0001,
    X_rld=0,
    X_rlq=0,
    X_d=1.68,
    X_q=1.609984,
    X′_d=0.232064,
    X′_q=0.32,
    X″_d=0.2,
    X″_q=0.2,
    X_ls=0.094976,
    T′_d0=5.89,
    T″_d0=0.0575,
    T′_q0=0.6,
    T″_q0=0.08,
    cosn=0.85,
    dkd=0,
    dpe=0,
    salientpole=0,
    dpu=0,
    addmt=0,
    xmdm=0,
)
secondary_parameters_gen3 = secondary_from_primary(; primary_parameters_gen3...)
allp_gen3 = Dict(pairs(primary_parameters_gen3)..., pairs(secondary_parameters_gen3)...)
# get rid of parametes which are not needed
for p in [:X_d, :X′_q, :T′_d0, :T′_q0, :X′_d, :T″_q0, :T″_d0, :X_q, :X_fd, :X_1d, :X_2q]
    delete!(allp_gen3, p)
end
renamedp_gen3 = Dict(map(s->Symbol("machine__", s), collect(keys(allp_gen3))) .=> Base.values(allp_gen3))


@named mtkbus1 = StandardBus(; renamedp_gen1...)
#@named mtkbus2 = StandardBus(; renamedp_gen2...)
@named mtkbus3 = StandardBus(; renamedp_gen3...)
@named mtkbus4 = MTKBus()
@named mtkbus5 = LoadBus(;load__Pset=-1.25, load__Qset=-0.5)
@named mtkbus6 = LoadBus(;load__Pset=-0.90, load__Qset=-0.3)
@named mtkbus7 = MTKBus()
@named mtkbus8 = LoadBus(;load__Pset=-1.0, load__Qset=-0.35)
@named mtkbus9 = MTKBus()

components = [] # vector to collect dynamical components
controleqs = Equation[] # equations connecting avr and gov to machine
@named machine = Library.StandardModel_pf_testneu(; allp_gen2...)
@named avr = AVRTypeIS(
            # vref = p.Vref # let this free for initialization
            Ka = 25, Ke = -0.044, Kf = 0.0805,
            Ta = 0.2, Tf = 0.35, Te = 0.5, Tr = 0.06,
            vr_min = -1, vr_max = 1,
            A=0.0016, B=1.465,
        )
append!(controleqs, [connect(machine.v_mag_out, avr.vh), connect(avr.vf, machine.vf_in)])
comp = CompositeInjector(
        [machine, avr],
        controleqs,
        name=:ctrld_gen
    )
push!(components, comp)

@named bus2 = Bus(MTKBus(components...); vidx=2, pf=pfPV(V=1.025, P=1.63))


# generate the dynamic component functions
@named bus1 = Bus(mtkbus1; vidx=1, pf=pfSlack(V=1.04))
#@named bus2 = Bus(mtkbus2; vidx=2, pf=pfPV(V=1.025, P=1.63))
@named bus3 = Bus(mtkbus3; vidx=3, pf=pfPV(V=1.025, P=0.85))
@named bus4 = Bus(mtkbus4; vidx=4)
@named bus5 = Bus(mtkbus5; vidx=5, pf=pfPQ(P=-1.25, Q=-0.5))
@named bus6 = Bus(mtkbus6; vidx=6, pf=pfPQ(P=-0.9, Q=-0.3))
@named bus7 = Bus(mtkbus7; vidx=7)
@named bus8 = Bus(mtkbus8; vidx=8, pf=pfPQ(P=-1.0, Q=-0.35))
@named bus9 = Bus(mtkbus9; vidx=9)
#@named l45 = Line(piline(; R=0.0100, X=0.0850, B=0.1760), src=4, dst=5)

@named l46 = Line(piline(; R=0.0170, X=0.0920, B=0.1580), src=4, dst=6)
@named l69 = Line(piline(; R=0.0390, X=0.1700, B=0.3580), src=6, dst=9)
@named l78 = Line(piline(; R=0.0085, X=0.0720, B=0.1490), src=7, dst=8)
@named l89 = Line(piline(; R=0.0119, X=0.1008, B=0.2090), src=8, dst=9)
@named t14 = Line(transformer(; R=0, X=0.0576), src=1, dst=4)
@named t27 = Line(transformer(; R=0, X=0.0625), src=2, dst=7)
@named t39 = Line(transformer(; R=0, X=0.0586), src=3, dst=9)
@named l57 = Line(piline_shortcircuit(; R=0.0320, X=0.1610, B=0.3060, pos=0.99), src=5, dst=7) #S_b = 100 MVA, U_b = 230 kV; 2 Ω also G_fault=0.003781

#test implementation for line and transformer with PowerFactory Data
S_b = 100000000 #100MVA -> Woher? Ist das überhaupt korrekt? (in 9-Bus-System wird es so gemacht, aber warum?)
V_b = 230000 #Busspannung #230000 #230kV
ω_b = 60*2*π
line_length = 1
R_l_perkm = 5.29
X_l_perkm = 44.965
B_l_perkm = 332.7 * 10^(-6)#ω_b * 0.0095491* 10^(-6) #ω*C
R_pu = (R_l_perkm * line_length) * S_b/V_b^2
X_pu = (X_l_perkm * line_length) * S_b/V_b^2
B_pu = (B_l_perkm * line_length) * V_b^2/S_b
@named l45 = Line(piline(; R=R_pu, X=X_pu, B=B_pu), src=4, dst=5)

S_b = 100000000 #100MVA -> Bezugsscheinleistung
S_n = 200000000 #200MVA -> Bemessungsleistung
#V_OS = 230000
#V_US = 180000
X_pu_PF = 0.1250
R_pu_PF = 0
X_pu_baseS = X_pu_PF * S_b/S_n
R_pu_baseS = R_pu_PF * S_b/S_n
@named t27 = Line(transformer(; R=R_pu_baseS, X=X_pu_baseS), src=2, dst=7)

# build network
vertexfs = [bus1, bus2, bus3, bus4, bus5, bus6, bus7, bus8, bus9];
edgefs = [l45, l46, l69, l78, l89, t14, t27, t39, l57];
nw = Network(vertexfs, edgefs)

# solve powerflow and initialize
OpPoDyn.solve_powerflow!(nw)
#initialize_component!(bus1)
#dump_initial_state(bus1)
OpPoDyn.initialize!(nw)

# get state for actual calculation
u0 = NWState(nw)

# create faults
affect1! = (integrator) -> begin
    if integrator.t == 1
        @info "Short circuit on line 57 at t = $(integrator.t)"
        p = NWParameter(integrator)
        p.e[6, :pibranch₊shortcircuit] = 1
        auto_dt_reset!(integrator)
        save_parameters!(integrator)
    else
        error("Should not be reached.")
    end
end
cb_shortcircuit = PresetTimeCallback([1.0], affect1!)

affect2! = (integrator) -> begin
    if integrator.t == 1.05
        @info "Deactivate line 57 at t = $(integrator.t)"
        p = NWParameter(integrator)
        p.e[6, :pibranch₊active] = 0
        auto_dt_reset!(integrator)
        save_parameters!(integrator)
    else
        error("Should not be reached.")
    end
end
cb_deactivate = PresetTimeCallback([1.05], affect2!)

cb_set = CallbackSet(cb_shortcircuit, cb_deactivate)
prob = ODEProblem(nw, uflat(u0), (0,5), copy(pflat(u0)) ; callback=cb_set)
sol = solve(prob, Rodas5P(), dtmax=0.0001);
nothing

break # stop execution of script here

inspect(sol)

#### Voltage Magnitude
ref = CSV.read("Bus5-7_standardModelPF_avrAmplidyne.csv", DataFrame; header=2, decimal=',')
fig = Figure();
ax = Axis(fig[1, 1]; title="Bus voltage magnitude (Power Factory Standard Model)")
ts = range(sol.t[begin],sol.t[end],length=10000)
umag5 = sqrt.(sol(ts; idxs=VIndex(5, :busbar₊u_r)).^2 + sol(ts; idxs=VIndex(5, :busbar₊u_i)).^2)
umag7 = sqrt.(sol(ts; idxs=VIndex(7, :busbar₊u_r)).^2 + sol(ts; idxs=VIndex(7, :busbar₊u_i)).^2)
lines!(ax, ts, umag5.u; label="Bus5")
lines!(ax, ref."Zeitpunkt in s", ref."u1, Betrag in p.u._1", color=Cycled(1), linestyle=:dash, label="Bus 5 ref")
lines!(ax, ts, umag7.u; label="Bus7")
lines!(ax, ref."Zeitpunkt in s", ref."u1, Betrag in p.u.", color=Cycled(2), linestyle=:dash, label="Bus 7 ref")
axislegend(ax; position=:rb)
xlims!(ax, 0, 5)
fig



# Bus 2
#### id and iq generator
ref = CSV.read("gen2_data_avrAmplidyne.csv", DataFrame; header=2, decimal=',')
fig = Figure();
ax = Axis(fig[1, 1]; title="stator current gen 2")
ts = range(sol.t[begin],sol.t[end],length=10000)
id = sol(ts; idxs=VIndex(2, :ctrld_gen₊machine₊I_d))
iq = sol(ts; idxs=VIndex(2, :ctrld_gen₊machine₊I_q))
lines!(ax, ts, id.u; label="i_d")
lines!(ax, ref."Zeitpunkt in s", ref."Ständerstrom, d-Achse in p.u.", color=Cycled(1), linestyle=:dash, label="i_d ref")
lines!(ax, ts, iq.u; label="i_q")
lines!(ax, ref."Zeitpunkt in s", ref."Ständerstrom, q-Achse in p.u.", color=Cycled(2), linestyle=:dash, label="i_q ref")
axislegend(ax; position=:rt)
xlims!(ax, 0.9, 2)
fig

#### ud and uq generator
ref = CSV.read("gen2_data_avrAmplidyne.csv", DataFrame; header=2, decimal=',')
fig = Figure();
ax = Axis(fig[1, 1]; title="voltage at generator 2")
ts = range(sol.t[begin],sol.t[end],length=1000)
vd = sol(ts; idxs=VIndex(2, :ctrld_gen₊machine₊V_d))
vq = sol(ts; idxs=VIndex(2, :ctrld_gen₊machine₊V_q))
lines!(ax, ts, vd.u; label="u_d")
lines!(ax, ref."Zeitpunkt in s", ref."Spannung, d-Achse in p.u.", color=Cycled(1), linestyle=:dash, label="u_d ref")
lines!(ax, ts, vq.u; label="u_q")
lines!(ax, ref."Zeitpunkt in s", ref."Spannung, q-Achse in p.u.", color=Cycled(2), linestyle=:dash, label="u_q ref")
axislegend(ax; position=:rt)
xlims!(ax, 0.9, 2)
fig

#vr in OpPoDyn
ref = CSV.read("Gen2_standardModelPF_avrdata.csv", DataFrame; header=2, decimal=',', delim=';')
fig = Figure();
ax = Axis(fig[1, 1]; title="vr")
ts = range(sol.t[begin],sol.t[end],length=10000)
vr = sol(ts; idxs=VIndex(2, :ctrld_gen₊avr₊vr))
vfout = sol(ts; idxs=VIndex(2, :ctrld_gen₊avr₊vfout))
lines!(ax, ts, vr.u; label="vr")
lines!(ax, ref."Zeitpunkt in s", ref."o1", color=Cycled(1), linestyle=:dash, label="vr in PowerFactory")
#lines!(ax, ts, vfout.u; label="vfout")
axislegend(ax; position=:rb)
xlims!(ax, 0.9, 5)
fig

#vh.u
ref = CSV.read("Gen2_standardModelPF_avrdata.csv", DataFrame; header=2, decimal=',', delim=';')
fig = Figure();
ax = Axis(fig[1, 1]; title="vh.u")
ts = range(sol.t[begin],sol.t[end],length=10000)
vh = sol(ts; idxs=VIndex(2, :ctrld_gen₊machine₊v_mag))
lines!(ax, ts, vh.u; label="vh.u in OpPoDyn")
lines!(ax, ref."Zeitpunkt in s", ref."u", color=Cycled(1), linestyle=:dash, label="u in PowerFactory")
axislegend(ax; position=:rb)
xlims!(ax, 0.9, 5)
fig

#vref passt
ref = CSV.read("Gen2_standardModelPF_avrdata.csv", DataFrame; header=2, decimal=',', delim=';')
ref.summe = ref."upss" .+ ref."voel" .+ ref."vuel" .+ ref."avrref"
fig = Figure();
ax = Axis(fig[1, 1]; title="v_ref")
ts = range(sol.t[begin],sol.t[end],length=10000)
vref = sol(ts; idxs=VIndex(2, :ctrld_gen₊avr₊vref)) ##möglich, da eig _vref, aber vref_input=false
lines!(ax, ts, vref.u; label="v_ref in OpPoDyn")
lines!(ax, ref."Zeitpunkt in s", ref.summe, color=Cycled(1), linestyle=:dash, label="v_ref aus PowerFactory")
axislegend(ax; position=:rt)
xlims!(ax, 0, 5)
fig

#output
ref = CSV.read("Gen2_standardModelPF_avrdata.csv", DataFrame; header=2, decimal=',', delim=';')
fig = Figure();
ax = Axis(fig[1, 1]; title="vfout")
ts = range(sol.t[begin],sol.t[end],length=10000)
vfout = sol(ts; idxs=VIndex(2, :ctrld_gen₊avr₊vfout))
lines!(ax, ts, vfout.u; label="vfout")
lines!(ax, ref."Zeitpunkt in s", ref."uerrs", color=Cycled(1), linestyle=:dash, label="vfout in PowerFactory")
axislegend(ax; position=:rt)
xlims!(ax, 0.9, 5)
fig