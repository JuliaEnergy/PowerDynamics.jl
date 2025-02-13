using OpPoDyn
using OpPoDyn.Library
using ModelingToolkit
using NetworkDynamics
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

@mtkmodel ClassicBus begin
    @components begin
        machine = Library.ClassicalMachine_powerfactory(;
            p_m_input=false,
            S_b=100,
            V_b=18,
            ω_b=2π*60,
            X′_d,
            R_s=0.0026,
            vf_set=nothing,
            p_m_set=nothing,
            H,
             D=0
        )
        busbar = BusBar()
    end
    @equations begin
        connect(machine.terminal, busbar.terminal)
    end
end

@mtkmodel StandardBus begin
    @components begin
        machine = Library.StandardModel_pf_testneu(;
            S_b,
            V_b,
            Sn=S_b,
            Vn=V_b,
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
            pt,
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

#calculate parameters externally
function create_standardgenerator(; ω_b, H, S_b, V_b, D, R_s, X_rld, X_rlq, X_d, X_q, X′_d, X′_q=0.0001, X″_d, X″_q, X_ls, T′_d0, T″_d0, T′_q0=0.0001, T″_q0, cosn, dkd, dpe, salientpole=1, pt, dpu=0, addmt=0, xmdm=0)
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
    #X″_d = X_ad + X_ls - (k_1d + k_fd) * X_ad #??
    k_1qs = X_aq / (X_aq + X_rlq + X_1q) #salient pole
        #k_2qs = 0 #salient pole
    #X″_qs = X_aq + X_ls - k_1qs * X_aq #salient pole
    k_1qr = (X_aq * X_2q)/((X_aq + X_rlq) * (X_2q + X_1q) + X_2q * X_1q) #round rotor
    k_2qr = (X_aq * X_1q)/((X_aq + X_rlq) * (X_2q + X_1q) + X_2q * X_1q) #round rotor
    #X″_qr = X_aq + X_ls - (k_2qr + k_1qr) * X_aq #round rotor
    k_1q = salientpole * k_1qs + (1-salientpole) * k_1qr
    k_2q = (1-salientpole) * k_2qr
    #X″_q = salientpole * X″_qs + (1-salientpole)* X″_qr

    #(69), (71)
    X_det_d = (X_ad + X_rld) * (X_1d + X_fd) + X_fd * X_1d
    X_det_q = (X_aq + X_rlq) * (X_2q + X_1q) + X_2q * X_1q
    X_fd_loop = X_ad + X_rld + X_fd
    X_1d_loop = X_ad + X_rld + X_1d
    X_1q_loop = X_aq + X_rlq + X_1q
    X_2q_loop = X_aq + X_rlq + X_2q
    return ω_b, X_rld, X_rlq, X″_d, X″_q, X_ls, X_ad, X_aq, X_det_d, X_det_q, X_fd_loop, X_1d_loop, X_1q_loop, X_2q_loop, k_fd, k_1d, k_1q, k_2q, X_1q, R_1q, X_2q, R_2q, X_fd, R_1d, R_fd, X_1d, R_s, H, S_b, V_b, D, cosn, dkd, dpe, salientpole, pt, dpu, addmt, xmdm
end

# generate all MTK bus models
(ω_b, X_rld, X_rlq, X″_d, X″_q, X_ls, X_ad, X_aq, X_det_d, X_det_q, X_fd_loop, X_1d_loop, X_1q_loop, X_2q_loop, k_fd, k_1d, k_1q, k_2q, X_1q, R_1q, X_2q, R_2q, X_fd, R_1d, R_fd, X_1d, R_s, H, S_b, V_b, D, cosn, dkd, dpe, salientpole, pt, dpu, addmt, xmdm) =
    create_standardgenerator(; ω_b=2π*60, H=23.64, S_b=247.5, V_b=16.5, D=0, R_s=0, X_rld=0, X_rlq=0, X_d=0.36135, X_q=0.2398275, X′_d=0.15048, X′_q=0.0001, X″_d=0.1, X″_q=0.1, X_ls=0.08316, T′_d0=8.96, T″_d0=0.075, T′_q0=0.0001, T″_q0=0.15, cosn=1, dkd=0, dpe=0, salientpole=1, pt=0.2895, dpu=0, addmt=0, xmdm=0)
@named mtkbus1 = StandardBus(; machine__S_b=S_b, machine__V_b=V_b, machine__H=H, machine__D=D, machine__R_s=R_s, machine__X_rld=X_rld, machine__X_rlq=X_rlq, machine__X″_d=X″_d, machine__X″_q=X″_q, machine__X_ls=X_ls, machine__X_ad=X_ad, machine__X_aq=X_aq, machine__X_1q=X_1q, machine__X_det_d=X_det_d, machine__X_det_q=X_det_q, machine__X_fd_loop=X_fd_loop, machine__X_1d_loop=X_1d_loop, machine__X_1q_loop=X_1q_loop, machine__X_2q_loop=X_2q_loop, machine__k_fd=k_fd, machine__k_1d=k_1d, machine__k_1q=k_1q, machine__k_2q=k_2q, machine__R_fd=R_fd, machine__R_1d=R_1d, machine__R_1q=R_1q, machine__R_2q=R_2q, machine__cosn=cosn, machine__dkd=dkd, machine__dpe=dpe, machine__salientpole=salientpole, machine__dpu=dpu, machine__addmt=addmt, machine__xmdm=xmdm)

(ω_b, X_rld, X_rlq, X″_d, X″_q, X_ls, X_ad, X_aq, X_det_d, X_det_q, X_fd_loop, X_1d_loop, X_1q_loop, X_2q_loop, k_fd, k_1d, k_1q, k_2q, X_1q, R_1q, X_2q, R_2q, X_fd, R_1d, R_fd, X_1d, R_s, H, S_b, V_b, D, cosn, dkd, dpe, salientpole, pt, dpu, addmt, xmdm) =
    create_standardgenerator(; ω_b=2π*60, H=6.40, S_b=192, V_b=18, D=0, R_s=0.005, X_rld=0, X_rlq=0, X_d=1.719936, X_q=1.65984, X′_d=0.230016, X′_q=0.378048, X″_d=0.2, X″_q=0.2, X_ls=0.100032, T′_d0=6, T″_d0=0.0575, T′_q0=0.535, T″_q0=0.0945, cosn=0.85, dkd=0, dpe=0, salientpole=0, pt=1.003, dpu=0, addmt=0, xmdm=0)
@named mtkbus2 = StandardBus(; machine__S_b=S_b, machine__V_b=V_b, machine__H=H, machine__D=D, machine__R_s=R_s, machine__X_rld=X_rld, machine__X_rlq=X_rlq, machine__X″_d=X″_d, machine__X″_q=X″_q, machine__X_ls=X_ls, machine__X_ad=X_ad, machine__X_aq=X_aq, machine__X_1q=X_1q, machine__X_det_d=X_det_d, machine__X_det_q=X_det_q, machine__X_fd_loop=X_fd_loop, machine__X_1d_loop=X_1d_loop, machine__X_1q_loop=X_1q_loop, machine__X_2q_loop=X_2q_loop, machine__k_fd=k_fd, machine__k_1d=k_1d, machine__k_1q=k_1q, machine__k_2q=k_2q, machine__R_fd=R_fd, machine__R_1d=R_1d, machine__R_1q=R_1q, machine__R_2q=R_2q, machine__cosn=cosn, machine__dkd=dkd, machine__dpe=dpe, machine__salientpole=salientpole, machine__dpu=dpu, machine__addmt=addmt, machine__xmdm=xmdm)

(ω_b, X_rld, X_rlq, X″_d, X″_q, X_ls, X_ad, X_aq, X_det_d, X_det_q, X_fd_loop, X_1d_loop, X_1q_loop, X_2q_loop, k_fd, k_1d, k_1q, k_2q, X_1q, R_1q, X_2q, R_2q, X_fd, R_1d, R_fd, X_1d, R_s, H, S_b, V_b, D, cosn, dkd, dpe, salientpole, pt, dpu, addmt, xmdm) =
    create_standardgenerator(; ω_b=2π*60, H=3.01, S_b=128, V_b=13.8, D=0, R_s=0.0001, X_rld=0, X_rlq=0, X_d=1.68, X_q=1.609984, X′_d=0.232064, X′_q=0.32, X″_d=0.2, X″_q=0.2, X_ls=0.094976, T′_d0=5.89, T″_d0=0.0575, T′_q0=0.6, T″_q0=0.08, cosn=0.85, dkd=0, dpe=0, salientpole=0, pt=0.7813, dpu=0, addmt=0, xmdm=0)
@named mtkbus3 = StandardBus(; machine__S_b=S_b, machine__V_b=V_b, machine__H=H, machine__D=D, machine__R_s=R_s, machine__X_rld=X_rld, machine__X_rlq=X_rlq, machine__X″_d=X″_d, machine__X″_q=X″_q, machine__X_ls=X_ls, machine__X_ad=X_ad, machine__X_aq=X_aq, machine__X_1q=X_1q, machine__X_det_d=X_det_d, machine__X_det_q=X_det_q, machine__X_fd_loop=X_fd_loop, machine__X_1d_loop=X_1d_loop, machine__X_1q_loop=X_1q_loop, machine__X_2q_loop=X_2q_loop, machine__k_fd=k_fd, machine__k_1d=k_1d, machine__k_1q=k_1q, machine__k_2q=k_2q, machine__R_fd=R_fd, machine__R_1d=R_1d, machine__R_1q=R_1q, machine__R_2q=R_2q, machine__cosn=cosn, machine__dkd=dkd, machine__dpe=dpe, machine__salientpole=salientpole, machine__dpu=dpu, machine__addmt=addmt, machine__xmdm=xmdm)
#@named mtkbus1 = ClassicBus(; machine__H=23.64, machine__X′_d=0.0608)
#@named mtkbus2 = ClassicBus(; machine__H= 6.40, machine__X′_d=0.1198)
#@named mtkbus3 = ClassicBus(; machine__H= 3.01, machine__X′_d=0.1813)
@named mtkbus4 = MTKBus()
@named mtkbus5 = LoadBus(;load__Pset=-1.25, load__Qset=-0.5)
@named mtkbus6 = LoadBus(;load__Pset=-0.90, load__Qset=-0.3)
@named mtkbus7 = MTKBus()
@named mtkbus8 = LoadBus(;load__Pset=-1.0, load__Qset=-0.35)
@named mtkbus9 = MTKBus()


# generate the dynamic component functions
@named bus1 = Bus(mtkbus1; vidx=1, pf=pfSlack(V=1.04))
@named bus2 = Bus(mtkbus2; vidx=2, pf=pfPV(V=1.025, P=1.63))
@named bus3 = Bus(mtkbus3; vidx=3, pf=pfPV(V=1.025, P=0.85))
@named bus4 = Bus(mtkbus4; vidx=4)
@named bus5 = Bus(mtkbus5; vidx=5, pf=pfPQ(P=-1.25, Q=-0.5))
@named bus6 = Bus(mtkbus6; vidx=6, pf=pfPQ(P=-0.9, Q=-0.3))
@named bus7 = Bus(mtkbus7; vidx=7)
@named bus8 = Bus(mtkbus8; vidx=8, pf=pfPQ(P=-1.0, Q=-0.35))
@named bus9 = Bus(mtkbus9; vidx=9)

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


@named l45 = Line(piline(; R=0.0100, X=0.0850, B=0.1760), src=4, dst=5)
@named l46 = Line(piline(; R=0.0170, X=0.0920, B=0.1580), src=4, dst=6)
@named l69 = Line(piline(; R=0.0390, X=0.1700, B=0.3580), src=6, dst=9)
@named l78 = Line(piline(; R=0.0085, X=0.0720, B=0.1490), src=7, dst=8)
@named l89 = Line(piline(; R=0.0119, X=0.1008, B=0.2090), src=8, dst=9)
@named t14 = Line(transformer(; R=0, X=0.0576), src=1, dst=4)
@named t27 = Line(transformer(; R=0, X=0.0625), src=2, dst=7)
@named t39 = Line(transformer(; R=0, X=0.0586), src=3, dst=9)
@named l57 = Line(piline_shortcircuit(; R=0.0320, X=0.1610, B=0.3060, pos=0.99), src=5, dst=7) #S_b = 100 MVA, U_b = 230 kV; 2 Ω also G_fault=0.003781

# build network
vertexfs = [bus1, bus2, bus3, bus4, bus5, bus6, bus7, bus8, bus9];
edgefs = [l45, l46, l69, l78, l89, t14, t27, t39, l57];
nw = Network(vertexfs, edgefs)

# solve powerflow and initialize
OpPoDyn.solve_powerflow!(nw) #hier: ┌ Warning: Potential Rank Deficient Matrix Detected. Attempting to solve using Pivoted QR Factorization. └ @ NonlinearSolve C:\Users\smmywael\.julia\packages\NonlinearSolve\sETeN\src\internal\linear_solve.jl:162
OpPoDyn.initialize!(nw)

# get state for actual calculation
u0 = NWState(nw)

# create faults
affect1! = (integrator) -> begin
    if integrator.t == 1.0
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

cb_set = CallbackSet(CallbackSet(cb_shortcircuit, cb_deactivate))
prob = ODEProblem(nw, uflat(u0), (0,2), copy(pflat(u0)) ; callback=cb_set)
sol = solve(prob, Rodas5P());
nothing

break # stop execution of script here

#### Machine Angle
ref = CSV.read("RotorAngle_standardModelPF_shortcircuit.csv", DataFrame; header=2, decimal=',')
fig = Figure();
ax = Axis(fig[1, 1]; title="Läuferwinkel (Power Factory Standard Model)")
ts = range(sol.t[begin],sol.t[end],length=1000)
δ2 = sol(ts, idxs=VIndex(2, :machine₊δ)).u - sol(ts, idxs=VIndex(1, :machine₊δ)).u
δ3 = sol(ts, idxs=VIndex(3, :machine₊δ)).u - sol(ts, idxs=VIndex(1, :machine₊δ)).u
lines!(ax, ts, rad2deg.(δ2); label="Bus 2")
lines!(ax, ref."Zeitpunkt in s", ref."firel in deg", color=Cycled(1), linestyle=:dash, label="Bus 2 ref")
lines!(ax, ts, rad2deg.(δ3); label="Bus 3")
lines!(ax, ref."Zeitpunkt in s", ref."firel in deg_1", color=Cycled(2), linestyle=:dash, label="Bus 3 ref")
axislegend(ax; position=:lt)
fig


#### Voltage Magnitude
ref = CSV.read("Bus5-7_standardModelPF_shortcircuit.csv", DataFrame; header=2, decimal=',')
fig = Figure();
ax = Axis(fig[1, 1]; title="Bus voltage magnitude (Power Factory Standard Model)")
ts = range(sol.t[begin],sol.t[end],length=1000)
umag5 = sqrt.(sol(ts; idxs=VIndex(5, :busbar₊u_r)).^2 + sol(ts; idxs=VIndex(5, :busbar₊u_i)).^2)
umag7 = sqrt.(sol(ts; idxs=VIndex(7, :busbar₊u_r)).^2 + sol(ts; idxs=VIndex(7, :busbar₊u_i)).^2)
lines!(ax, ts, umag5.u; label="Bus5")
lines!(ax, ref."Zeitpunkt in s", ref."u1, Betrag in p.u._1", color=Cycled(1), linestyle=:dash, label="Bus 5 ref")
lines!(ax, ts, umag7.u; label="Bus7")
lines!(ax, ref."Zeitpunkt in s", ref."u1, Betrag in p.u.", color=Cycled(2), linestyle=:dash, label="Bus 7 ref")
axislegend(ax; position=:rb)
xlims!(ax, 0.9, 2)
fig


#=
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
=#
