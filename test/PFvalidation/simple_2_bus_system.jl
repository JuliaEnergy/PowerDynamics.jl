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
using LinearAlgebra

# load reference data Data_all_PF_shortcircuit
gen_dat, linea_dat, lineb_dat, load_dat = let
    file = joinpath(pkgdir(OpPoDyn),"2bustest","Data_all_PF_shortcircuit_withQ.csv")
    headers = split(readline(file), ';')
    gen_cols = findall(c -> c ∈ ("Alle Berechnungsarten", "gen1"), headers)
    linea_cols = findall(c -> c ∈ ("Alle Berechnungsarten", "L1-2a"), headers)
    lineb_cols = findall(c -> c ∈ ("Alle Berechnungsarten", "L1-2b"), headers)
    load_cols = findall(c -> c ∈ ("Alle Berechnungsarten", "load1"), headers)

    all_data = CSV.read(file, DataFrame; delim=';', decimal=',', header=2)
    all_data[:, gen_cols], all_data[:, linea_cols], all_data[:, lineb_cols], all_data[:, load_cols]
end;

@mtkmodel LoadBus begin
    @components begin
        busbar = BusBar()
        load = ConstantYLoad(Pset, Qset, Vset=nothing)
        # load = ConstantYLoad(Pset, Vset=nothing)
    end
    @equations begin
        connect(load.terminal, busbar.terminal)
    end
end

@mtkmodel StandardBus begin
    @components begin
        machine = Library.StandardModel_pf_testneu(;
            S_b,
            Sn,
            V_b,
            Vn,
            ω_b,
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
            # X_fd,
            # X_1d,
            # X_2q,
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
            xmdm,
            #τ_m_test_in,
            #vf_test_in
        )
        busbar = BusBar()
    end
    @equations begin
        connect(machine.terminal, busbar.terminal)
    end
end

#calculate parameters externally
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
primary_parameters = (;
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
    pt=0.20307028,
    dpu=0,
    addmt=0,
    xmdm=0,
)
secondary_parameters = secondary_from_primary(; primary_parameters...)

@warn "Symbols $(keys(primary_parameters) ∩ keys(secondary_parameters)) are overwrittenby secondary parameters, but didn't really change..."
#primary_parameters.X″_d - secondary_parameters.X″_d
#primary_parameters.X″_q - secondary_parameters.X″_q

allp = Dict(pairs(primary_parameters)..., pairs(secondary_parameters)...)
# get rid of parametes which are not needed
for p in [:X_d, :X′_q, :T′_d0, :T′_q0, :X′_d, :T″_q0, :T″_d0, :X_q, :X_fd, :X_1d, :X_2q]
    delete!(allp, p)
end
renamedp = Dict(map(s->Symbol("machine__", s), collect(keys(allp))) .=> Base.values(allp))


# (ω_b, X_rld, X_rlq, X″_d, X″_q, X_ls, X_ad, X_aq, X_det_d, X_det_q, X_fd_loop, X_1d_loop, X_1q_loop, X_2q_loop, k_fd, k_1d, k_1q, k_2q, X_1q, R_1q, X_2q, R_2q, X_fd, R_1d, R_fd, X_1d, R_s, H, S_b, V_b, D, cosn, dkd, dpe, salientpole, pt, dpu, addmt, xmdm, X′_d, X′_q, X_d, X_q) =
#     create_standardgenerator(; ω_b=2π*60, H=9.551516, S_b=247.5, V_b=16.5, D=0, R_s=0, X_rld=0, X_rlq=0, X_d=0.36135, X_q=0.2398275, X′_d=0.15048, X′_q=0.0001, X″_d=0.1, X″_q=0.1, X_ls=0.08316, T′_d0=8.96, T″_d0=0.075, T′_q0=0.0001, T″_q0=0.15, cosn=1, dkd=0, dpe=0, salientpole=1, pt=0.2895, dpu=0, addmt=0, xmdm=0)
# variable_names = [:X_ad, :X_aq, :X_ls, :R_s, :X_rld, :X_rlq, :X″_d, :X″_q, :X′_d, :X′_q, :X_d, :X_q, :X_fd, :X_1d, :X_1q, :X_2q, :R_fd, :R_1d, :R_1q, :R_2q, :X_det_d, :X_det_q, :X_fd_loop, :X_1d_loop, :X_1q_loop, :X_2q_loop, :k_fd, :k_1d, :k_1q, :k_2q]
# # values = [X_ad, X_aq, X_ls, R_s, X_rld, X_rlq, X″_d, X″_q, X′_d, X′_q, X_d, X_q, X_fd, X_1d, X_1q, X_2q, R_fd, R_1d, R_1q, R_2q, X_det_d, X_det_q, X_fd_loop, X_1d_loop, X_1q_loop, X_2q_loop, k_fd, k_1d, k_1q, k_2q]
# df = DataFrame(Variable=variable_names, Value=values)
# formatted_df = DataFrame(
#     Variable=df.Variable,
#     Value=map(x -> replace(string(round(x, digits=8)), "." => ","), df.Value)
# )
# transposed_df = permutedims(formatted_df)
# # CSV.write(joinpath("2bustest", "gen1_parameter_oppodyn_2bus_basecase.csv"), transposed_df; delim=';')

# Branches

function piline_parallel(; Ra, Xa, Ba, posa, Ga_fault=0, Ba_fault=0, Rb, Xb, Bb)
    @named pibranch1 = PiLine_fault(;R=Ra, X=Xa, B_src=Ba/2, B_dst=Ba/2, G_src=0, G_dst=0, G_fault=Ga_fault, B_fault=Ba_fault, pos=posa)#, faultimp)
    @named pibranch2 = PiLine(; R=Rb, X=Xb, B_src=Bb/2, B_dst=Bb/2, G_src=0, G_dst=0)
    MTKLine(pibranch1, pibranch2)
end

# define line based on reference data...
Za_from_csv = 100*linea_dat[1, "Impedanz (bus1), mag in p.u."]*exp(im*deg2rad(linea_dat[1, "Impedanz (bus1), phi in deg"])) + 100*linea_dat[1, "Impedanz (bus2), mag in p.u."]*exp(im*deg2rad(linea_dat[1, "Impedanz (bus2), phi in deg"]))
#@named l12a = Line(piline_shortcircuit(; R=real(Za_from_csv), X=imag(Za_from_csv), B=0, pos=0.5), src=1, dst=2)
Zb_from_csv = 100*lineb_dat[1, "Impedanz, mag in p.u."]*exp(im*deg2rad(lineb_dat[1, "Impedanz, phi in deg"])) #TODO factor 100??
#@named l12b = Line(piline(; R=real(Zb_from_csv), X=imag(Zb_from_csv), B=0), src=1, dst=2) #in umgekehrter Reihenfolge von src und dst bekomme ich Fehler
@named l12 = Line(piline_parallel(; Ra=real(Za_from_csv), Xa=imag(Za_from_csv), Ba=0, posa=0.5, Ga_fault=0, Ba_fault=0, Rb=real(Zb_from_csv), Xb=imag(Zb_from_csv), Bb=0), src=1, dst=2)

@named mtkbus1 = StandardBus(; renamedp...)
@named mtkbus2 = LoadBus(;load__Pset=-0.5, load__Qset=-0.1) #50MW bzw 10MVar bezogen auf 100MVA

# generate the dynamic component functions
@named bus1 = Bus(mtkbus1; vidx=1, pf=pfSlack(V=1.0))
@named bus2 = Bus(mtkbus2; vidx=2, pf=pfPQ(P=-0.5, Q=-0.1))

# build network
vertexfs = [bus1, bus2];
edgefs = [l12]
nw = Network(vertexfs, edgefs)

# solve powerflow and initialize
OpPoDyn.solve_powerflow!(nw)
OpPoDyn.initialize!(nw)

####
#### funktionierende states
####
udq_pf = [gen_dat[1, "Spannung, d-Achse in p.u."],  gen_dat[1, "Spannung, q-Achse in p.u."]]
udq = get_initial_state(bus1, [:machine₊V_d,:machine₊V_q])
# WORKS

ψdq_pf = [gen_dat[1, "Ständerfluss, d-Achse"], gen_dat[1, "Ständerfluss, q-Achse in p.u."]]
ψdq = get_initial_state(bus1, [:machine₊ψ_d,:machine₊ψ_q])
# works!

V″dq_pf = [gen_dat[1, "Subtransiente Spannung, d-Achse in p.u."], gen_dat[1, "Subtransiente Spannung, q-Achse in p.u."]]
V″dq = get_initial_state(bus1, [:machine₊V″_d,:machine₊V″_q])
# close!

ψ″dq_pf = [gen_dat[1, "Subtransienter Fluss, d-Achse in p.u."], gen_dat[1, "Subtransienter Fluss, q-Achse in p.u."]]
ψ″dq = get_initial_state(bus1, [:machine₊ψ″_d,:machine₊ψ″_q])
# close!

# mechanische leistung
Pm_pf = gen_dat[1, "Turbinenleistung in p.u."]
Pm = get_initial_state(bus1, :machine₊τ_m)/get_initial_state(bus1, :machine₊n)

# mechanical torque
τ_m_pf = gen_dat[1, "Mechanisches Moment in p.u."]
τ_m = get_initial_state(bus1, :machine₊τ_m)

# electrical torque
τ_e_pf = gen_dat[1, "Elektrisches Moment in p.u."]
τ_e = get_initial_state(bus1, :machine₊τ_e)

# flux through excitation
ψ_fd_pf = gen_dat[1, "Fluss in Erregerwicklung in p.u."]
ψ_fd = get_initial_state(bus1, :machine₊ψ_fd)

# excitation current
I_fd_pf = gen_dat[1, "Strom in Erregerwicklung in p.u."]
I_fd = get_initial_state(bus1, :machine₊I_fd)

# damper currents
Idamp_pf = [gen_dat[1,"Strom in 1d-Dämpferwicklung in p.u."], gen_dat[1,"Strom in 1q-Dämpferwicklung in p.u."], gen_dat[1,"Strom in 2q-Dämpferwicklung in p.u."]]
Idamp = get_initial_state(bus1, [:machine₊I_1d,:machine₊I_1q,:machine₊I_2q])

# damper flux
ψdamp_pf = [gen_dat[1,"Fluss in 1d-Dämpferwicklung in p.u."], gen_dat[1,"Fluss in 1q-Dämpferwicklung in p.u."], gen_dat[1,"Fluss in 2q-Dämpferwicklung in p.u."]]
ψdamp = get_initial_state(bus1, [:machine₊ψ_1d,:machine₊ψ_1q,:machine₊ψ_2q])


####
#### problematische states
####
idq_pf = [gen_dat[1, "Ständerstrom, d-Achse in p.u."], gen_dat[1, "Ständerstrom, q-Achse in p.u."]]
idq = get_initial_state(bus1, [:machine₊I_d,:machine₊I_q])
# works now

# erregerspannung
vfd_pf = gen_dat[1, "Erregerspannung in p.u."] / gen_dat[1, "Konvertierungsfaktor (=xadu/rfd)"]
vf = get_initial_state(bus1, :machine₊vf)
# FIXME: something wrong: PF factor 1000 bigger than OpPoDyn

# parameter
gen_dat[1, "Reaktanz der Erregerwicklung in p.u."]
get_initial_state(bus1, :machine₊X_fd)
# FIXME: X_fd tauch nicht  mehr auf

gen_dat[1, "Widerstand der Erregerwicklung in p.u."]
get_initial_state(bus1, :machine₊R_fd)
# correct


# dump all initial conditions
for n in names(gen_dat)
    contains(n, "rre") || continue
    print(n, " = ")
    println(gen_dat[1, n])
end


u0 = NWState(nw)

affect1! = (integrator) -> begin
    if integrator.t == 10
        @info "Short circuit on line 1-2a at t = $(integrator.t)"
        p = NWParameter(integrator)
        p.e[1, :pibranch1₊shortcircuit] = 1
        auto_dt_reset!(integrator)
        save_parameters!(integrator)
    else
        error("Should not be reached.")
    end
end
cb_shortcircuit = PresetTimeCallback([10], affect1!)

affect2! = (integrator) -> begin
    if integrator.t == 12.5
        @info "Deactivate line 1-2a at t = $(integrator.t)"
        p = NWParameter(integrator)
        p.e[1, :pibranch1₊active] = 0
        auto_dt_reset!(integrator)
        save_parameters!(integrator)
    else
        error("Should not be reached.")
    end
end
cb_deactivate = PresetTimeCallback([12.5], affect2!)

cb_set = CallbackSet(cb_shortcircuit, cb_deactivate)
prob = ODEProblem(nw, uflat(u0), (0,12.9), copy(pflat(u0)) ; callback=cb_set)
sol = solve(prob, Rodas5P());
nothing

break

freq = @obsex (VIndex(1,:machine₊n) * 60)  #@obsex (VIndex(1,:machine₊n) *VIndex(1,:machine₊H) + VIndex(2,:machine₊n)*VIndex(2,:machine₊H)) / (VIndex(1,:machine₊H) + VIndex(2,:machine₊H))
plot(sol, idxs=freq; label="OpPoDyn")

ref = CSV.read("2bustest/frequency_PF_shortcircuit_withQ.csv", DataFrame; header=2, decimal=',')
fig = Figure();
ax = Axis(fig[1, 1]; title="Frequency")
ts = range(sol.t[begin],sol.t[end],length=1000)
f_oppodyn = round.(sol(ts; idxs=VIndex(1, :machine₊n)).*60, digits=8)
#lines!(ax, sol;idxs=@obsex (VIndex(1,:machine₊n) * 60), label="OpPoDyn")
lines!(ax, ts, f_oppodyn.u; label="OpPoDyn")
lines!(ax, ref."Zeitpunkt in s", ref."Elektrische Frequenz in Hz", color=Cycled(1), linestyle=:dash, label="Power Factory")
axislegend(ax; position=:rb)
xlims!(ax, 9.9, 12.9)
ylims!(ax, 59.9, 60.3)
fig

#### id and iq generator
ref = CSV.read("2bustest/Data_all_PF_shortcircuit_withQ.csv", DataFrame; header=2, decimal=',')
fig = Figure();
ax = Axis(fig[1, 1]; title="stator current")
ts = range(sol.t[begin],sol.t[end],length=1000)
id = sol(ts; idxs=VIndex(1, :machine₊I_d))
iq = sol(ts; idxs=VIndex(1, :machine₊I_q))
lines!(ax, ts, id.u; label="i_d")
lines!(ax, ref."Zeitpunkt in s", ref."Ständerstrom, d-Achse in p.u.", color=Cycled(1), linestyle=:dash, label="i_d ref")
lines!(ax, ts, iq.u; label="i_q")
lines!(ax, ref."Zeitpunkt in s", ref."Ständerstrom, q-Achse in p.u.", color=Cycled(2), linestyle=:dash, label="i_q ref")
axislegend(ax; position=:rt)
xlims!(ax, 9.9, 12.9)
ylims!(ax, 0, 4)
fig

#### ud and uq generator
ref = CSV.read("2bustest/Data_all_PF_shortcircuit_withQ.csv", DataFrame; header=2, decimal=',')
fig = Figure();
ax = Axis(fig[1, 1]; title="voltage at generator")
ts = range(sol.t[begin],sol.t[end],length=1000)
vd = sol(ts; idxs=VIndex(1, :machine₊V_d))
vq = sol(ts; idxs=VIndex(1, :machine₊V_q))
lines!(ax, ts, vd.u; label="u_d")
lines!(ax, ref."Zeitpunkt in s", ref."Spannung, d-Achse in p.u.", color=Cycled(1), linestyle=:dash, label="u_d ref")
lines!(ax, ts, vq.u; label="u_q")
lines!(ax, ref."Zeitpunkt in s", ref."Spannung, q-Achse in p.u.", color=Cycled(2), linestyle=:dash, label="u_q ref")
axislegend(ax; position=:rt)
xlims!(ax, 9.9, 12.9)
ylims!(ax, 0, 1.1)
fig

#### generator torque
ref = CSV.read("2bustest/Data_all_PF_shortcircuit_withQ.csv", DataFrame; header=2, decimal=',')
fig = Figure();
ax = Axis(fig[1, 1]; title="torque")
ts = range(sol.t[begin],sol.t[end],length=1000)
τ_m = sol(ts; idxs=VIndex(1, :machine₊τ_m))
τ_e = sol(ts; idxs=VIndex(1, :machine₊τ_e))
lines!(ax, ts, τ_m.u; label="τ_m")
lines!(ax, ref."Zeitpunkt in s", ref."Mechanisches Moment in p.u.", color=Cycled(1), linestyle=:dash, label="τ_m ref")
lines!(ax, ts, τ_e.u; label="τ_e")
lines!(ax, ref."Zeitpunkt in s", ref."Elektrisches Moment in p.u.", color=Cycled(2), linestyle=:dash, label="τ_e ref")
axislegend(ax; position=:lb)
xlims!(ax, 9, 12.9)
#ylims!(ax, 0, 1.1)
fig


#### Voltage Magnitude
ref = CSV.read("2bustest/voltage_PF_shortcircuit_withQ.csv", DataFrame; header=2, decimal=',')
fig = Figure();
ax = Axis(fig[1, 1]; title="Bus voltage magnitude")
ts = range(sol.t[begin],sol.t[end],length=1000)
umag1 = sqrt.(sol(ts; idxs=VIndex(1, :busbar₊u_r)).^2 + sol(ts; idxs=VIndex(1, :busbar₊u_i)).^2)
umag2 = sqrt.(sol(ts; idxs=VIndex(2, :busbar₊u_r)).^2 + sol(ts; idxs=VIndex(2, :busbar₊u_i)).^2)
lines!(ax, ts, umag1.u; label="Bus1")
lines!(ax, ref."Zeitpunkt in s", ref."u1, Betrag in p.u.", color=Cycled(1), linestyle=:dash, label="Bus 1 ref")
lines!(ax, ts, umag2.u; label="Bus2")
lines!(ax, ref."Zeitpunkt in s", ref."u1, Betrag in p.u._1", color=Cycled(2), linestyle=:dash, label="Bus 2 ref")
axislegend(ax; position=:rt)
xlims!(ax, 9, 12.9)
fig

#### Voltage angle
ref = CSV.read("2bustest/voltage_PF_shortcircuit_withQ.csv", DataFrame; header=2, decimal=',')
fig = Figure();
ax = Axis(fig[1, 1]; title="Bus voltage angle")
ts = range(sol.t[begin],sol.t[end],length=1000)
uang1 = rad2deg.(atan.(sol(ts; idxs=VIndex(1, :busbar₊u_i)), sol(ts; idxs=VIndex(1, :busbar₊u_r))))
uang2 = rad2deg.(atan.(sol(ts; idxs=VIndex(2, :busbar₊u_i)), sol(ts; idxs=VIndex(2, :busbar₊u_r))))
lines!(ax, ts, uang1.u; label="Bus1")
lines!(ax, ref."Zeitpunkt in s", ref."U, Winkel in deg", color=Cycled(1), linestyle=:dash, label="Bus 1 ref")
lines!(ax, ts, uang2.u; label="Bus2")
lines!(ax, ref."Zeitpunkt in s", ref."U, Winkel in deg_1", color=Cycled(2), linestyle=:dash, label="Bus 2 ref")
axislegend(ax; position=:lt)
xlims!(ax, 0, 12.9)
fig