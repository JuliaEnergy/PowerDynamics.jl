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

@mtkmodel StandardBus begin
    @components begin
        machine = Library.StandardModel_pf_testneu(;
            S_b,
            V_b,
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
            X_fd,
            X_1d,
            X_2q,
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
    X″_d = X_ad + X_ls - (k_1d + k_fd) * X_ad #??
    k_1qs = X_aq / (X_aq + X_rlq + X_1q) #salient pole
        #k_2qs = 0 #salient pole
    X″_qs = X_aq + X_ls - k_1qs * X_aq #salient pole
    k_1qr = (X_aq * X_2q)/((X_aq + X_rlq) * (X_2q + X_1q) + X_2q * X_1q) #round rotor
    k_2qr = (X_aq * X_1q)/((X_aq + X_rlq) * (X_2q + X_1q) + X_2q * X_1q) #round rotor
    X″_qr = X_aq + X_ls - (k_2qr + k_1qr) * X_aq #round rotor
    k_1q = salientpole * k_1qs + (1-salientpole) * k_1qr
    k_2q = (1-salientpole) * k_2qr
    X″_q = salientpole * X″_qs + (1-salientpole)* X″_qr

    #(69), (71)
    X_det_d = (X_ad + X_rld) * (X_1d + X_fd) + X_fd * X_1d 
    X_det_q = (X_aq + X_rlq) * (X_2q + X_1q) + X_2q * X_1q
    X_fd_loop = X_ad + X_rld + X_fd
    X_1d_loop = X_ad + X_rld + X_1d
    X_1q_loop = X_aq + X_rlq + X_1q 
    X_2q_loop = X_aq + X_rlq + X_2q
    return ω_b, X_rld, X_rlq, X″_d, X″_q, X_ls, X_ad, X_aq, X_det_d, X_det_q, X_fd_loop, X_1d_loop, X_1q_loop, X_2q_loop, k_fd, k_1d, k_1q, k_2q, X_1q, R_1q, X_2q, R_2q, X_fd, R_1d, R_fd, X_1d, R_s, H, S_b, V_b, D, cosn, dkd, dpe, salientpole, pt, dpu, addmt, xmdm, X′_d, X′_q, X_d, X_q
end 

# generate all MTK bus models
#aus PF Implementierung
(ω_b, X_rld, X_rlq, X″_d, X″_q, X_ls, X_ad, X_aq, X_det_d, X_det_q, X_fd_loop, X_1d_loop, X_1q_loop, X_2q_loop, k_fd, k_1d, k_1q, k_2q, X_1q, R_1q, X_2q, R_2q, X_fd, R_1d, R_fd, X_1d, R_s, H, S_b, V_b, D, cosn, dkd, dpe, salientpole, pt, dpu, addmt, xmdm, X′_d, X′_q, X_d, X_q) = 
    create_standardgenerator(; ω_b=2π*60, H=9.551516, S_b=247.5, V_b=16.5, D=0, R_s=0, X_rld=0, X_rlq=0, X_d=0.36135, X_q=0.2398275, X′_d=0.15048, X′_q=0.0001, X″_d=0.1, X″_q=0.1, X_ls=0.08316, T′_d0=8.96, T″_d0=0.075, T′_q0=0.0001, T″_q0=0.15, cosn=1, dkd=0, dpe=0, salientpole=1, pt=0.2895, dpu=0, addmt=0, xmdm=0)
variable_names = [:X_ad, :X_aq, :X_ls, :R_s, :X_rld, :X_rlq, :X″_d, :X″_q, :X′_d, :X′_q, :X_d, :X_q, :X_fd, :X_1d, :X_1q, :X_2q, :R_fd, :R_1d, :R_1q, :R_2q, :X_det_d, :X_det_q, :X_fd_loop, :X_1d_loop, :X_1q_loop, :X_2q_loop, :k_fd, :k_1d, :k_1q, :k_2q]
values = [X_ad, X_aq, X_ls, R_s, X_rld, X_rlq, X″_d, X″_q, X′_d, X′_q, X_d, X_q, X_fd, X_1d, X_1q, X_2q, R_fd, R_1d, R_1q, R_2q, X_det_d, X_det_q, X_fd_loop, X_1d_loop, X_1q_loop, X_2q_loop, k_fd, k_1d, k_1q, k_2q]
df = DataFrame(Variable=variable_names, Value=values)
formatted_df = DataFrame(
    Variable=df.Variable, 
    Value=map(x -> replace(string(round(x, digits=8)), "." => ","), df.Value)
)
transposed_df = permutedims(formatted_df)
CSV.write(joinpath("2bustest", "gen1_parameter_oppodyn_2bus_basecase.csv"), transposed_df; delim=';')

@named mtkbus1 = StandardBus(; machine__S_b=S_b, machine__V_b=V_b, machine__H=H, machine__D=D, machine__R_s=R_s, machine__X_rld=X_rld, machine__X_rlq=X_rlq, machine__X″_d=X″_d, machine__X″_q=X″_q, machine__X_ls=X_ls, machine__X_ad=X_ad, machine__X_aq=X_aq, machine__X_1d=X_1d, machine__X_1q=X_1q, machine__X_2q=X_2q, machine__X_fd=X_fd, machine__X_det_d=X_det_d, machine__X_det_q=X_det_q, machine__X_fd_loop=X_fd_loop, machine__X_1d_loop=X_1d_loop, machine__X_1q_loop=X_1q_loop, machine__X_2q_loop=X_2q_loop, machine__k_fd=k_fd, machine__k_1d=k_1d, machine__k_1q=k_1q, machine__k_2q=k_2q, machine__R_fd=R_fd, machine__R_1d=R_1d, machine__R_1q=R_1q, machine__R_2q=R_2q, machine__cosn=cosn, machine__dkd=dkd, machine__dpe=dpe, machine__salientpole=salientpole, machine__pt=pt, machine__dpu=dpu, machine__addmt=addmt, machine__xmdm=xmdm)
@named mtkbus2 = LoadBus(;load__Pset=-0.5, load__Qset=0)

# generate the dynamic component functions
@named bus1 = Bus(mtkbus1; vidx=1, pf=pfSlack(V=1.0)) 
@named bus2 = Bus(mtkbus2; vidx=2,  pf=pfPQ(P=-0.5, Q=0))

# Branches
function piline(; R, X, B)
    @named pibranch = PiLine(;R, X, B_src=B/2, B_dst=B/2, G_src=0, G_dst=0)
    MTKLine(pibranch)
end

function piline_shortcircuit(; R, X, B, pos, G_fault=0, B_fault=0)
    @named pibranch = PiLine_fault(;R, X, B_src=B/2, B_dst=B/2, G_src=0, G_dst=0, G_fault, B_fault, pos)
    MTKLine(pibranch)
end

function transformer(; R, X)
    @named transformer = PiLine(;R, X, B_src=0, B_dst=0, G_src=0, G_dst=0)
    MTKLine(transformer)
end

S_b = 100000000 #100MVA -> Woher? Ist das überhaupt korrekt? (in 9-Bus-System wird es so gemacht, aber warum?)
V_b = 230000 #230kV
ω_b = 60*2*π
line_length = 1
R_l_perkm = 5.29
X_l_perkm = 44.965
B_l_perkm = 332.7 * 10^(-6)#ω_b * 0.0095491* 10^(-6) #ω*C
R_pu = (R_l_perkm * line_length) * S_b/V_b^2
X_pu = (X_l_perkm * line_length) * S_b/V_b^2
B_pu = (B_l_perkm * line_length) * V_b^2/S_b
@named l12 = Line(piline(; R=R_pu, X=X_pu, B=B_pu), src=1, dst=2)


# build network
vertexfs = [bus1, bus2];
edgefs = [l12];
nw = Network(vertexfs, edgefs)

# solve powerflow and initialize
OpPoDyn.solve_powerflow!(nw)
OpPoDyn.initialize!(nw)

# get state for actual calculation
u0 = NWState(nw)

# create faults
affect1! = (integrator) -> begin
    if integrator.t == 0.0
        @info "Short circuit on line 57 at t = $(integrator.t)"
        p = NWParameter(integrator)
        p.e[6, :pibranch₊shortcircuit] = 1
    else
        error("Should not be reached.")
    end
end
cb_shortcircuit = PresetTimeCallback([0.0], affect1!) 

affect2! = (integrator) -> begin
    if integrator.t == 0.05
        @info "Deactivate line 57 at t = $(integrator.t)"
        p = NWParameter(integrator)
        p.e[6, :pibranch₊active] = 0
    else
        error("Should not be reached.")
    end
end
cb_deactivate = PresetTimeCallback([0.05], affect2!)

cb_set = CallbackSet() #CallbackSet(cb_shortcircuit, cb_deactivate)
prob = ODEProblem(nw, uflat(u0), (0,1), copy(pflat(u0)) ; callback=cb_set)
sol = solve(prob, Rodas5P());
nothing

break # stop execution of script here



ts=0
n=1
variable_names = [:machine₊τ_m, :machine₊τ_e, :machine₊I_d, :machine₊I_q, :machine₊I_fd, :machine₊I_1d, :machine₊I_1q, :machine₊I_2q, :machine₊ψ_d, :machine₊ψ_q, :machine₊ψ″_d, :machine₊ψ″_q, :machine₊ψ_fd, :machine₊ψ_1d, :machine₊ψ_1q, :machine₊ψ_2q, :machine₊V_d, :machine₊V_q, :machine₊V″_d, :machine₊V″_q, :machine₊vf]
values = [sol(ts, idxs=VIndex(n, :machine₊τ_m)), sol(ts, idxs=VIndex(n, :machine₊τ_e)), sol(ts, idxs=VIndex(n, :machine₊I_d)), sol(ts, idxs=VIndex(n, :machine₊I_q)), sol(ts, idxs=VIndex(n, :machine₊I_fd)), sol(ts, idxs=VIndex(n, :machine₊I_1d)), sol(ts, idxs=VIndex(n, :machine₊I_1q)), sol(ts, idxs=VIndex(n, :machine₊I_2q)), sol(ts, idxs=VIndex(n, :machine₊ψ_d)), sol(ts, idxs=VIndex(n, :machine₊ψ_q)), sol(ts, idxs=VIndex(n, :machine₊ψ″_d)), sol(ts, idxs=VIndex(n, :machine₊ψ″_q)), sol(ts, idxs=VIndex(n, :machine₊ψ_fd)), sol(ts, idxs=VIndex(n, :machine₊ψ_1d)), sol(ts, idxs=VIndex(n, :machine₊ψ_1q)), sol(ts, idxs=VIndex(n, :machine₊ψ_2q)), sol(ts, idxs=VIndex(n, :machine₊V_d)), sol(ts, idxs=VIndex(n, :machine₊V_q)), sol(ts, idxs=VIndex(n, :machine₊V″_d)), sol(ts, idxs=VIndex(n, :machine₊V″_q)), sol(ts, idxs=VIndex(n, :machine₊vf)) ]
df = DataFrame(Variable=variable_names, Value=values)
formatted_df = DataFrame(
    Variable=df.Variable,
    Value=map(x -> replace(string(round(x, digits=8)), "." => ","), df.Value))
transposed_df = permutedims(formatted_df)
CSV.write(joinpath("2bustest", "gen1_variableresults_oppodyn_2bus_basecase.csv"), transposed_df; delim=";")
