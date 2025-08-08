using PowerDynamics
using NetworkDynamics
using ModelingToolkit
using ModelingToolkit: D_nounits as Dt, t_nounits as t
using CSV
using SteadyStateDiffEq

ω0    = 2π*50    # rad/s
Sbase = 300      # MW
Vbase = 110      # kV

Rline = 1        # Ω
Lline = (1/100π) # H
Cline = (2e-6)   # F
Pload = -300     # MW

Rline_pu = Rline / Zbase(Sbase, Vbase)
Lline_pu = Lline / Zbase(Sbase, Vbase)
Cline_pu = Cline / Ybase(Sbase, Vbase)
Pload_pu = Pload / Sbase


slackbus = Bus(pfSlack(; V=1), vidx=1)

@mtkmodel DynamicShunt begin
    @components begin
        terminal = Terminal()
    end
    @variables begin
        u_r(t), [guess=1, description="Real part of voltage"]
        u_i(t), [guess=0, description="Imaginary part of voltage"]
        ## voltage in abc as observable
        u_a(t), [description="Voltage in a phase"]
        u_b(t), [description="Voltage in b phase"]
        u_c(t), [description="Voltage in c phase"]
    end
    @parameters begin
        C, [description="Capacitance"]
        ω0, [description="Angular frequency of dq Frame"]
    end
    begin
        Tdqinv(δ) =  #=√(2/3)*=# [cos(δ)       -sin(δ)
                             cos(δ-2pi/3) -sin(δ-2pi/3)
                             cos(δ+2pi/3) -sin(δ+2pi/3)]
    end
    @equations begin
        Dt(u_r) ~  ω0*u_i + 1/C * terminal.i_r
        Dt(u_i) ~ -ω0*u_r + 1/C * terminal.i_i
        ## connection to terminal
        terminal.u_r ~ u_r
        terminal.u_i ~ u_i
        ## abc voltages
        [u_a, u_b, u_c] ~ Tdqinv(ω0*t) * [u_r, u_i]
    end
end


@named load = PQLoad(Pset=-Pload_pu, Qset=0)
@named shunt = DynamicShunt(C=Cline_pu, ω0=ω0)
loadbus = Bus(
    MTKBus(load, shunt);
    vidx=2
)

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
        Dt(i_r) ~  ω0 * i_i  - R/L * i_r + 1/L*(dst.u_r - src.u_r)
        Dt(i_i) ~ -ω0 * i_r  - R/L * i_i + 1/L*(dst.u_i - src.u_i)
        ## connection to terminal
        src.i_r ~ -i_r
        src.i_i ~ -i_i
        dst.i_r ~ i_r
        dst.i_i ~ i_i
    end
end
@named branch = DynamicRLBranch(; R=Rline_pu, L=Lline_pu, ω0=ω0)
line_model = Line(
    MTKLine(branch);
    src=1, dst=2
)

nw = Network([slackbus, loadbus], line_model)
try #hide #md
s0 = find_fixpoint(nw; alg=DynamicSS(Rodas5P()))
catch e #hide #md
    @error e #hide #md
end #hide #md

#=
well thats a pity, initialiation of those systems is not os easy.
This time we have to reach reach deep to finde a model suitable for initialization
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
        Dt(i_r) ~ 1e3*(Pset * terminal.u_r/(terminal.u_r^2 + terminal.u_i^2) - i_r)
        Dt(i_i) ~ 1e3*(Pset * terminal.u_i/(terminal.u_r^2 + terminal.u_i^2) - i_i)
        terminal.i_r ~ i_r
        terminal.i_i ~ i_i
    end
end
@named less_stiff_load = LessStiffPQLoad(Pset=-Pload_pu)
less_stiff_loadbus = Bus(
    MTKBus(less_stiff_load, shunt);
    vidx=2
)
less_stiff_nw = Network([slackbus, less_stiff_loadbus], line_model)
less_stiff_s0 = find_fixpoint(less_stiff_nw; alg=DynamicSS(Rodas5P()))

#=
perfect the trick with the less stiff load worked!
No lets us this a a starting point for the system we actually want to solve
=#

s0guess = NWState(nw)
s0guess[VIndex(2, :busbar₊u_i)] = less_stiff_s0[VIndex(2, :busbar₊u_i)]
s0guess[VIndex(2, :busbar₊u_r)] = less_stiff_s0[VIndex(2, :busbar₊u_r)]
s0guess[EIndex(1, :branch₊i_i)] = less_stiff_s0[EIndex(1, :branch₊i_i)]
s0guess[EIndex(1, :branch₊i_r)] = less_stiff_s0[EIndex(1, :branch₊i_r)]
s0 = find_fixpoint(nw, s0guess; alg=DynamicSS(Rodas5P()))

#=
yay, workd

no for the perturbation, diable load at 0.1s
=#
disable_load_affect = ComponentAffect([], [:load₊Pset]) do u, p, ctx
    println("Disabling load affect at time $(ctx.t)")
    p[:load₊Pset] = 0
end
set_callback!(loadbus, PresetTimeComponentCallback(0.1, disable_load_affect))
loadbus #hide #md

#=
now we can simualte
=#
prob = ODEProblem(nw, uflat(s0), (0.0, 0.124), copy(pflat(s0)); callback=get_callbacks(nw))
sol = solve(prob, Rodas5P());


fig = let
    fig = Figure()
    ax = Axis(fig[1,1])
    ts = range(0.09, 0.124; length=2000)

    df = CSV.read(
        joinpath(pkgdir(PowerDynamics),"docs","examples", "emt_data_minimal.csv.gz"),
        DataFrame
    )
    lines!(df.t, df.u_2_a; label="PowerFactory A", color=:lightgray,  linewidth=5)
    lines!(df.t, df.u_2_b; label="PowerFactory B", color=:lightgray,  linewidth=5)
    lines!(df.t, df.u_2_c; label="PowerFactory C", color=:lightgray,  linewidth=5)

    a = sol(ts, idxs=VIndex(2, :shunt₊u_a)).u
    b = sol(ts, idxs=VIndex(2, :shunt₊u_b)).u
    c = sol(ts, idxs=VIndex(2, :shunt₊u_c)).u
    lines!(ts, a, label="a phase")
    lines!(ts, b, label="a phase")
    lines!(ts, c, label="a phase")
    xlims!(ax, ts[begin], ts[end])
    fig
end

#=
Lets zoom in for comparison
=#

xlims!(0.0995,0.105)
fig
