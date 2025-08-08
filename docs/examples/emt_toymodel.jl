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


@mtkmodel DynamicRLBranch begin
    @components begin
        src = Terminal()
        dst = Terminal()
    end
    @variables begin
        i_r(t)=0, [description="Current in real part"]
        i_i(t)=0, [description="Current in imaginary part"]
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

@mtkmodel DynamicPQLoad begin
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

@named branch = DynamicRLBranch(; R=Rline_pu, L=Lline_pu, ω0=ω0)
line_model = Line(
    MTKLine(branch);
    src=1, dst=2
)

# @named load = DynamicPQLoad(Pset=-Pload_pu)
@named load = PQLoad(Pset=-Pload_pu, Qset=0)
@named shunt = DynamicShunt(C=Cline_pu, ω0=ω0)
loadbus = Bus(
    MTKBus(load, shunt);
    vidx=2
)
disable_load_affect = ComponentAffect([], [:load₊Pset]) do u, p, ctx
    println("Disabling load affect at time $(ctx.t)")
    p[:load₊Pset] = 0
end
set_callback!(loadbus, PresetTimeComponentCallback(0.1, disable_load_affect))
loadbus #hide #md


slackbus = Bus(pfSlack(; V=1), vidx=1)

nw = Network([slackbus, loadbus], line_model)
s0guess = NWState(nw)
# for +P
# s0guess[VIndex(2, :busbar₊u_i)] = 0.02413382124823399
# s0guess[VIndex(2, :busbar₊u_r)] = 1.0242809498380367
# s0guess[EIndex(1, :branch₊i_i)] = -0.002967095041392488
# s0guess[EIndex(1, :branch₊i_r)] = 0.9763645487283649

  s0guess[VIndex(2, :busbar₊u_i)] = -0.025390487676966115
  s0guess[VIndex(2, :busbar₊u_r)] = 0.9745092559815217
  s0guess[EIndex(1, :branch₊i_i)] = 0.00202183746961876
  s0guess[EIndex(1, :branch₊i_r)] = -1.026104840439586

s0 = find_fixpoint(nw, s0guess; alg=DynamicSS(Rodas5P()))


prob = ODEProblem(nw, uflat(s0), (0.0, 0.124), copy(pflat(s0)); callback=get_callbacks(nw))
sol = solve(prob, Rodas5P());

df = CSV.read(joinpath(pkgdir(PowerDynamics),"docs","examples","emt_data", "Test_EMT.csv"), DataFrame, skipto=3,
              header=[:t, :u_1_a, :u_1_b, :u_1_c, :u_2_a, :u_2_b, :u_2_c])

let
    fig = Figure()
    ax = Axis(fig[1,1])
    ts = range(0.09, 0.124; length=2000)
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
end |> display
