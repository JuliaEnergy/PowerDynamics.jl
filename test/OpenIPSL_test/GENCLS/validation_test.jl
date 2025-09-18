using PowerDynamics
PowerDynamics.load_pdtesting()
using Main.PowerDynamicsTesting

using PowerDynamics.Library
using ModelingToolkit
using OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqNonlinearSolve


using CSV
using DataFrames
using CairoMakie


ref = CSV.read(
    joinpath(pkgdir(PowerDynamics),"test","OpenIPSL_test","GENCLS","modelica_results.csv.gz"),
    DataFrame
)

# bus 1 is provided from outside
# TEMP: for now, hardcoded for GENCLS testing
GENCLS_BUS = let
    # GENCLS generator parameters from OpenIPSL test
    S_b = 100e6
    H = 6.0
    M_b = 100e6
    P_0 = 40e6
    Q_0 = 5.416582e6
    v_0 = 1.0
    angle_0 = 0.070492225331847
    X_d = 0.2
    D = 0
    V_b = 400e3
    ω_b = 2π*50

    @named gencls = PSSE_GENCLS(; S_b, V_b, ω_b, H, M_b, P_0, Q_0, v_0, angle_0, X_d, D)
    busmodel = MTKBus(gencls; name=:GEN1)
    # compile_bus(busmodel, pf=pfSlack(V=v_0, δ=angle_0))
    compile_bus(busmodel, pf=pfPQ(P=P_0/S_b, Q=Q_0/S_b))
    # compile_bus(busmodel)
end

sol = OpenIPSL_SMIB(GENCLS_BUS)

fig = let
    fig = Figure(size=(1200, 800))
    ts = sol.t

    # Plot 1: Angular frequency ω
    ax1 = Axis(fig[1,1]; xlabel="Time [s]", ylabel="ω [pu]", title="Angular Frequency")
    lines!(ax1, ts, sol(ts, idxs=VIndex(1, :gencls₊ω)).u; label="PD", color=:blue)
    lines!(ax1, ref.time, ref."gENCLS1.omega"; label="OpenIPSL", color=:red, linestyle=:dash)
    axislegend(ax1)

    # Plot 2: Rotor angle δ
    ax2 = Axis(fig[1,2]; xlabel="Time [s]", ylabel="δ [rad]", title="Rotor Angle")
    lines!(ax2, ts, sol(ts, idxs=VIndex(1, :gencls₊δ)).u; label="PD", color=:blue)
    lines!(ax2, ref.time, ref."gENCLS1.delta"; label="OpenIPSL", color=:red, linestyle=:dash)
    axislegend(ax2)

    # Plot 3: Active power P
    ax3 = Axis(fig[2,1]; xlabel="Time [s]", ylabel="P [pu]", title="Active Power")
    lines!(ax3, ts, sol(ts, idxs=VIndex(1, :gencls₊P)).u; label="PD", color=:blue)
    lines!(ax3, ref.time, ref."gENCLS1.P"; label="OpenIPSL", color=:red, linestyle=:dash)
    axislegend(ax3)

    # Plot 4: Reactive power Q
    ax4 = Axis(fig[2,2]; xlabel="Time [s]", ylabel="Q [pu]", title="Reactive Power")
    lines!(ax4, ts, sol(ts, idxs=VIndex(1, :gencls₊Q)).u; label="PD", color=:blue)
    lines!(ax4, ref.time, ref."gENCLS1.Q"; label="OpenIPSL", color=:red, linestyle=:dash)
    axislegend(ax4)

    # Plot 5: Voltage magnitude V
    ax5 = Axis(fig[3,1:2]; xlabel="Time [s]", ylabel="V [pu]", title="Voltage Magnitude")
    lines!(ax5, ts, sol(ts, idxs=VIndex(1, :gencls₊V)).u; label="PD", color=:blue)
    lines!(ax5, ref.time, ref."gENCLS1.V"; label="OpenIPSL", color=:red, linestyle=:dash)
    axislegend(ax5)

    fig
end
