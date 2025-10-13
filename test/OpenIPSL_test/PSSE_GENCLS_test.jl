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
    DataFrame;
    drop=(i,name) -> contains(string(name), "nrows="),
    silencewarnings=true
)

# bus 1 is provided from outside
# TEMP: for now, hardcoded for GENCLS testing
GENCLS_BUS = let
    # GENCLS generator parameters from OpenIPSL test
    S_b = 100e6
    H = 6.0
    M_b = 100e6
    X_d = 0.2
    D = 0
    # V_b = 400e3
    ω_b = 2π*50

    # powerflow results, used to set up pfmodel
    # P_0 = 40e6
    # Q_0 = 5.416582e6
    v_0 = 1.0
    angle_0 = 0.070492225331847

    @named gencls = PSSE_GENCLS(; S_b, ω_b, H, M_b, #=P_0,=# X_d, D)
    busmodel = MTKBus(gencls; name=:GEN1)
    compile_bus(busmodel, pf=pfSlack(V=v_0, δ=angle_0))
end

sol = OpenIPSL_SMIB(GENCLS_BUS);

## perform tests
@test ref_rms_error(sol, ref, VIndex(:GEN1, :gencls₊ω), "gENCLS1.omega") < 1e-5
@test ref_rms_error(sol, ref, VIndex(:GEN1, :gencls₊δ), "gENCLS1.delta") < 1e-3
@test ref_rms_error(sol, ref, VIndex(:GEN1, :gencls₊P), "gENCLS1.P") < 1e-3
@test ref_rms_error(sol, ref, VIndex(:GEN1, :gencls₊Q), "gENCLS1.Q") < 1e-4

# debug code blow
if isdefined(Main, :EXPORT_FIGURES) && Main.EXPORT_FIGURES
    fig = let
        fig = Figure(size=(1200, 800))
        ts = refine_timeseries(sol.t)

        # Plot 1: Angular frequency ω
        ax1 = Axis(fig[1,1]; xlabel="Time [s]", ylabel="ω [pu]", title="Generator: Angular Frequency")
        lines!(ax1, ref.time, ref[!, Symbol("gENCLS1.omega")]; label="OpenIPSL", color=Cycled(1), linewidth=2, alpha=0.5)
        lines!(ax1, ts, sol(ts, idxs=VIndex(1, :gencls₊ω)).u; label="PD", color=Cycled(1), linewidth=2, linestyle=:dash)
        axislegend(ax1)

        # Plot 2: Rotor angle δ
        ax2 = Axis(fig[1,2]; xlabel="Time [s]", ylabel="δ [rad]", title="Generator: Rotor Angle")
        lines!(ax2, ref.time, ref[!, Symbol("gENCLS1.delta")]; label="OpenIPSL", color=Cycled(1), linewidth=2, alpha=0.5)
        lines!(ax2, ts, sol(ts, idxs=VIndex(1, :gencls₊δ)).u; label="PD", color=Cycled(1), linewidth=2, linestyle=:dash)
        axislegend(ax2)

        # Plot 3: Active power P
        ax3 = Axis(fig[2,1]; xlabel="Time [s]", ylabel="P [pu]", title="Generator: Active Power")
        lines!(ax3, ref.time, ref[!, Symbol("gENCLS1.P")]; label="OpenIPSL", color=Cycled(1), linewidth=2, alpha=0.5)
        lines!(ax3, ts, sol(ts, idxs=VIndex(1, :gencls₊P)).u; label="PD", color=Cycled(1), linewidth=2, linestyle=:dash)
        axislegend(ax3)

        # Plot 4: Reactive power Q
        ax4 = Axis(fig[2,2]; xlabel="Time [s]", ylabel="Q [pu]", title="Generator: Reactive Power")
        lines!(ax4, ref.time, ref[!, Symbol("gENCLS1.Q")]; label="OpenIPSL", color=Cycled(1), linewidth=2, alpha=0.5)
        lines!(ax4, ts, sol(ts, idxs=VIndex(1, :gencls₊Q)).u; label="PD", color=Cycled(1), linewidth=2, linestyle=:dash)
        axislegend(ax4)

        # Plot 5: Voltage magnitude V
        ax5 = Axis(fig[3,1:2]; xlabel="Time [s]", ylabel="V [pu]", title="Generator: Voltage Magnitude")
        lines!(ax5, ref.time, ref[!, Symbol("gENCLS1.V")]; label="OpenIPSL", color=Cycled(1), linewidth=2, alpha=0.5)
        lines!(ax5, ts, sol(ts, idxs=VIndex(1, :gencls₊V)).u; label="PD", color=Cycled(1), linewidth=2, linestyle=:dash)
        axislegend(ax5)

        fig
    end
    save(joinpath(pkgdir(PowerDynamics),"docs","src","assets","OpenIPSL_valid","GENCLS.png"), fig)
end

#=
# Load extended reference data for bus and power comparisons
ref_extended = CSV.read(
    joinpath(pkgdir(PowerDynamics),"test","OpenIPSL_test","GENCLS","modelica_results_extended.csv.gz"),
    DataFrame;
    drop=(i,name) -> contains(string(name), "nrows="),
    silencewarnings=true
)
# Separate plot for Voltage Magnitude across all buses
fig_V = let
    fig = Figure(size=(1200, 600))
    ts = refine_timeseries(sol.t)
    ax = Axis(fig[1,1]; xlabel="Time [s]", ylabel="Voltage Magnitude [pu]", title="Voltage Magnitude Comparison - All Buses")

    for (i, name) in enumerate([:GEN1, :GEN2, :LOAD, :FAULT, :SHUNT])
        lines!(ax, ref_extended.time, ref_extended[!, "$name.v"]; label="OpenIPSL: $name Bus", color=Cycled(i), linewidth=2, alpha=0.5)
        lines!(ax, ts, sol(ts, idxs=VIndex(name, :busbar₊u_mag)).u; label="PD: $name Bus", linewidth=2, color=Cycled(i), linestyle=:dash)
    end

    axislegend(ax; position=:rt)
    fig
end

# Separate plot for Voltage Angles across all buses
fig_angles = let
    fig = Figure(size=(1200, 600))
    ts = refine_timeseries(sol.t)
    ax = Axis(fig[1,1]; xlabel="Time [s]", ylabel="Voltage Angle [rad]", title="Voltage Angle Comparison - All Buses")

    for (i, name) in enumerate([:GEN1, :GEN2, :LOAD, #=:FAULT,=# :SHUNT])
        lines!(ax, ref_extended.time, ref_extended[!, "$name.angle"]; label="OpenIPSL: $name Bus θ", color=Cycled(i), linewidth=2, alpha=0.5)
        lines!(ax, ts, sol(ts, idxs=VIndex(name, :busbar₊u_arg)).u; label="PD: $name Bus θ", linewidth=2, color=Cycled(i), linestyle=:dash)
    end

    axislegend(ax; position=:rt)
    fig
end

# Power comparison plots for injectors
fig_P = let
    fig = Figure(size=(1200, 600))
    ts = refine_timeseries(sol.t)
    ax = Axis(fig[1,1]; xlabel="Time [s]", ylabel="Active Power [pu]", title="Active Power Comparison - Injectors")

    # GENCLS1 machine (main generator)
    lines!(ax, ref.time, ref[!, "gENCLS1.P"]; label="OpenIPSL: GENCLS1", color=Cycled(1), linewidth=2, alpha=0.5)
    lines!(ax, ts, sol(ts, idxs=VIndex(:GEN1, :gencls₊P)).u; label="PD: GENCLS1", linewidth=2, color=Cycled(1), linestyle=:dash)

    # Constant Load
    lines!(ax, ref_extended.time, ref_extended[!, "constantLoad.P"]; label="OpenIPSL: Load", color=Cycled(2), linewidth=2, alpha=0.5)
    lines!(ax, ts, sol(ts, idxs=VIndex(:LOAD, :constantLoad₊P)).u; label="PD: Load", linewidth=2, color=Cycled(2), linestyle=:dash)

    # Infinite bus generator (gENCLS)
    lines!(ax, ref_extended.time, ref_extended[!, "gENCLS.P"]; label="OpenIPSL: Inf Bus", color=Cycled(3), linewidth=2, alpha=0.5)
    lines!(ax, ts, sol(ts, idxs=VIndex(:GEN2, :gencls_inf₊P)).u; label="PD: Inf Bus", linewidth=2, color=Cycled(3), linestyle=:dash)

    axislegend(ax; position=:rt)
    fig
end

fig_Q = let
    fig = Figure(size=(1200, 600))
    ts = refine_timeseries(sol.t)
    ax = Axis(fig[1,1]; xlabel="Time [s]", ylabel="Reactive Power [pu]", title="Reactive Power Comparison - Injectors")

    # GENCLS1 machine (main generator)
    lines!(ax, ref.time, ref[!, "gENCLS1.Q"]; label="OpenIPSL: GENCLS1", color=Cycled(1), linewidth=2, alpha=0.5)
    lines!(ax, ts, sol(ts, idxs=VIndex(:GEN1, :gencls₊Q)).u; label="PD: GENCLS1", linewidth=2, color=Cycled(1), linestyle=:dash)

    # Constant Load
    lines!(ax, ref_extended.time, ref_extended[!, "constantLoad.Q"]; label="OpenIPSL: Load", color=Cycled(2), linewidth=2, alpha=0.5)
    lines!(ax, ts, sol(ts, idxs=VIndex(:LOAD, :constantLoad₊Q)).u; label="PD: Load", linewidth=2, color=Cycled(2), linestyle=:dash)

    # Infinite bus generator (gENCLS)
    lines!(ax, ref_extended.time, ref_extended[!, "gENCLS.Q"]; label="OpenIPSL: Inf Bus", color=Cycled(3), linewidth=2, alpha=0.5)
    lines!(ax, ts, sol(ts, idxs=VIndex(:GEN2, :gencls_inf₊Q)).u; label="PD: Inf Bus", linewidth=2, color=Cycled(3), linestyle=:dash)

    axislegend(ax; position=:rt)
    fig
end
=#
