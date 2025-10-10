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
    joinpath(pkgdir(PowerDynamics),"test","OpenIPSL_test","GENROU","modelica_results.csv.gz"),
    DataFrame;
    drop=(i,name) -> contains(string(name), "nrows="),
    silencewarnings=true
)

# bus 1 is provided from outside
# GENROU generator parameters from OpenIPSL test
GENROU_BUS = let
    # Machine parameters from OpenIPSL GENROU test case
    S_b = 100e6
    M_b = 100e6
    H = 4.28
    D = 0
    V_b = 400e3
    ω_b = 2π*50

    # Machine reactances and time constants
    Xd = 1.84
    Xq = 1.75
    Xpd = 0.41
    Xpq = 0.6
    Xppd = 0.2
    Xppq = 0.2
    Xl = 0.12
    Tpd0 = 5.0
    Tpq0 = 0.9
    Tppd0 = 0.07
    Tppq0 = 0.09

    # Saturation parameters
    S10 = 0.11
    S12 = 0.39

    # Resistance
    R_a = 0

    # Powerflow results (same as GENCLS test)
    v_0 = 1.0
    angle_0 = 0.070492225331847

    @named genrou = PSSE_GENROU(;
        S_b, H, M_b, D,
        Xd, Xq, Xpd, Xpq, Xppd, Xppq, Xl,
        Tpd0, Tpq0, Tppd0, Tppq0,
        S10, S12, R_a,
        pmech_input=false, efd_input=false,
        # delta = ref."gENROU.delta"[1]
    )
    busmodel = MTKBus(genrou; name=:GEN1)
    compile_bus(busmodel, pf=pfSlack(V=v_0, δ=angle_0))
end

sol = OpenIPSL_SMIB(GENROU_BUS);

## perform tests for all variables of interest
# Core machine variables
@test ref_rms_error(sol, ref, VIndex(:GEN1, :genrou₊w), "gENROU.w") < 1e-5
@test ref_rms_error(sol, ref, VIndex(:GEN1, :genrou₊delta), "gENROU.delta") < 1e-3
@test ref_rms_error(sol, ref, VIndex(:GEN1, :genrou₊P), "gENROU.P") < 1e-3
@test ref_rms_error(sol, ref, VIndex(:GEN1, :genrou₊Q), "gENROU.Q") < 1e-4
@test ref_rms_error(sol, ref, VIndex(:GEN1, :genrou₊Vt), "gENROU.Vt") < 1e-4

# State variables
@test ref_rms_error(sol, ref, VIndex(:GEN1, :genrou₊Epd), "gENROU.Epd") < 1e-3
@test ref_rms_error(sol, ref, VIndex(:GEN1, :genrou₊Epq), "gENROU.Epq") < 1e-4
@test ref_rms_error(sol, ref, VIndex(:GEN1, :genrou₊PSIkd), "gENROU.PSIkd") < 1e-3
@test ref_rms_error(sol, ref, VIndex(:GEN1, :genrou₊PSIkq), "gENROU.PSIkq") < 1e-3

# Field and torque
@test ref_rms_error(sol, ref, VIndex(:GEN1, :genrou₊XadIfd), "gENROU.XadIfd") < 1e-3
@test ref_rms_error(sol, ref, VIndex(:GEN1, :genrou₊Te), "gENROU.Te") < 1e-3

# Current and voltage components
@test ref_rms_error(sol, ref, VIndex(:GEN1, :genrou₊id), "gENROU.id") < 1e-3
@test ref_rms_error(sol, ref, VIndex(:GEN1, :genrou₊iq), "gENROU.iq") < 1e-3
@test ref_rms_error(sol, ref, VIndex(:GEN1, :genrou₊ud), "gENROU.ud") < 1e-3
@test ref_rms_error(sol, ref, VIndex(:GEN1, :genrou₊uq), "gENROU.uq") < 1e-3

# debug code below
if isdefined(Main, :EXPORT_FIGURES) && Main.EXPORT_FIGURES
    # Create comprehensive comparison plot
    fig = let
        fig = Figure(size=(1400, 1000))
        ts = refine_timeseries(sol.t)

        # Plot 1: Angular frequency ω
        ax1 = Axis(fig[1,1]; xlabel="Time [s]", ylabel="ω [pu]", title="Generator: Angular Frequency")
        lines!(ax1, ref.time, ref[!, Symbol("gENROU.w")]; label="OpenIPSL", color=Cycled(1), linewidth=2, alpha=0.7)
        lines!(ax1, ts, sol(ts, idxs=VIndex(:GEN1, :genrou₊w)).u; label="PowerDynamics", color=Cycled(1), linewidth=2, linestyle=:dash)
        axislegend(ax1)

        # Plot 2: Rotor angle δ
        ax2 = Axis(fig[1,2]; xlabel="Time [s]", ylabel="δ [rad]", title="Generator: Rotor Angle")
        lines!(ax2, ref.time, ref[!, Symbol("gENROU.delta")]; label="OpenIPSL", color=Cycled(1), linewidth=2, alpha=0.7)
        lines!(ax2, ts, sol(ts, idxs=VIndex(:GEN1, :genrou₊delta)).u; label="PowerDynamics", color=Cycled(1), linewidth=2, linestyle=:dash)
        axislegend(ax2)

        # Plot 3: Active power P
        ax3 = Axis(fig[2,1]; xlabel="Time [s]", ylabel="P [pu]", title="Generator: Active Power")
        lines!(ax3, ref.time, ref[!, Symbol("gENROU.P")]; label="OpenIPSL", color=Cycled(2), linewidth=2, alpha=0.7)
        lines!(ax3, ts, sol(ts, idxs=VIndex(:GEN1, :genrou₊P)).u; label="PowerDynamics", color=Cycled(2), linewidth=2, linestyle=:dash)
        axislegend(ax3)

        # Plot 4: Reactive power Q
        ax4 = Axis(fig[2,2]; xlabel="Time [s]", ylabel="Q [pu]", title="Generator: Reactive Power")
        lines!(ax4, ref.time, ref[!, Symbol("gENROU.Q")]; label="OpenIPSL", color=Cycled(2), linewidth=2, alpha=0.7)
        lines!(ax4, ts, sol(ts, idxs=VIndex(:GEN1, :genrou₊Q)).u; label="PowerDynamics", color=Cycled(2), linewidth=2, linestyle=:dash)
        axislegend(ax4)

        # Plot 5: Terminal voltage Vt
        ax5 = Axis(fig[3,1]; xlabel="Time [s]", ylabel="Vt [pu]", title="Generator: Terminal Voltage")
        lines!(ax5, ref.time, ref[!, Symbol("gENROU.Vt")]; label="OpenIPSL", color=Cycled(3), linewidth=2, alpha=0.7)
        lines!(ax5, ts, sol(ts, idxs=VIndex(:GEN1, :genrou₊Vt)).u; label="PowerDynamics", color=Cycled(3), linewidth=2, linestyle=:dash)
        axislegend(ax5)

        # Plot 6: Electrical torque Te
        ax6 = Axis(fig[3,2]; xlabel="Time [s]", ylabel="Te [pu]", title="Generator: Electrical Torque")
        lines!(ax6, ref.time, ref[!, Symbol("gENROU.Te")]; label="OpenIPSL", color=Cycled(3), linewidth=2, alpha=0.7)
        lines!(ax6, ts, sol(ts, idxs=VIndex(:GEN1, :genrou₊Te)).u; label="PowerDynamics", color=Cycled(3), linewidth=2, linestyle=:dash)
        axislegend(ax6)

        # Plot 7: State variables Epd and Epq
        ax7 = Axis(fig[4,1]; xlabel="Time [s]", ylabel="[pu]", title="Generator: Transient EMF")
        lines!(ax7, ref.time, ref[!, Symbol("gENROU.Epd")]; label="OpenIPSL Epd", color=Cycled(4), linewidth=2, alpha=0.7)
        lines!(ax7, ts, sol(ts, idxs=VIndex(:GEN1, :genrou₊Epd)).u; label="PowerDynamics Epd", color=Cycled(4), linewidth=2, linestyle=:dash)
        lines!(ax7, ref.time, ref[!, Symbol("gENROU.Epq")]; label="OpenIPSL Epq", color=Cycled(5), linewidth=2, alpha=0.7)
        lines!(ax7, ts, sol(ts, idxs=VIndex(:GEN1, :genrou₊Epq)).u; label="PowerDynamics Epq", color=Cycled(5), linewidth=2, linestyle=:dash)
        axislegend(ax7)

        # Plot 8: Field current XadIfd
        ax8 = Axis(fig[4,2]; xlabel="Time [s]", ylabel="XadIfd [pu]", title="Generator: Field Current")
        lines!(ax8, ref.time, ref[!, Symbol("gENROU.XadIfd")]; label="OpenIPSL", color=Cycled(6), linewidth=2, alpha=0.7)
        lines!(ax8, ts, sol(ts, idxs=VIndex(:GEN1, :genrou₊XadIfd)).u; label="PowerDynamics", color=Cycled(6), linewidth=2, linestyle=:dash)
        axislegend(ax8)

        fig
    end
    save(joinpath(pkgdir(PowerDynamics),"docs","src","assets","OpenIPSL_valid","GENROU.png"), fig)
end

#=
# debug code below
# Load extended reference data for bus and power comparisons
ref_extended = CSV.read(
    joinpath(pkgdir(PowerDynamics),"test","OpenIPSL_test","GENROU","modelica_results_extended.csv.gz"),
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

    # GENROU machine (main generator)
    lines!(ax, ref.time, ref[!, "gENROU.P"]; label="OpenIPSL: GENROU", color=Cycled(1), linewidth=2, alpha=0.5)
    lines!(ax, ts, sol(ts, idxs=VIndex(:GEN1, :genrou₊P)).u; label="PD: GENROU", linewidth=2, color=Cycled(1), linestyle=:dash)

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

    # GENROU machine (main generator)
    lines!(ax, ref.time, ref[!, "gENROU.Q"]; label="OpenIPSL: GENROU", color=Cycled(1), linewidth=2, alpha=0.5)
    lines!(ax, ts, sol(ts, idxs=VIndex(:GEN1, :genrou₊Q)).u; label="PD: GENROU", linewidth=2, color=Cycled(1), linestyle=:dash)

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
