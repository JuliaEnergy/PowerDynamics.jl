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
    joinpath(pkgdir(PowerDynamics),"test","OpenIPSL_test","GENSAE","modelica_results.csv.gz"),
    DataFrame;
    drop=(i,name) -> contains(string(name), "nrows="),
    silencewarnings=true
)

# bus 1 is provided from outside
# GENSAE generator parameters from OpenIPSL test
GENSAE_BUS = let
    # GENSAE generator parameters from OpenIPSL test case
    S_b = 100e6
    H = 4.28
    M_b = 100e6
    D = 0
    # Machine reactances
    Xd = 1.84
    Xq = 1.75
    Xpd = 0.41
    Xppd = 0.2
    Xppq = 0.2
    Xl = 0.12
    # Time constants
    Tpd0 = 5.0
    Tppd0 = 0.07
    Tppq0 = 0.09
    # Saturation
    S10 = 0.11
    S12 = 0.39
    # Armature resistance
    R_a = 0.0
    # V_b = 400e3
    # ω_b = 2π*50

    # powerflow results, used to set up pfmodel
    # P_0 = 40e6
    # Q_0 = 5.416582e6
    v_0 = 1.0
    angle_0 = 0.070492225331847

    @named gensae = PSSE_GENSAE(;
        S_b, H, M_b, D,
        Xd, Xq, Xpd, Xppd, Xppq, Xl,
        Tpd0, Tppd0, Tppq0,
        S10, S12, R_a,
        pmech_input=false,
        efd_input=false,
    )
    busmodel = MTKBus(gensae; name=:GEN1)
    compile_bus(busmodel, pf=pfSlack(V=v_0, δ=angle_0))
end

sol = OpenIPSL_SMIB(GENSAE_BUS);

## perform tests for specified variables
# Machine currents
@test ref_rms_error(sol, ref, VIndex(:GEN1, :gensae₊iq), "gENSAE.iq") < 7e-4
@test ref_rms_error(sol, ref, VIndex(:GEN1, :gensae₊id), "gENSAE.id") < 1.3e-3

# State variables (3 states for salient pole)
@test ref_rms_error(sol, ref, VIndex(:GEN1, :gensae₊PSIkd), "gENSAE.PSIkd") < 5e-4
@test ref_rms_error(sol, ref, VIndex(:GEN1, :gensae₊Epq), "gENSAE.Epq") < 2e-4
@test ref_rms_error(sol, ref, VIndex(:GEN1, :gensae₊PSIppq), "gENSAE.PSIppq") < 1e-3

# Rotor dynamics
@test ref_rms_error(sol, ref, VIndex(:GEN1, :gensae₊delta), "gENSAE.delta") < 1.7e-3
@test ref_rms_error(sol, ref, VIndex(:GEN1, :gensae₊w), "gENSAE.w") < 8e-6

# Observables
@test ref_rms_error(sol, ref, VIndex(:GEN1, :gensae₊P), "gENSAE.P") < 1.3e-3
@test ref_rms_error(sol, ref, VIndex(:GEN1, :gensae₊Q), "gENSAE.Q") < 2e-4
@test ref_rms_error(sol, ref, VIndex(:GEN1, :gensae₊anglev), "gENSAE.anglev") < 5e-4
@test ref_rms_error(sol, ref, VIndex(:GEN1, :gensae₊Vt), "gENSAE.Vt") < 4e-5

# Create comprehensive comparison plot
if isdefined(Main, :EXPORT_FIGURES) && Main.EXPORT_FIGURES
    fig = let
        fig = Figure(size=(1400, 1000))
        ts = refine_timeseries(sol.t)

        # Plot 1: Angular frequency ω
        ax1 = Axis(fig[1,1]; xlabel="Time [s]", ylabel="ω [pu]", title="Generator: Angular Frequency")
        lines!(ax1, ref.time, ref[!, Symbol("gENSAE.w")]; label="OpenIPSL", color=Cycled(1), linewidth=2, alpha=0.7)
        lines!(ax1, ts, sol(ts, idxs=VIndex(:GEN1, :gensae₊w)).u; label="PowerDynamics", color=Cycled(1), linewidth=2, linestyle=:dash)
        axislegend(ax1)

        # Plot 2: Rotor angle δ
        ax2 = Axis(fig[1,2]; xlabel="Time [s]", ylabel="δ [rad]", title="Generator: Rotor Angle")
        lines!(ax2, ref.time, ref[!, Symbol("gENSAE.delta")]; label="OpenIPSL", color=Cycled(1), linewidth=2, alpha=0.7)
        lines!(ax2, ts, sol(ts, idxs=VIndex(:GEN1, :gensae₊delta)).u; label="PowerDynamics", color=Cycled(1), linewidth=2, linestyle=:dash)
        axislegend(ax2)

        # Plot 3: Active power P
        ax3 = Axis(fig[2,1]; xlabel="Time [s]", ylabel="P [pu]", title="Generator: Active Power")
        lines!(ax3, ref.time, ref[!, Symbol("gENSAE.P")]; label="OpenIPSL", color=Cycled(2), linewidth=2, alpha=0.7)
        lines!(ax3, ts, sol(ts, idxs=VIndex(:GEN1, :gensae₊P)).u; label="PowerDynamics", color=Cycled(2), linewidth=2, linestyle=:dash)
        axislegend(ax3)

        # Plot 4: Reactive power Q
        ax4 = Axis(fig[2,2]; xlabel="Time [s]", ylabel="Q [pu]", title="Generator: Reactive Power")
        lines!(ax4, ref.time, ref[!, Symbol("gENSAE.Q")]; label="OpenIPSL", color=Cycled(2), linewidth=2, alpha=0.7)
        lines!(ax4, ts, sol(ts, idxs=VIndex(:GEN1, :gensae₊Q)).u; label="PowerDynamics", color=Cycled(2), linewidth=2, linestyle=:dash)
        axislegend(ax4)

        # Plot 5: Terminal voltage Vt
        ax5 = Axis(fig[3,1]; xlabel="Time [s]", ylabel="Vt [pu]", title="Generator: Terminal Voltage")
        lines!(ax5, ref.time, ref[!, Symbol("gENSAE.Vt")]; label="OpenIPSL", color=Cycled(3), linewidth=2, alpha=0.7)
        lines!(ax5, ts, sol(ts, idxs=VIndex(:GEN1, :gensae₊Vt)).u; label="PowerDynamics", color=Cycled(3), linewidth=2, linestyle=:dash)
        axislegend(ax5)

        # Plot 6: Voltage angle
        ax6 = Axis(fig[3,2]; xlabel="Time [s]", ylabel="anglev [rad]", title="Generator: Voltage Angle")
        lines!(ax6, ref.time, ref[!, Symbol("gENSAE.anglev")]; label="OpenIPSL", color=Cycled(3), linewidth=2, alpha=0.7)
        lines!(ax6, ts, sol(ts, idxs=VIndex(:GEN1, :gensae₊anglev)).u; label="PowerDynamics", color=Cycled(3), linewidth=2, linestyle=:dash)
        axislegend(ax6)

        # Plot 7: State variables PSIkd and Epq
        ax7 = Axis(fig[4,1]; xlabel="Time [s]", ylabel="[pu]", title="Generator: State Variables")
        lines!(ax7, ref.time, ref[!, Symbol("gENSAE.PSIkd")]; label="OpenIPSL PSIkd", color=Cycled(4), linewidth=2, alpha=0.7)
        lines!(ax7, ts, sol(ts, idxs=VIndex(:GEN1, :gensae₊PSIkd)).u; label="PowerDynamics PSIkd", color=Cycled(4), linewidth=2, linestyle=:dash)
        lines!(ax7, ref.time, ref[!, Symbol("gENSAE.Epq")]; label="OpenIPSL Epq", color=Cycled(5), linewidth=2, alpha=0.7)
        lines!(ax7, ts, sol(ts, idxs=VIndex(:GEN1, :gensae₊Epq)).u; label="PowerDynamics Epq", color=Cycled(5), linewidth=2, linestyle=:dash)
        axislegend(ax7)

        # Plot 8: Currents id and iq
        ax8 = Axis(fig[4,2]; xlabel="Time [s]", ylabel="Current [pu]", title="Generator: dq Currents")
        lines!(ax8, ref.time, ref[!, Symbol("gENSAE.id")]; label="OpenIPSL id", color=Cycled(6), linewidth=2, alpha=0.7)
        lines!(ax8, ts, sol(ts, idxs=VIndex(:GEN1, :gensae₊id)).u; label="PowerDynamics id", color=Cycled(6), linewidth=2, linestyle=:dash)
        lines!(ax8, ref.time, ref[!, Symbol("gENSAE.iq")]; label="OpenIPSL iq", color=Cycled(7), linewidth=2, alpha=0.7)
        lines!(ax8, ts, sol(ts, idxs=VIndex(:GEN1, :gensae₊iq)).u; label="PowerDynamics iq", color=Cycled(7), linewidth=2, linestyle=:dash)
        axislegend(ax8)

        fig
    end
    save(joinpath(pkgdir(PowerDynamics),"docs","src","assets","OpenIPSL_valid","GENSAE.png"), fig)
end
