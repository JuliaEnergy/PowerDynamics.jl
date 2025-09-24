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

# Load reference data from IEEET1 simulation
ref = CSV.read(
    joinpath(pkgdir(PowerDynamics),"test","OpenIPSL_test","IEEET1","modelica_results.csv.gz"),
    DataFrame;
    drop=(i,name) -> contains(string(name), "nrows="),
    silencewarnings=true
)

# IEEET1 exciter parameters from OpenIPSL IEEET1 test case
BUS = let
    GENROE = PowerDynamicsTesting.default_controller_smib_genroe()

    IEEET1_EXCITER = let
        # IEEET1 exciter parameters (matching OpenIPSL test exactly, lines 27-39)
        T_R = 0.02
        K_A = 200
        T_A = 0.001
        T_E = 0.55
        K_F = 0.06
        E_1 = 2.85
        S_EE_1 = 0.3
        E_2 = 3.8
        S_EE_2 = 0.6
        V_RMAX = 2
        V_RMIN = -2
        K_E = 0.1

        PSSE_IEEET1(;
            T_R, K_A, T_A, V_RMAX, V_RMIN,
            K_E, T_E, K_F, T_F=1,  # T_F default from implementation
            E_1, S_EE_1, E_2, S_EE_2,
            name=:ex
        )
    end

    angle_0 = 0.070492225331847
    v_0 = 1
    con = [
        connect(GENROE.ETERM_out, IEEET1_EXCITER.ECOMP_in)
        connect(IEEET1_EXCITER.EFD_out, GENROE.EFD_in)
    ]
    busmodel = MTKBus([GENROE, IEEET1_EXCITER], con; name=:GEN1)
    compile_bus(busmodel, pf=pfSlack(V=v_0, δ=angle_0))
end

sol = OpenIPSL_SMIB(BUS);

## perform tests for all variables of interest
# Core generator variables
@test ref_rms_error(sol, ref, VIndex(:GEN1, :machine₊w), "gENROE.w") < 1.5e-5
@test ref_rms_error(sol, ref, VIndex(:GEN1, :machine₊delta), "gENROE.delta") < 8e-4
@test ref_rms_error(sol, ref, VIndex(:GEN1, :machine₊P), "gENROE.P") < 6e-4
@test ref_rms_error(sol, ref, VIndex(:GEN1, :machine₊Vt), "gENROE.Vt") < 8e-5
@test ref_rms_error(sol, ref, VIndex(:GEN1, :machine₊Epd), "gENROE.Epd") < 2.5e-4
@test ref_rms_error(sol, ref, VIndex(:GEN1, :machine₊Epq), "gENROE.Epq") < 8e-5
@test ref_rms_error(sol, ref, VIndex(:GEN1, :machine₊XadIfd), "gENROE.XadIfd") < 6.5e-4

# IEEET1 exciter variables
@test ref_rms_error(sol, ref, VIndex(:GEN1, :ex₊exciter₊EFD), "iEEET1.EFD") < 1e-3
@test ref_rms_error(sol, ref, VIndex(:GEN1, :ex₊amplifier₊out), "iEEET1.simpleLagLim.y") < 1.5e-3
@test ref_rms_error(sol, ref, VIndex(:GEN1, :ex₊derivative_lag₊out), "iEEET1.derivativeLag.y") < 2.5e-5
@test ref_rms_error(sol, ref, VIndex(:GEN1, :ex₊transducer₊out), "iEEET1.TransducerDelay.y") < 5.5e-5

#=
fig_ieeet1 = let
    fig = Figure(size=(1400, 1000))
    ts = refine_timeseries(sol.t)

    # Plot 1: Angular frequency ω
    ax1 = Axis(fig[1,1]; xlabel="Time [s]", ylabel="ω [pu]", title="Generator: Angular Frequency")
    lines!(ax1, ref.time, ref[!, Symbol("gENROE.w")]; label="OpenIPSL", color=Cycled(1), linewidth=2, alpha=0.7)
    lines!(ax1, ts, sol(ts, idxs=VIndex(:GEN1, :machine₊w)).u; label="PowerDynamics", color=Cycled(1), linewidth=2, linestyle=:dash)
    axislegend(ax1)

    # Plot 2: Rotor angle δ
    ax2 = Axis(fig[1,2]; xlabel="Time [s]", ylabel="δ [rad]", title="Generator: Rotor Angle")
    lines!(ax2, ref.time, ref[!, Symbol("gENROE.delta")]; label="OpenIPSL", color=Cycled(1), linewidth=2, alpha=0.7)
    lines!(ax2, ts, sol(ts, idxs=VIndex(:GEN1, :machine₊delta)).u; label="PowerDynamics", color=Cycled(1), linewidth=2, linestyle=:dash)
    axislegend(ax2)

    # Plot 3: Active power P
    ax3 = Axis(fig[2,1]; xlabel="Time [s]", ylabel="P [pu]", title="Generator: Active Power")
    lines!(ax3, ref.time, ref[!, Symbol("gENROE.P")]; label="OpenIPSL", color=Cycled(2), linewidth=2, alpha=0.7)
    lines!(ax3, ts, sol(ts, idxs=VIndex(:GEN1, :machine₊P)).u; label="PowerDynamics", color=Cycled(2), linewidth=2, linestyle=:dash)
    axislegend(ax3)

    # Plot 4: Terminal voltage Vt
    ax4 = Axis(fig[2,2]; xlabel="Time [s]", ylabel="Vt [pu]", title="Generator: Terminal Voltage")
    lines!(ax4, ref.time, ref[!, Symbol("gENROE.Vt")]; label="OpenIPSL", color=Cycled(3), linewidth=2, alpha=0.7)
    lines!(ax4, ts, sol(ts, idxs=VIndex(:GEN1, :machine₊Vt)).u; label="PowerDynamics", color=Cycled(3), linewidth=2, linestyle=:dash)
    axislegend(ax4)

    # Plot 5: Transient EMF (Epd and Epq)
    ax5 = Axis(fig[3,1]; xlabel="Time [s]", ylabel="[pu]", title="Generator: Transient EMF")
    lines!(ax5, ref.time, ref[!, Symbol("gENROE.Epd")]; label="OpenIPSL Epd", color=Cycled(4), linewidth=2, alpha=0.7)
    lines!(ax5, ts, sol(ts, idxs=VIndex(:GEN1, :machine₊Epd)).u; label="PowerDynamics Epd", color=Cycled(4), linewidth=2, linestyle=:dash)
    lines!(ax5, ref.time, ref[!, Symbol("gENROE.Epq")]; label="OpenIPSL Epq", color=Cycled(5), linewidth=2, alpha=0.7)
    lines!(ax5, ts, sol(ts, idxs=VIndex(:GEN1, :machine₊Epq)).u; label="PowerDynamics Epq", color=Cycled(5), linewidth=2, linestyle=:dash)
    axislegend(ax5)

    # Plot 6: Field current XadIfd
    ax6 = Axis(fig[3,2]; xlabel="Time [s]", ylabel="XadIfd [pu]", title="Generator: Field Current")
    lines!(ax6, ref.time, ref[!, Symbol("gENROE.XadIfd")]; label="OpenIPSL", color=Cycled(6), linewidth=2, alpha=0.7)
    lines!(ax6, ts, sol(ts, idxs=VIndex(:GEN1, :machine₊XadIfd)).u; label="PowerDynamics", color=Cycled(6), linewidth=2, linestyle=:dash)
    axislegend(ax6)

    # Plot 7: Exciter field voltage EFD
    ax7 = Axis(fig[4,1]; xlabel="Time [s]", ylabel="EFD [pu]", title="Exciter: Field Voltage")
    lines!(ax7, ref.time, ref[!, Symbol("iEEET1.EFD")]; label="OpenIPSL", color=Cycled(7), linewidth=2, alpha=0.7)
    lines!(ax7, ts, sol(ts, idxs=VIndex(:GEN1, :ex₊exciter₊EFD)).u; label="PowerDynamics", color=Cycled(7), linewidth=2, linestyle=:dash)
    axislegend(ax7)

    # Plot 8: Exciter internal signals
    ax8 = Axis(fig[4,2]; xlabel="Time [s]", ylabel="[pu]", title="Exciter: Internal Signals")
    lines!(ax8, ref.time, ref[!, Symbol("iEEET1.simpleLagLim.y")]; label="OpenIPSL Amplifier", color=Cycled(8), linewidth=2, alpha=0.7)
    lines!(ax8, ts, sol(ts, idxs=VIndex(:GEN1, :ex₊amplifier₊out)).u; label="PowerDynamics Amplifier", color=Cycled(8), linewidth=2, linestyle=:dash)
    lines!(ax8, ref.time, ref[!, Symbol("iEEET1.derivativeLag.y")]; label="OpenIPSL Derivative", color=Cycled(9), linewidth=2, alpha=0.7)
    lines!(ax8, ts, sol(ts, idxs=VIndex(:GEN1, :ex₊derivative_lag₊out)).u; label="PowerDynamics Derivative", color=Cycled(9), linewidth=2, linestyle=:dash)
    axislegend(ax8)

    fig
end
=#
