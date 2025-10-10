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
    joinpath(pkgdir(PowerDynamics),"test","OpenIPSL_test","ESST1A","modelica_results.csv.gz"),
    DataFrame;
    drop=(i,name) -> contains(string(name), "nrows="),
    silencewarnings=true
)

# GENROE+ESST1A bus model parameters from OpenIPSL ESST1A test case
BUS = let
    # Use default GENROE parameters (matches OpenIPSL ESST1A test exactly)
    GENROE = PowerDynamicsTesting.default_controller_smib_genroe()

    # Create ESST1A exciter with parameters from OpenIPSL ESST1A.mo (lines 26-44)
    @named esst1a = PSSE_ESST1A(;
        V_IMAX=0.3,
        V_IMIN=-0.3,
        T_C=2,
        T_B=10,
        T_C1=0.08,
        T_B1=0.083,
        K_A=300,
        V_AMAX=7,
        V_AMIN=-7,
        V_RMAX=5.2,
        V_RMIN=-5.2,
        K_C=0.38,
        K_F=1,
        T_F=1,
        K_LR=1,
        I_LR=0,
        T_A=0.1,
        T_R=0.1,
        VOS=1,
        UEL=1,
    )

    con = [
        # Connect exciter output to machine field voltage
        connect(esst1a.EFD_out, GENROE.EFD_in)

        # Connect machine field current to exciter
        connect(GENROE.XADIFD_out, esst1a.XADIFD_in)

        # Connect terminal voltage measurements
        connect(GENROE.ETERM_out, esst1a.ECOMP_in)
        connect(GENROE.ETERM_out, esst1a.VT_in)
    ]

    angle_0 = 0.070492225331847
    v_0 = 1
    # Create bus model with proper connections
    busmodel = MTKBus([GENROE, esst1a], con; name=:GEN1)
    bm = compile_bus(busmodel, pf=pfSlack(V=v_0, δ=angle_0))
end

sol = OpenIPSL_SMIB(BUS);

## Validation tests for GENROE machine variables (core 5 variables)
@test ref_rms_error(sol, ref, VIndex(:GEN1, :machine₊w), "gENROE.w") < 5e-6    # Actual: 1.66e-6
@test ref_rms_error(sol, ref, VIndex(:GEN1, :machine₊delta), "gENROE.delta") < 1e-3   # Actual: 5.38e-4
@test ref_rms_error(sol, ref, VIndex(:GEN1, :machine₊P), "gENROE.P") < 1e-3    # Actual: 4.34e-4
@test ref_rms_error(sol, ref, VIndex(:GEN1, :machine₊Q), "gENROE.Q") < 1e-4    # Actual: 6.60e-5
@test ref_rms_error(sol, ref, VIndex(:GEN1, :machine₊Vt), "gENROE.Vt") < 5e-5  # Actual: 1.95e-5

## Validation tests for ESST1A exciter variables (signal flow order)
# 1. First lead-lag output
@test ref_rms_error(sol, ref, VIndex(:GEN1, :esst1a₊leadlag1₊out), "eSST1A.imLeadLag.y") < 1e-5     # Actual: 3.71e-6
# 2. Second lead-lag output
@test ref_rms_error(sol, ref, VIndex(:GEN1, :esst1a₊leadlag2₊out), "eSST1A.imLeadLag1.y") < 1e-5    # Actual: 3.71e-6
# 3. Amplifier/voltage regulator output
@test ref_rms_error(sol, ref, VIndex(:GEN1, :esst1a₊amplifier₊out), "eSST1A.simpleLagLim.y") < 2e-3  # Actual: 1.10e-3
# 4. Derivative feedback output
@test ref_rms_error(sol, ref, VIndex(:GEN1, :esst1a₊derivative_feedback₊out), "eSST1A.imDerivativeLag.y") < 5e-5  # Actual: 1.63e-5
# 5. Final field voltage output
@test ref_rms_error(sol, ref, VIndex(:GEN1, :esst1a₊EFD_out₊u), "eSST1A.EFD") < 1e-3                # Actual: 5.53e-4

# Create focused comparison plot following signal flow
if isdefined(Main, :EXPORT_FIGURES) && Main.EXPORT_FIGURES
    fig = let
        fig = Figure(size=(1400, 1000))
        ts = refine_timeseries(sol.t)

        # Plot 1: Generator Terminal Voltage
        ax1 = Axis(fig[1,1]; xlabel="Time [s]", ylabel="Vt [pu]", title="Generator: Terminal Voltage")
        lines!(ax1, ref.time, ref[!, Symbol("gENROE.Vt")]; label="OpenIPSL", color=Cycled(1), linewidth=2, alpha=0.7)
        lines!(ax1, ts, sol(ts, idxs=VIndex(:GEN1, :machine₊Vt)).u; label="PowerDynamics", color=Cycled(1), linewidth=2, linestyle=:dash)
        axislegend(ax1)

        # Plot 2: Generator Rotor Angle
        ax2 = Axis(fig[1,2]; xlabel="Time [s]", ylabel="δ [rad]", title="Generator: Rotor Angle")
        lines!(ax2, ref.time, ref[!, Symbol("gENROE.delta")]; label="OpenIPSL", color=Cycled(1), linewidth=2, alpha=0.7)
        lines!(ax2, ts, sol(ts, idxs=VIndex(:GEN1, :machine₊delta)).u; label="PowerDynamics", color=Cycled(1), linewidth=2, linestyle=:dash)
        axislegend(ax2)

        # Plot 3: ESST1A First Lead-Lag Output
        ax3 = Axis(fig[2,1]; xlabel="Time [s]", ylabel="[pu]", title="ESST1A: First Lead-Lag Output")
        lines!(ax3, ref.time, ref[!, Symbol("eSST1A.imLeadLag.y")]; label="OpenIPSL", color=Cycled(2), linewidth=2, alpha=0.7)
        lines!(ax3, ts, sol(ts, idxs=VIndex(:GEN1, :esst1a₊leadlag1₊out)).u; label="PowerDynamics", color=Cycled(2), linewidth=2, linestyle=:dash)
        axislegend(ax3)

        # Plot 4: ESST1A Second Lead-Lag Output
        ax4 = Axis(fig[2,2]; xlabel="Time [s]", ylabel="[pu]", title="ESST1A: Second Lead-Lag Output")
        lines!(ax4, ref.time, ref[!, Symbol("eSST1A.imLeadLag1.y")]; label="OpenIPSL", color=Cycled(3), linewidth=2, alpha=0.7)
        lines!(ax4, ts, sol(ts, idxs=VIndex(:GEN1, :esst1a₊leadlag2₊out)).u; label="PowerDynamics", color=Cycled(3), linewidth=2, linestyle=:dash)
        axislegend(ax4)

        # Plot 5: ESST1A Amplifier/Voltage Regulator Output
        ax5 = Axis(fig[3,1]; xlabel="Time [s]", ylabel="[pu]", title="ESST1A: Amplifier Output")
        lines!(ax5, ref.time, ref[!, Symbol("eSST1A.simpleLagLim.y")]; label="OpenIPSL", color=Cycled(4), linewidth=2, alpha=0.7)
        lines!(ax5, ts, sol(ts, idxs=VIndex(:GEN1, :esst1a₊amplifier₊out)).u; label="PowerDynamics", color=Cycled(4), linewidth=2, linestyle=:dash)
        axislegend(ax5)

        # Plot 6: ESST1A Derivative Feedback Output
        ax6 = Axis(fig[3,2]; xlabel="Time [s]", ylabel="[pu]", title="ESST1A: Derivative Feedback")
        lines!(ax6, ref.time, ref[!, Symbol("eSST1A.imDerivativeLag.y")]; label="OpenIPSL", color=Cycled(5), linewidth=2, alpha=0.7)
        lines!(ax6, ts, sol(ts, idxs=VIndex(:GEN1, :esst1a₊derivative_feedback₊out)).u; label="PowerDynamics", color=Cycled(5), linewidth=2, linestyle=:dash)
        axislegend(ax6)

        # Plot 7: ESST1A Final Field Voltage Output
        ax7 = Axis(fig[4,1]; xlabel="Time [s]", ylabel="EFD [pu]", title="ESST1A: Final Field Voltage")
        lines!(ax7, ref.time, ref[!, Symbol("eSST1A.EFD")]; label="OpenIPSL", color=Cycled(6), linewidth=2, alpha=0.7)
        lines!(ax7, ts, sol(ts, idxs=VIndex(:GEN1, :esst1a₊EFD_out₊u)).u; label="PowerDynamics", color=Cycled(6), linewidth=2, linestyle=:dash)
        axislegend(ax7)

        fig
    end
    save(joinpath(pkgdir(PowerDynamics),"docs","src","assets","OpenIPSL_valid","ESST1A.png"), fig)
end
