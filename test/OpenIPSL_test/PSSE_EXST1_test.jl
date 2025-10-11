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

# Load reference data from EXST1 simulation
ref = CSV.read(
    joinpath(pkgdir(PowerDynamics),"test","OpenIPSL_test","EXST1","modelica_results.csv.gz"),
    DataFrame;
    drop=(i,name) -> contains(string(name), "nrows="),
    silencewarnings=true
)

# EXST1 exciter parameters from OpenIPSL EXST1 test case
BUS = let
    GENROE = PowerDynamicsTesting.default_controller_smib_genroe()

    EXST1_EXCITER = let
        # EXST1 exciter parameters (matching OpenIPSL test exactly, lines 27-38)
        T_R = 0.02
        V_IMAX = 10
        V_IMIN = -10
        T_C = 0.1
        T_B = 1
        K_A = 80
        T_A = 0.05
        V_RMAX = 8
        V_RMIN = -3
        K_C = 0.2
        K_F = 0.1
        T_F = 1

        PSSE_EXST1(;
            T_R, V_IMAX, V_IMIN, T_C, T_B,
            K_A, T_A, V_RMAX, V_RMIN, K_C, K_F, T_F,
            name=:ex
        )
    end

    angle_0 = 0.070492225331847
    v_0 = 1
    con = [
        connect(GENROE.ETERM_out, EXST1_EXCITER.ECOMP_in)
        connect(GENROE.XADIFD_out, EXST1_EXCITER.XADIFD_in)
        connect(EXST1_EXCITER.EFD_out, GENROE.EFD_in)
    ]
    busmodel = MTKBus([GENROE, EXST1_EXCITER], con; name=:GEN1)
    compile_bus(busmodel, pf=pfSlack(V=v_0, δ=angle_0))
end

sol = OpenIPSL_SMIB(BUS);

## perform tests for all variables of interest
# Core generator variables
@test ref_rms_error(sol, ref, VIndex(:GEN1, :machine₊w), "gENROE.w") < 1e-5
@test ref_rms_error(sol, ref, VIndex(:GEN1, :machine₊delta), "gENROE.delta") < 8e-4
@test ref_rms_error(sol, ref, VIndex(:GEN1, :machine₊P), "gENROE.P") < 6e-4
@test ref_rms_error(sol, ref, VIndex(:GEN1, :machine₊Q), "gENROE.Q") < 2e-4
@test ref_rms_error(sol, ref, VIndex(:GEN1, :machine₊Vt), "gENROE.Vt") < 8e-5

# EXST1 exciter variables
@test ref_rms_error(sol, ref, VIndex(:GEN1, :ex₊EFD_out₊u), "eXST1.EFD") < 8e-4
@test ref_rms_error(sol, ref, VIndex(:GEN1, :ex₊transducer₊out), "eXST1.TransducerDelay.y") < 7e-5
@test ref_rms_error(sol, ref, VIndex(:GEN1, :ex₊leadlag₊out), "eXST1.imLeadLag.y") < 1e-5
@test ref_rms_error(sol, ref, VIndex(:GEN1, :ex₊amplifier₊out), "eXST1.Vm1.y") < 8e-4
@test ref_rms_error(sol, ref, VIndex(:GEN1, :ex₊derivative_feedback₊out), "eXST1.imDerivativeLag.y") < 6e-5

if isdefined(Main, :EXPORT_FIGURES) && Main.EXPORT_FIGURES
    fig = let
        fig = Figure(size=(1400, 1200))
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

        ax3 = Axis(fig[2,1]; xlabel="Time [s]", ylabel="EFD [pu]", title="ESST4B: Final Field Voltage")
        lines!(ax3, ref.time, ref[!, Symbol("eXST1.EFD")]; label="OpenIPSL", color=Cycled(8), linewidth=2, alpha=0.7)
        lines!(ax3, ts, sol(ts, idxs=VIndex(:GEN1, :ex₊EFD_out₊u)).u; label="PowerDynamics", color=Cycled(8), linewidth=2, linestyle=:dash)
        axislegend(ax3)
        fig

        ax4 = Axis(fig[2,2]; xlabel="Time [s]", ylabel="[pu]", title="ESST1: Transducer Output")
        lines!(ax4, ref.time, ref[!, Symbol("eXST1.TransducerDelay.y")]; label="OpenIPSL", color=Cycled(2), linewidth=2, alpha=0.7)
        lines!(ax4, ts, sol(ts, idxs=VIndex(:GEN1, :ex₊transducer₊out)).u; label="PowerDynamics", color=Cycled(2), linewidth=2, linestyle=:dash)
        axislegend(ax4)

        ax5 = Axis(fig[3,1]; xlabel="Time [s]", ylabel="[pu]", title="ESST1: IM Lead-Lag Output")
        lines!(ax5, ref.time, ref[!, Symbol("eXST1.imLeadLag.y")]; label="OpenIPSL", color=Cycled(3), linewidth=2, alpha=0.7)
        lines!(ax5, ts, sol(ts, idxs=VIndex(:GEN1, :ex₊leadlag₊out)).u; label="PowerDynamics", color=Cycled(3), linewidth=2, linestyle=:dash)
        axislegend(ax5)

        ax6 = Axis(fig[3,2]; xlabel="Time [s]", ylabel="[pu]", title="ESST1: Amplifier Output")
        lines!(ax6, ref.time, ref[!, Symbol("eXST1.Vm1.y")]; label="OpenIPSL", color=Cycled(4), linewidth=2, alpha=0.7)
        lines!(ax6, ts, sol(ts, idxs=VIndex(:GEN1, :ex₊amplifier₊out)).u; label="PowerDynamics", color=Cycled(4), linewidth=2, linestyle=:dash)
        axislegend(ax6)

        fig
    end
    save(joinpath(pkgdir(PowerDynamics),"docs","src","assets","OpenIPSL_valid","EXST1.png"), fig)
end
