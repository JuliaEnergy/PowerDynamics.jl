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
    joinpath(pkgdir(PowerDynamics),"test","OpenIPSL_test","IEEEST","modelica_results.csv.gz"),
    DataFrame;
    drop=(i,name) -> contains(string(name), "nrows="),
    silencewarnings=true
)

# GENROE+ESST1A+IEEEST bus model parameters from OpenIPSL IEEEST test case
BUS = let
    # Use default GENROE parameters (matches OpenIPSL IEEEST test exactly)
    GENROE = PowerDynamicsTesting.default_controller_smib_genroe()

    # Create ESST1A exciter with parameters from OpenIPSL IEEEST test case
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
        vothsg_input=true,
    )

    # Create IEEEST PSS with parameters from OpenIPSL IEEEST.mo (lines 48-66)
    @named ieeest = PSSE_IEEEST(;
        A1=48.7435,
        A2=4.7488,
        A3=0.0,
        A4=0.0,
        A5=-85.7761,
        A6=0.0459,
        T1=0.7361,
        T2=1.5868,
        T3=0.0,
        T4=0.02,
        T5=13.8921,
        T6=0.1057,
        Ks=0.0099,
        Lsmax=0.1,
        Lsmin=-0.1,
        # deviation from OpenIPSL: they deactivate on 0/0
        Vcu=Inf,
        Vcl=-Inf,
    )

    con = [
        # Connect exciter output to machine field voltage
        connect(esst1a.EFD_out, GENROE.EFD_in)

        # Connect machine field current to exciter
        connect(GENROE.XADIFD_out, esst1a.XADIFD_in)

        # Connect terminal voltage measurements
        connect(GENROE.ETERM_out, esst1a.ECOMP_in)
        connect(GENROE.ETERM_out, esst1a.VT_in)
        connect(GENROE.ETERM_out, ieeest.V_CT_in)

        # Connect PSS input (electrical power) and output
        connect(GENROE.PELEC_out, ieeest.INPUT_in)
        connect(ieeest.VOTHSG_out, esst1a.VOTHSG_in)
    ]

    angle_0 = 0.070492225331847
    v_0 = 1
    # Create bus model with proper connections
    busmodel = MTKBus([GENROE, esst1a, ieeest], con; name=:GEN1)
    bm = compile_bus(busmodel, pf=pfSlack(V=v_0, δ=angle_0))
end

sol = OpenIPSL_SMIB(BUS);

## Validation tests for GENROE machine variables (core 5 variables)
@test ref_rms_error(sol, ref, VIndex(:GEN1, :machine₊w), "gENROE.w") < 2e-6
@test ref_rms_error(sol, ref, VIndex(:GEN1, :machine₊delta), "gENROE.delta") < 7e-4
@test ref_rms_error(sol, ref, VIndex(:GEN1, :machine₊P), "gENROE.P") < 7e-4
@test ref_rms_error(sol, ref, VIndex(:GEN1, :machine₊Q), "gENROE.Q") < 1e-4
@test ref_rms_error(sol, ref, VIndex(:GEN1, :machine₊Vt), "gENROE.Vt") < 5e-5

## Validation tests for ESST1A exciter variables (signal flow order)
@test ref_rms_error(sol, ref, VIndex(:GEN1, :esst1a₊leadlag1₊out), "eSST1A.imLeadLag.y") < 7e-5
@test ref_rms_error(sol, ref, VIndex(:GEN1, :esst1a₊leadlag2₊out), "eSST1A.imLeadLag1.y") < 7e-5
@test ref_rms_error(sol, ref, VIndex(:GEN1, :esst1a₊amplifier₊out), "eSST1A.simpleLagLim.y") < 4e-3
@test ref_rms_error(sol, ref, VIndex(:GEN1, :esst1a₊derivative_feedback₊out), "eSST1A.imDerivativeLag.y") < 5e-4
@test ref_rms_error(sol, ref, VIndex(:GEN1, :esst1a₊EFD_out₊u), "eSST1A.EFD") < 7e-4

## Validation tests for IEEEST PSS variables (requested outputs)
# Filter2 output (input to lead-lag 1)
@test ref_rms_error(sol, ref, VIndex(:GEN1, :ieeest₊filter3₊in), "iEEEST.T_1_T_2.u") < 5e-4
# Filter3 output (lead-lag 1)
@test ref_rms_error(sol, ref, VIndex(:GEN1, :ieeest₊filter3₊out), "iEEEST.T_1_T_2.y") < 5e-4
# Filter4 output (lead-lag 2, input to derivative)
@test ref_rms_error(sol, ref, VIndex(:GEN1, :ieeest₊filter4₊out), "iEEEST.T_3_T_4.y") < 5e-4
# Overall PSS output
@test ref_rms_error(sol, ref, VIndex(:GEN1, :ieeest₊VOTHSG_out₊u), "iEEEST.VOTHSG") < 5e-5

#=
# Create comparison plot for IEEEST PSS validation
fig_ieeest = let
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

    # Plot 3: IEEEST Filter2 Output (Input to Lead-Lag 1)
    ax3 = Axis(fig[2,1]; xlabel="Time [s]", ylabel="[pu]", title="IEEEST: Filter2 Output (T_1_T_2.u)")
    lines!(ax3, ref.time, ref[!, Symbol("iEEEST.T_1_T_2.u")]; label="OpenIPSL", color=Cycled(2), linewidth=2, alpha=0.7)
    lines!(ax3, ts, sol(ts, idxs=VIndex(:GEN1, :ieeest₊filter3₊in)).u; label="PowerDynamics", color=Cycled(2), linewidth=2, linestyle=:dash)
    axislegend(ax3)

    # Plot 4: IEEEST Filter3 Output (Lead-Lag 1)
    ax4 = Axis(fig[2,2]; xlabel="Time [s]", ylabel="[pu]", title="IEEEST: Filter3 Output (T_1_T_2.y)")
    lines!(ax4, ref.time, ref[!, Symbol("iEEEST.T_1_T_2.y")]; label="OpenIPSL", color=Cycled(3), linewidth=2, alpha=0.7)
    lines!(ax4, ts, sol(ts, idxs=VIndex(:GEN1, :ieeest₊filter3₊out)).u; label="PowerDynamics", color=Cycled(3), linewidth=2, linestyle=:dash)
    axislegend(ax4)

    # Plot 5: IEEEST Filter4 Output (Lead-Lag 2)
    ax5 = Axis(fig[3,1]; xlabel="Time [s]", ylabel="[pu]", title="IEEEST: Filter4 Output (T_3_T_4.y)")
    lines!(ax5, ref.time, ref[!, Symbol("iEEEST.T_3_T_4.y")]; label="OpenIPSL", color=Cycled(4), linewidth=2, alpha=0.7)
    lines!(ax5, ts, sol(ts, idxs=VIndex(:GEN1, :ieeest₊filter4₊out)).u; label="PowerDynamics", color=Cycled(4), linewidth=2, linestyle=:dash)
    axislegend(ax5)

    # Plot 6: IEEEST Overall PSS Output
    ax6 = Axis(fig[3,2]; xlabel="Time [s]", ylabel="VOTHSG [pu]", title="IEEEST: Overall PSS Output")
    lines!(ax6, ref.time, ref[!, Symbol("iEEEST.VOTHSG")]; label="OpenIPSL", color=Cycled(5), linewidth=2, alpha=0.7)
    lines!(ax6, ts, sol(ts, idxs=VIndex(:GEN1, :ieeest₊VOTHSG_out₊u)).u; label="PowerDynamics", color=Cycled(5), linewidth=2, linestyle=:dash)
    axislegend(ax6)

    # Plot 7: ESST1A Final Field Voltage Output
    ax7 = Axis(fig[4,1]; xlabel="Time [s]", ylabel="EFD [pu]", title="ESST1A: Final Field Voltage")
    lines!(ax7, ref.time, ref[!, Symbol("eSST1A.EFD")]; label="OpenIPSL", color=Cycled(6), linewidth=2, alpha=0.7)
    lines!(ax7, ts, sol(ts, idxs=VIndex(:GEN1, :esst1a₊EFD_out₊u)).u; label="PowerDynamics", color=Cycled(6), linewidth=2, linestyle=:dash)
    axislegend(ax7)

    # Plot 8: Generator Speed Deviation
    ax8 = Axis(fig[4,2]; xlabel="Time [s]", ylabel="ω [pu]", title="Generator: Speed Deviation")
    lines!(ax8, ref.time, ref[!, Symbol("gENROE.w")]; label="OpenIPSL", color=Cycled(7), linewidth=2, alpha=0.7)
    lines!(ax8, ts, sol(ts, idxs=VIndex(:GEN1, :machine₊w)).u; label="PowerDynamics", color=Cycled(7), linewidth=2, linestyle=:dash)
    axislegend(ax8)

    fig
end
=#
