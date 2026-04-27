using PowerDynamics
PowerDynamics.load_pdtesting()
using Main.PowerDynamicsTesting

using PowerDynamics.Library
using ModelingToolkitBase
using OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqNonlinearSolve

using CSV
using DataFrames
using CairoMakie

# Load reference data from SCRX simulation
ref = CSV.read(
    joinpath(pkgdir(PowerDynamics),"test","OpenIPSL_test","SCRX","modelica_results.csv.gz"),
    DataFrame;
    drop=(i,name) -> contains(string(name), "nrows="),
    silencewarnings=true
)

# SCRX test system parameters from OpenIPSL SCRX test case
BUS = let
    # Power flow from OpenIPSL SCRX test
    v_0 = 1.0
    angle_0 = 0.070492225331847

    # Create GENROU machine with EFD input enabled for exciter control
    @named genrou = PSSE_GENROU(;
        # System parameters
        S_b = 100e6,
        M_b = 100e6,
        # GENROU machine parameters from OpenIPSL SCRX.mo (lines 6-28)
        H = 4.28,
        D = 0,
        Xd = 1.84,
        Xq = 1.75,
        Xpd = 0.41,
        Xpq = 0.6,
        Xppd = 0.2,
        Xppq = 0.2,
        Xl = 0.12,
        Tpd0 = 5,
        Tppd0 = 0.07,
        Tpq0 = 0.9,
        Tppq0 = 0.09,
        S10 = 0.11,
        S12 = 0.39,
        R_a = 0,
        pmech_input = false,
        efd_input = true  # EFD controlled by SCRX
    )

    # Create SCRX exciter with parameters matching OpenIPSL test (lines 30-38)
    @named scrx = PSSE_SCRX(;
        T_AT_B = 0.1,
        T_B = 1,
        K = 100,
        T_E = 0.005,
        E_MIN = -10,
        E_MAX = 10,
        r_cr_fd = 0,
        C_SWITCH = false,  # Bus fed mode
        vothsg_input = false,  # No other signal input
        vuel_input = false,    # No under-excitation limiter
        voel_input = false     # No over-excitation limiter
    )

    con = [
        # Connect exciter output to machine field voltage
        connect(scrx.EFD_out, genrou.EFD_in)

        # Connect machine field current to exciter
        connect(genrou.XADIFD_out, scrx.XADIFD_in)

        # Connect terminal voltage measurement
        connect(genrou.ETERM_out, scrx.ECOMP_in)
    ]

    # Create bus model with proper connections
    busmodel = MTKBus([genrou, scrx], con; name=:GEN1)
    bm = compile_bus(busmodel, pf=pfSlack(V=v_0, δ=angle_0), mtkcompile=:compare)
end

sol = OpenIPSL_SMIB(BUS);

## Validation tests for GENROU machine variables (core 3 variables)
@test ref_rms_error(sol, ref, VIndex(:GEN1, :genrou₊w), "gENROU.w") < 1e-5
@test ref_rms_error(sol, ref, VIndex(:GEN1, :genrou₊delta), "gENROU.delta") < 1e-3
@test ref_rms_error(sol, ref, VIndex(:GEN1, :genrou₊Vt), "gENROU.Vt") < 1e-4

## Validation tests for SCRX exciter variables (4 control path variables)
@test ref_rms_error(sol, ref, VIndex(:GEN1, :scrx₊leadlag₊out), "sCRX.imLeadLag.y") < 5e-5
@test ref_rms_error(sol, ref, VIndex(:GEN1, :scrx₊amplifier₊out), "sCRX.simpleLagLim.y") < 5e-3
@test ref_rms_error(sol, ref, VIndex(:GEN1, :genrou₊XADIFD_out₊u), "gENROU.XADIFD") < 1e-3
@test ref_rms_error(sol, ref, VIndex(:GEN1, :scrx₊EFD_out₊u), "sCRX.EFD") < 2e-3

if isdefined(Main, :EXPORT_FIGURES) && Main.EXPORT_FIGURES
    fig = let
        fig = Figure(size=(1400, 1200))
        ts = refine_timeseries(sol.t)

        # Plot 1: Generator Rotor Angle
        ax1 = Axis(fig[1,1]; xlabel="Time [s]", ylabel="δ [rad]", title="Generator: Rotor Angle")
        lines!(ax1, ref.time, ref[!, Symbol("gENROU.delta")]; label="OpenIPSL", color=Cycled(1), linewidth=2, alpha=0.7)
        lines!(ax1, ts, sol(ts, idxs=VIndex(:GEN1, :genrou₊delta)).u; label="PowerDynamics", color=Cycled(1), linewidth=2, linestyle=:dash)
        axislegend(ax1)

        # Plot 2: Generator Speed Deviation
        ax2 = Axis(fig[1,2]; xlabel="Time [s]", ylabel="ω [pu]", title="Generator: Speed Deviation")
        lines!(ax2, ref.time, ref[!, Symbol("gENROU.w")]; label="OpenIPSL", color=Cycled(2), linewidth=2, alpha=0.7)
        lines!(ax2, ts, sol(ts, idxs=VIndex(:GEN1, :genrou₊w)).u; label="PowerDynamics", color=Cycled(2), linewidth=2, linestyle=:dash)
        axislegend(ax2)

        # Plot 3: Generator Terminal Voltage
        ax3 = Axis(fig[2,1]; xlabel="Time [s]", ylabel="Vt [pu]", title="Generator: Terminal Voltage")
        lines!(ax3, ref.time, ref[!, Symbol("gENROU.Vt")]; label="OpenIPSL", color=Cycled(3), linewidth=2, alpha=0.7)
        lines!(ax3, ts, sol(ts, idxs=VIndex(:GEN1, :genrou₊Vt)).u; label="PowerDynamics", color=Cycled(3), linewidth=2, linestyle=:dash)
        axislegend(ax3)

        # Plot 4: SCRX LeadLag Output
        ax4 = Axis(fig[2,2]; xlabel="Time [s]", ylabel="[pu]", title="SCRX: LeadLag Output")
        lines!(ax4, ref.time, ref[!, Symbol("sCRX.imLeadLag.y")]; label="OpenIPSL", color=Cycled(4), linewidth=2, alpha=0.7)
        lines!(ax4, ts, sol(ts, idxs=VIndex(:GEN1, :scrx₊leadlag₊out)).u; label="PowerDynamics", color=Cycled(4), linewidth=2, linestyle=:dash)
        axislegend(ax4)

        # Plot 5: SCRX Amplifier Output
        ax5 = Axis(fig[3,1]; xlabel="Time [s]", ylabel="[pu]", title="SCRX: Amplifier Output")
        lines!(ax5, ref.time, ref[!, Symbol("sCRX.simpleLagLim.y")]; label="OpenIPSL", color=Cycled(5), linewidth=2, alpha=0.7)
        lines!(ax5, ts, sol(ts, idxs=VIndex(:GEN1, :scrx₊amplifier₊out)).u; label="PowerDynamics", color=Cycled(5), linewidth=2, linestyle=:dash)
        axislegend(ax5)

        # Plot 6: SCRX Final EFD Output
        ax6 = Axis(fig[3,2]; xlabel="Time [s]", ylabel="EFD [pu]", title="SCRX: Final Field Voltage")
        lines!(ax6, ref.time, ref[!, Symbol("sCRX.EFD")]; label="OpenIPSL", color=Cycled(6), linewidth=2, alpha=0.7)
        lines!(ax6, ts, sol(ts, idxs=VIndex(:GEN1, :scrx₊EFD_out₊u)).u; label="PowerDynamics", color=Cycled(6), linewidth=2, linestyle=:dash)
        axislegend(ax6)

        fig
    end
    save(joinpath(pkgdir(PowerDynamics),"docs","src","assets","OpenIPSL_valid","SCRX.png"), fig)
end
