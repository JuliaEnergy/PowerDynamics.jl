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

# Load reference data from HYGOV simulation
ref = CSV.read(
    joinpath(pkgdir(PowerDynamics),"test","OpenIPSL_test","HYGOV","modelica_results.csv.gz"),
    DataFrame;
    drop=(i,name) -> contains(string(name), "nrows="),
    silencewarnings=true
)

#=
Once again we see a problem with the govenor initialization. The core culprit are the two
itegral states for the
- limited valve position
- water flow
both of which can be initialized at various points (when n_ref and h0 are free)

So essentially we need to pin on of the (nref, valve position) and one of the (h0, water int)
states to get a unique initialization! In this example, i've chose the same as in OpenIPSL

We take h0=1 from open ipsl and fix the vale position from the csv
=#

# HYGOV test system parameters from OpenIPSL HYGOV test case
BUS = let
    # Power flow from OpenIPSL HYGOV test (from HYGOV.mo lines 20-22)
    v_0 = 1.0
    angle_0 = 0.07068583470577

    # Create GENSAL machine with EFD and PMECH inputs enabled
    @named gensal = PSSE_GENSAL(;
        # System parameters
        S_b = 100e6,
        M_b = 100e6,
        # GENSAL machine parameters from OpenIPSL HYGOV.mo (lines 4-19)
        H = 4.41,
        D = 0,
        Xd = 1.22,
        Xq = 0.76,
        Xpd = 0.297,
        Xppd = 0.2,
        Xppq = 0.2,
        Xl = 0.12,
        Tpd0 = 6.7,
        Tppd0 = 0.028,
        Tppq0 = 0.0358,
        S10 = 0.186,
        S12 = 0.802,
        R_a = 0,
        pmech_input = true,  # PMECH controlled by HYGOV
        efd_input = true     # EFD controlled by SCRX
    )

    # Create SCRX exciter with parameters matching OpenIPSL test (lines 39-47)
    @named scrx = PSSE_SCRX(;
        T_B = 10,
        K = 100,
        T_E = 0.05,
        E_MIN = 0,
        E_MAX = 5,
        r_cr_fd = 0,
        C_SWITCH = false,  # Bus fed mode
        vothsg_input = false,  # No other signal input
        vuel_input = false,    # No under-excitation limiter
        voel_input = false     # No over-excitation limiter
    )

    # Create HYGOV governor with parameters matching OpenIPSL test (lines 23-35)
    @named hygov = PSSE_HYGOV(;
        VELM = 0.02,
        G_MAX = 0.415,
        R = 0.05,
        r = 0.3,
        T_r = 5,
        T_f = 0.05,
        T_g = 0.5,
        G_MIN = 0,
        T_w = 1.25,
        A_t = 1.2,
        D_turb = 0.2,
        q_NL = 0.08
    )

    con = [
        # Exciter connections (SCRX to GENSAL)
        connect(scrx.EFD_out, gensal.EFD_in)
        connect(gensal.XADIFD_out, scrx.XADIFD_in)
        connect(gensal.ETERM_out, scrx.ECOMP_in)

        # Governor connections (HYGOV to GENSAL)
        connect(hygov.PMECH_out, gensal.PMECH_in)
        connect(gensal.SPEED_out, hygov.SPEED_in)
    ]

    # Create bus model with proper connections
    busmodel = MTKBus([gensal, scrx, hygov], con; name=:GEN1)
    bm = compile_bus(busmodel, pf=pfSlack(V=v_0, δ=angle_0))
    # set_default!(bm, :hygov₊n_ref, 0.0206)
    set_default!(bm, :hygov₊position_limiter₊out, ref[1, "hYGOV.Position_Limiter.y"])
    # set_default!(bm, :hygov₊water_integrator₊out, ref[1, "hYGOV.q.y"])
    # set_default!(bm, :hygov₊filter₊out, 0)
    # set_default!(bm, :hygov₊water_integrator₊out, 0.41)
    # c = @initconstraint(:hygov₊PMECH_out₊u/(:hygov₊A_t*:hygov₊h0) + :hygov₊q_NL - :hygov₊water_flow)
    # add_initconstraint!(bm, c)
    bm
end

sol = OpenIPSL_SMIB(BUS, tol=1e-4, nwtol=1e-4);

## Validation tests for GENSAL machine variables (core 3 variables)
@test ref_rms_error(sol, ref, VIndex(:GEN1, :gensal₊w), "gENSAL.w") < 1e-5
@test ref_rms_error(sol, ref, VIndex(:GEN1, :gensal₊delta), "gENSAL.delta") < 5e-4
@test ref_rms_error(sol, ref, VIndex(:GEN1, :gensal₊Vt), "gENSAL.Vt") < 3e-5

## Validation tests for SCRX exciter variables
@test ref_rms_error(sol, ref, VIndex(:GEN1, :scrx₊EFD_out₊u), "sCRX.EFD") < 2e-4

## Validation tests for HYGOV governor variables (main states)
@test ref_rms_error(sol, ref, VIndex(:GEN1, :hygov₊filter₊out), "hYGOV.SimpleLag1.y") < 5e-5
@test ref_rms_error(sol, ref, VIndex(:GEN1, :hygov₊position_limiter₊out), "hYGOV.Position_Limiter.y") < 5e-4
@test ref_rms_error(sol, ref, VIndex(:GEN1, :hygov₊servo₊out), "hYGOV.g.y") < 3e-4
@test ref_rms_error(sol, ref, VIndex(:GEN1, :hygov₊water_integrator₊out), "hYGOV.q.y") < 3e-4
@test ref_rms_error(sol, ref, VIndex(:GEN1, :hygov₊turbine_head), "hYGOV.H") < 3e-4
@test ref_rms_error(sol, ref, VIndex(:GEN1, :hygov₊PMECH_out₊u), "hYGOV.PMECH") < 3e-4

#=
fig_hygov = let
    fig = Figure(size=(1600, 1200))
    ts = refine_timeseries(sol.t)

    # Plot 1: Generator Rotor Angle
    ax1 = Axis(fig[1,1]; xlabel="Time [s]", ylabel="δ [rad]", title="GENSAL: Rotor Angle")
    lines!(ax1, ref.time, ref[!, Symbol("gENSAL.delta")]; label="OpenIPSL", color=Cycled(1), linewidth=2, alpha=0.7)
    lines!(ax1, ts, sol(ts, idxs=VIndex(:GEN1, :gensal₊delta)).u; label="PowerDynamics", color=Cycled(1), linewidth=2, linestyle=:dash)
    axislegend(ax1)

    # Plot 2: Generator Speed Deviation
    ax2 = Axis(fig[1,2]; xlabel="Time [s]", ylabel="ω [pu]", title="GENSAL: Speed Deviation")
    lines!(ax2, ref.time, ref[!, Symbol("gENSAL.w")]; label="OpenIPSL", color=Cycled(2), linewidth=2, alpha=0.7)
    lines!(ax2, ts, sol(ts, idxs=VIndex(:GEN1, :gensal₊w)).u; label="PowerDynamics", color=Cycled(2), linewidth=2, linestyle=:dash)
    axislegend(ax2)

    # Plot 3: Generator Terminal Voltage
    ax3 = Axis(fig[1,3]; xlabel="Time [s]", ylabel="Vt [pu]", title="GENSAL: Terminal Voltage")
    lines!(ax3, ref.time, ref[!, Symbol("gENSAL.Vt")]; label="OpenIPSL", color=Cycled(3), linewidth=2, alpha=0.7)
    lines!(ax3, ts, sol(ts, idxs=VIndex(:GEN1, :gensal₊Vt)).u; label="PowerDynamics", color=Cycled(3), linewidth=2, linestyle=:dash)
    axislegend(ax3)

    # Plot 4: SCRX Field Voltage Output
    ax4 = Axis(fig[2,1]; xlabel="Time [s]", ylabel="EFD [pu]", title="SCRX: Field Voltage")
    lines!(ax4, ref.time, ref[!, Symbol("sCRX.EFD")]; label="OpenIPSL", color=Cycled(4), linewidth=2, alpha=0.7)
    lines!(ax4, ts, sol(ts, idxs=VIndex(:GEN1, :scrx₊EFD_out₊u)).u; label="PowerDynamics", color=Cycled(4), linewidth=2, linestyle=:dash)
    axislegend(ax4)

    # Plot 5: HYGOV Filter Output
    ax5 = Axis(fig[2,2]; xlabel="Time [s]", ylabel="[pu]", title="HYGOV: Filter Output")
    lines!(ax5, ref.time, ref[!, Symbol("hYGOV.SimpleLag1.y")]; label="OpenIPSL", color=Cycled(5), linewidth=2, alpha=0.7)
    lines!(ax5, ts, sol(ts, idxs=VIndex(:GEN1, :hygov₊filter₊out)).u; label="PowerDynamics", color=Cycled(5), linewidth=2, linestyle=:dash)
    axislegend(ax5)

    # Plot 6: HYGOV Desired Gate Position
    ax6 = Axis(fig[2,3]; xlabel="Time [s]", ylabel="[pu]", title="HYGOV: Desired Gate Position")
    lines!(ax6, ref.time, ref[!, Symbol("hYGOV.Position_Limiter.y")]; label="OpenIPSL", color=Cycled(6), linewidth=2, alpha=0.7)
    lines!(ax6, ts, sol(ts, idxs=VIndex(:GEN1, :hygov₊position_limiter₊out)).u; label="PowerDynamics", color=Cycled(6), linewidth=2, linestyle=:dash)
    axislegend(ax6)

    # Plot 7: HYGOV Actual Gate Position
    ax7 = Axis(fig[3,1]; xlabel="Time [s]", ylabel="G [pu]", title="HYGOV: Actual Gate Position")
    lines!(ax7, ref.time, ref[!, Symbol("hYGOV.g.y")]; label="OpenIPSL", color=Cycled(7), linewidth=2, alpha=0.7)
    lines!(ax7, ts, sol(ts, idxs=VIndex(:GEN1, :hygov₊servo₊out)).u; label="PowerDynamics", color=Cycled(7), linewidth=2, linestyle=:dash)
    axislegend(ax7)

    # Plot 8: HYGOV Water Flow
    ax8 = Axis(fig[3,2]; xlabel="Time [s]", ylabel="Q [pu]", title="HYGOV: Water Flow")
    lines!(ax8, ref.time, ref[!, Symbol("hYGOV.q.y")]; label="OpenIPSL", color=Cycled(8), linewidth=2, alpha=0.7)
    lines!(ax8, ts, sol(ts, idxs=VIndex(:GEN1, :hygov₊water_integrator₊out)).u; label="PowerDynamics", color=Cycled(8), linewidth=2, linestyle=:dash)
    axislegend(ax8)

    # Plot 9: HYGOV Turbine Head
    ax9 = Axis(fig[3,3]; xlabel="Time [s]", ylabel="H [pu]", title="HYGOV: Turbine Head")
    lines!(ax9, ref.time, ref[!, Symbol("hYGOV.H")]; label="OpenIPSL", color=Cycled(9), linewidth=2, alpha=0.7)
    lines!(ax9, ts, sol(ts, idxs=VIndex(:GEN1, :hygov₊turbine_head)).u; label="PowerDynamics", color=Cycled(9), linewidth=2, linestyle=:dash)
    axislegend(ax9)

    # Plot 10: HYGOV Mechanical Power Output
    ax10 = Axis(fig[4,1]; xlabel="Time [s]", ylabel="PMECH [pu]", title="HYGOV: Mechanical Power")
    lines!(ax10, ref.time, ref[!, Symbol("hYGOV.PMECH")]; label="OpenIPSL", color=Cycled(10), linewidth=2, alpha=0.7)
    lines!(ax10, ts, sol(ts, idxs=VIndex(:GEN1, :hygov₊PMECH_out₊u)).u; label="PowerDynamics", color=Cycled(10), linewidth=2, linestyle=:dash)
    axislegend(ax10)

    fig
end
=#
