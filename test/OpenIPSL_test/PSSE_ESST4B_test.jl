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
    joinpath(pkgdir(PowerDynamics),"test","OpenIPSL_test","ESST4B","modelica_results.csv.gz"),
    DataFrame;
    drop=(i,name) -> contains(string(name), "nrows="),
    silencewarnings=true
)

# GENROU+ESST4B bus model parameters from OpenIPSL ESST4B test case
BUS = let
    # System parameters
    S_b = 100e6
    M_b = 100e6

    # GENROU machine parameters from OpenIPSL ESST4B.mo (lines 5-25)
    H = 4.0
    D = 0
    Xd = 1.41
    Xq = 1.35
    Xpd = 0.3
    Xpq = 0.6
    Xppd = 0.2
    Xppq = 0.2
    Xl = 0.12
    Tpd0 = 5.0
    Tpq0 = 0.7
    Tppd0 = 0.05
    Tppq0 = 0.1
    S10 = 0.1
    S12 = 0.5
    R_a = 0

    # Power flow from OpenIPSL ESST4B test
    v_0 = 1.0
    angle_0 = 0.070620673811798

    # Create GENROU machine with EFD input enabled for exciter control
    @named genrou = PSSE_GENROU(;
        S_b, H, M_b, D,
        Xd, Xq, Xpd, Xpq, Xppd, Xppq, Xl,
        Tpd0, Tpq0, Tppd0, Tppq0,
        S10, S12, R_a,
        pmech_input=false, efd_input=true  # EFD controlled by ESST4B
    )

    # Create ESST4B exciter with default parameters
    @named esst4b = PSSE_ESST4B()

    con = [
        # Connect exciter output to machine field voltage
        connect(esst4b.EFD_out, genrou.EFD_in)

        # Connect machine terminal outputs to exciter terminal inputs
        connect(genrou.TERM_VR_out, esst4b.TERM_VR_in)
        connect(genrou.TERM_VI_out, esst4b.TERM_VI_in)
        connect(genrou.TERM_IR_out, esst4b.TERM_IR_in)
        connect(genrou.TERM_II_out, esst4b.TERM_II_in)

        # Connect machine field current to exciter
        connect(genrou.XADIFD_out, esst4b.XADIFD_in)

        # Connect terminal voltage measurement
        connect(genrou.ETERM_out, esst4b.ECOMP_in)
    ]

    # Create bus model with proper connections
    busmodel = MTKBus([genrou, esst4b], con; name=:GEN1)
    bm = compile_bus(busmodel, pf=pfSlack(V=v_0, δ=angle_0))
end

sol = OpenIPSL_SMIB(BUS);

## Validation tests for GENROU machine variables (core 5 variables)
@test ref_rms_error(sol, ref, VIndex(:GEN1, :genrou₊w), "gENROU.w") < 2e-5    # Actual: 1.32e-5
@test ref_rms_error(sol, ref, VIndex(:GEN1, :genrou₊delta), "gENROU.delta") < 1e-3   # Actual: 9.34e-4
@test ref_rms_error(sol, ref, VIndex(:GEN1, :genrou₊P), "gENROU.P") < 1e-3    # Actual: 7.16e-4
@test ref_rms_error(sol, ref, VIndex(:GEN1, :genrou₊Q), "gENROU.Q") < 1e-3    # Actual: 9.08e-4
@test ref_rms_error(sol, ref, VIndex(:GEN1, :genrou₊Vt), "gENROU.Vt") < 5e-4  # Actual: 4.23e-4

## Validation tests for ESST4B exciter variables (signal flow order)
# 1. Transducer output (first in signal chain)
@test ref_rms_error(sol, ref, VIndex(:GEN1, :esst4b₊transducer₊out), "eSST4B.TransducerDelay.y") < 5e-4  # Actual: 3.61e-4
# 2. Voltage PI controller output (proportional + integral + limits)
@test ref_rms_error(sol, ref, VIndex(:GEN1, :esst4b₊vr_out), "eSST4B.limiter.y") < 2e-3             # Actual: 1.48e-3
# 3. Voltage integrator output
@test ref_rms_error(sol, ref, VIndex(:GEN1, :esst4b₊voltage_int₊out), "eSST4B.VR1.y") < 2e-3        # Actual: 1.02e-3
# 4. Thyristor bridge output
@test ref_rms_error(sol, ref, VIndex(:GEN1, :esst4b₊va_out), "eSST4B.VA.y") < 2e-3                  # Actual: 1.48e-3
# 5. Current PI controller output (proportional + integral + limits)
@test ref_rms_error(sol, ref, VIndex(:GEN1, :esst4b₊vm_out), "eSST4B.limiter1.y") < 3e-3            # Actual: 1.93e-3
# 6. Current integrator output
@test ref_rms_error(sol, ref, VIndex(:GEN1, :esst4b₊current_int₊out), "eSST4B.VM1.y") < 1e-3        # Actual: 8.20e-4
# 7. Terminal voltage calculation for rectifier (known VE calculation difference)
@test ref_rms_error(sol, ref, VIndex(:GEN1, :esst4b₊VE), "eSST4B.VE") < 5e-2                        # Actual: 4.27e-2
# 8. Rectifier output before V_BMAX limiting (affected by VE difference)
@test ref_rms_error(sol, ref, VIndex(:GEN1, :esst4b₊vb_signal), "eSST4B.rectifierCommutationVoltageDrop.EFD") < 5e-2  # Actual: 4.27e-2
# 9. Final field voltage output
@test ref_rms_error(sol, ref, VIndex(:GEN1, :esst4b₊EFD_out₊u), "eSST4B.EFD") < 2e-2                # Actual: 1.44e-2

#=
# Create focused comparison plot following signal flow
fig_esst4b = let
    fig = Figure(size=(1400, 1200))
    ts = refine_timeseries(sol.t)

    # Plot 1: Generator Terminal Voltage
    ax1 = Axis(fig[1,1]; xlabel="Time [s]", ylabel="Vt [pu]", title="Generator: Terminal Voltage")
    lines!(ax1, ref.time, ref[!, Symbol("gENROU.Vt")]; label="OpenIPSL", color=Cycled(1), linewidth=2, alpha=0.7)
    lines!(ax1, ts, sol(ts, idxs=VIndex(:GEN1, :genrou₊Vt)).u; label="PowerDynamics", color=Cycled(1), linewidth=2, linestyle=:dash)
    axislegend(ax1)

    # Plot 2: Generator Rotor Angle
    ax2 = Axis(fig[1,2]; xlabel="Time [s]", ylabel="δ [rad]", title="Generator: Rotor Angle")
    lines!(ax2, ref.time, ref[!, Symbol("gENROU.delta")]; label="OpenIPSL", color=Cycled(1), linewidth=2, alpha=0.7)
    lines!(ax2, ts, sol(ts, idxs=VIndex(:GEN1, :genrou₊delta)).u; label="PowerDynamics", color=Cycled(1), linewidth=2, linestyle=:dash)
    axislegend(ax2)

    # Plot 4: ESST4B Voltage PI Controller Output (proportional + integral + limits)
    ax4 = Axis(fig[2,1]; xlabel="Time [s]", ylabel="[pu]", title="ESST4B: Voltage PI Controller")
    lines!(ax4, ref.time, ref[!, Symbol("eSST4B.limiter.y")]; label="OpenIPSL", color=Cycled(3), linewidth=2, alpha=0.7)
    lines!(ax4, ts, sol(ts, idxs=VIndex(:GEN1, :esst4b₊vr_out)).u; label="PowerDynamics", color=Cycled(3), linewidth=2, linestyle=:dash)
    axislegend(ax4)

    # Plot 5: ESST4B Voltage Integrator Component
    ax5 = Axis(fig[2,2]; xlabel="Time [s]", ylabel="[pu]", title="ESST4B: Voltage Integrator")
    lines!(ax5, ref.time, ref[!, Symbol("eSST4B.VR1.y")]; label="OpenIPSL", color=Cycled(4), linewidth=2, alpha=0.7)
    lines!(ax5, ts, sol(ts, idxs=VIndex(:GEN1, :esst4b₊voltage_int₊out)).u; label="PowerDynamics", color=Cycled(4), linewidth=2, linestyle=:dash)
    axislegend(ax5)

    # Plot 6: ESST4B Current PI Controller Output (proportional + integral + limits)
    ax6 = Axis(fig[3,1]; xlabel="Time [s]", ylabel="[pu]", title="ESST4B: Current PI Controller")
    lines!(ax6, ref.time, ref[!, Symbol("eSST4B.limiter1.y")]; label="OpenIPSL", color=Cycled(5), linewidth=2, alpha=0.7)
    lines!(ax6, ts, sol(ts, idxs=VIndex(:GEN1, :esst4b₊vm_out)).u; label="PowerDynamics", color=Cycled(5), linewidth=2, linestyle=:dash)
    axislegend(ax6)

    # Plot 7: ESST4B Current Integrator Component
    ax7 = Axis(fig[3,2]; xlabel="Time [s]", ylabel="[pu]", title="ESST4B: Current Integrator")
    lines!(ax7, ref.time, ref[!, Symbol("eSST4B.VM1.y")]; label="OpenIPSL", color=Cycled(6), linewidth=2, alpha=0.7)
    lines!(ax7, ts, sol(ts, idxs=VIndex(:GEN1, :esst4b₊current_int₊out)).u; label="PowerDynamics", color=Cycled(6), linewidth=2, linestyle=:dash)
    axislegend(ax7)

    # Plot 3: ESST4B Transducer Output (first in signal chain)
    ax3 = Axis(fig[4,1]; xlabel="Time [s]", ylabel="[pu]", title="ESST4B: Transducer Output")
    lines!(ax3, ref.time, ref[!, Symbol("eSST4B.TransducerDelay.y")]; label="OpenIPSL", color=Cycled(2), linewidth=2, alpha=0.7)
    lines!(ax3, ts, sol(ts, idxs=VIndex(:GEN1, :esst4b₊transducer₊out)).u; label="PowerDynamics", color=Cycled(2), linewidth=2, linestyle=:dash)
    axislegend(ax3)

    # Plot 8: ESST4B Terminal Voltage Calculation (VE)
    ax8 = Axis(fig[4,2]; xlabel="Time [s]", ylabel="VE [pu]", title="ESST4B: Terminal Voltage Calculation")
    lines!(ax8, ref.time, ref[!, Symbol("eSST4B.VE")]; label="OpenIPSL", color=Cycled(7), linewidth=2, alpha=0.7)
    lines!(ax8, ts, sol(ts, idxs=VIndex(:GEN1, :esst4b₊VE)).u; label="PowerDynamics", color=Cycled(7), linewidth=2, linestyle=:dash)
    axislegend(ax8)

    # Plot 9: ESST4B Final Field Voltage Output
    ax9 = Axis(fig[5,1]; xlabel="Time [s]", ylabel="EFD [pu]", title="ESST4B: Final Field Voltage")
    lines!(ax9, ref.time, ref[!, Symbol("eSST4B.EFD")]; label="OpenIPSL", color=Cycled(8), linewidth=2, alpha=0.7)
    lines!(ax9, ts, sol(ts, idxs=VIndex(:GEN1, :esst4b₊EFD_out₊u)).u; label="PowerDynamics", color=Cycled(8), linewidth=2, linestyle=:dash)
    axislegend(ax9)

    names(ref)

    fig
end

fig_genrou = let
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
=#
