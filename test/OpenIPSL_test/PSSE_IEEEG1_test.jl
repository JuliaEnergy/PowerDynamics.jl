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

# Load reference data from IEEEG1 simulation
ref = CSV.read(
    joinpath(pkgdir(PowerDynamics),"test","OpenIPSL_test","IEEEG1","modelica_results.csv.gz"),
    DataFrame;
    drop=(i,name) -> contains(string(name), "nrows="),
    silencewarnings=true
)

# GENROU+IEEET1+IEEEG1 bus model parameters from OpenIPSL IEEEG1 test case
BUS = let
    # System parameters
    S_b = 100e6
    M_b = 100e6

    # GENROU machine parameters from OpenIPSL IEEEG1.mo (lines 4-25)
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

    # Power flow from OpenIPSL IEEEG1 test
    v_0 = 1.0
    angle_0 = 0.07068583470577

    # Create GENROU machine with EFD and PMECH inputs enabled for controller operation
    @named genrou = PSSE_GENROU(;
        S_b, H, M_b, D,
        Xd, Xq, Xpd, Xpq, Xppd, Xppq, Xl,
        Tpd0, Tpq0, Tppd0, Tppq0,
        S10, S12, R_a,
        pmech_input=true, efd_input=true  # Both PMECH and EFD controlled by controllers
    )

    # IEEET1 exciter parameters (matching OpenIPSL test exactly, lines 27-39)
    @named ieeet1 = PSSE_IEEET1(;
        T_R=0.06, K_A=200, T_A=0.001, T_E=0.55, K_F=0.06,
        E_1=2.85, S_EE_1=0.3, E_2=3.8, S_EE_2=0.6,
        V_RMAX=2, V_RMIN=-2, K_E=0.1
    )

    # IEEEG1 governor parameters (matching OpenIPSL test exactly, lines 41-52)
    @named ieeeg1 = PSSE_IEEEG1(;
        K=20, T_1=0.15, T_3=0.2, T_4=0.25, K_1=0.25,
        T_5=7.5, K_3=0.25, T_6=0.4, K_5=0.5, T_7=9999
    )

    con = [
        # Exciter connections
        connect(ieeet1.EFD_out, genrou.EFD_in)
        connect(genrou.ETERM_out, ieeet1.ECOMP_in)
        # Governor connections
        connect(genrou.SPEED_out, ieeeg1.SPEED_HP_in)
        connect(ieeeg1.PMECH_HP_out, genrou.PMECH_in)
    ]

    # Create bus model with proper connections
    busmodel = MTKBus([genrou, ieeet1, ieeeg1], con; name=:GEN1)
    bm = compile_bus(busmodel, pf=pfSlack(V=v_0, δ=angle_0))
end

sol = OpenIPSL_SMIB(BUS);

# Core generator variables
@test ref_rms_error(sol, ref, VIndex(:GEN1, :genrou₊w), "gENROU.w") < 2e-6          # Actual: 1.08e-6
@test ref_rms_error(sol, ref, VIndex(:GEN1, :genrou₊delta), "gENROU.delta") < 3e-4  # Actual: 2.58e-4
@test ref_rms_error(sol, ref, VIndex(:GEN1, :genrou₊P), "gENROU.P") < 3e-4          # Actual: 2.30e-4
@test ref_rms_error(sol, ref, VIndex(:GEN1, :genrou₊Vt), "gENROU.Vt") < 1e-5        # Actual: 7.50e-6
@test ref_rms_error(sol, ref, VIndex(:GEN1, :genrou₊Epd), "gENROU.Epd") < 1e-4      # Actual: 8.13e-5
@test ref_rms_error(sol, ref, VIndex(:GEN1, :genrou₊Epq), "gENROU.Epq") < 2e-5      # Actual: 1.51e-5
@test ref_rms_error(sol, ref, VIndex(:GEN1, :genrou₊XadIfd), "gENROU.XadIfd") < 3e-4  # Actual: 2.20e-4

# IEEET1 exciter variables
@test ref_rms_error(sol, ref, VIndex(:GEN1, :ieeet1₊exciter₊EFD), "iEEET1.EFD") < 3e-4              # Actual: 2.15e-4
@test ref_rms_error(sol, ref, VIndex(:GEN1, :ieeet1₊amplifier₊out), "iEEET1.simpleLagLim.y") < 5e-4 # Actual: 3.83e-4
@test ref_rms_error(sol, ref, VIndex(:GEN1, :ieeet1₊derivative_lag₊out), "iEEET1.derivativeLag.y") < 5e-6  # Actual: 3.08e-6
@test ref_rms_error(sol, ref, VIndex(:GEN1, :ieeet1₊transducer₊out), "iEEET1.TransducerDelay.y") < 1e-5   # Actual: 7.22e-6

# IEEEG1 governor variables
@test ref_rms_error(sol, ref, VIndex(:GEN1, :ieeeg1₊PMECH_HP_out₊u), "iEEEG1.PMECH_HP") < 3e-4      # Actual: 2.22e-4
@test ref_rms_error(sol, ref, VIndex(:GEN1, :ieeeg1₊PMECH_LP_out₊u), "iEEEG1.PMECH_LP") < 1e-6      # Actual: 0.0 (perfect match)
@test ref_rms_error(sol, ref, VIndex(:GEN1, :ieeeg1₊hp_turbine₊out), "iEEEG1.imSimpleLag.y") < 3e-4 # Actual: 2.27e-4
@test ref_rms_error(sol, ref, VIndex(:GEN1, :ieeeg1₊ip_turbine₊out), "iEEEG1.imSimpleLag1.y") < 3e-4 # Actual: 2.21e-4
@test ref_rms_error(sol, ref, VIndex(:GEN1, :ieeeg1₊lp1_turbine₊out), "iEEEG1.imSimpleLag2.y") < 3e-4 # Actual: 2.20e-4
@test ref_rms_error(sol, ref, VIndex(:GEN1, :ieeeg1₊lp2_turbine₊out), "iEEEG1.imSimpleLag3.y") < 3e-4 # Actual: 2.14e-4
@test ref_rms_error(sol, ref, VIndex(:GEN1, :ieeeg1₊valve_integrator₊out), "iEEEG1.limIntegrator.y") < 3e-4 # Actual: 2.29e-4
@test ref_rms_error(sol, ref, VIndex(:GEN1, :ieeeg1₊leadlag₊out), "iEEEG1.imLeadLag.y") < 2e-5      # Actual: 1.53e-5

#=
fig_ieeeg1 = let
    fig = Figure(size=(1800, 1600))  # Larger figure for more plots
    # ts = refine_timeseries(sol.t)[1:100]
    ts = refine_timeseries(sol.t)

    # Plot 1: Generator - Angular frequency ω
    ax1 = Axis(fig[1,1]; xlabel="Time [s]", ylabel="ω [pu]", title="Generator: Angular Frequency")
    lines!(ax1, ref.time, ref[!, Symbol("gENROU.w")]; label="gENROU.w", color=Cycled(1), linewidth=2, alpha=0.7)
    lines!(ax1, ts, sol(ts, idxs=VIndex(:GEN1, :genrou₊w)).u; label="genrou.w", color=Cycled(1), linewidth=2, linestyle=:dash)
    axislegend(ax1)

    # Plot 2: Generator - Rotor angle δ
    ax2 = Axis(fig[1,2]; xlabel="Time [s]", ylabel="δ [rad]", title="Generator: Rotor Angle")
    lines!(ax2, ref.time, ref[!, Symbol("gENROU.delta")]; label="gENROU.delta", color=Cycled(1), linewidth=2, alpha=0.7)
    lines!(ax2, ts, sol(ts, idxs=VIndex(:GEN1, :genrou₊delta)).u; label="genrou.delta", color=Cycled(1), linewidth=2, linestyle=:dash)
    axislegend(ax2)

    # Plot 3: Governor - HP Power Output
    ax3 = Axis(fig[2,1]; xlabel="Time [s]", ylabel="PMECH_HP [pu]", title="Governor: HP Power Output")
    lines!(ax3, ref.time, ref[!, Symbol("iEEEG1.PMECH_HP")]; label="iEEEG1.PMECH_HP", color=Cycled(2), linewidth=2, alpha=0.7)
    lines!(ax3, ts, sol(ts, idxs=VIndex(:GEN1, :ieeeg1₊PMECH_HP_out₊u)).u; label="ieeeg1.PMECH_HP_out", color=Cycled(2), linewidth=2, linestyle=:dash)
    axislegend(ax3)

    # Plot 4: Governor - LP Power Output
    ax4 = Axis(fig[2,2]; xlabel="Time [s]", ylabel="PMECH_LP [pu]", title="Governor: LP Power Output")
    lines!(ax4, ref.time, ref[!, Symbol("iEEEG1.PMECH_LP")]; label="iEEEG1.PMECH_LP", color=Cycled(3), linewidth=2, alpha=0.7)
    lines!(ax4, ts, sol(ts, idxs=VIndex(:GEN1, :ieeeg1₊PMECH_LP_out₊u)).u; label="ieeeg1.PMECH_LP_out", color=Cycled(3), linewidth=2, linestyle=:dash)
    axislegend(ax4)

    # Plot 5: Governor - Turbine Stages
    ax5 = Axis(fig[3,1]; xlabel="Time [s]", ylabel="[pu]", title="Governor: Turbine Stages")
    lines!(ax5, ref.time, ref[!, Symbol("iEEEG1.imSimpleLag.y")]; label="iEEEG1.imSimpleLag.y (HP)", color=Cycled(4), linewidth=2, alpha=0.7)
    lines!(ax5, ts, sol(ts, idxs=VIndex(:GEN1, :ieeeg1₊hp_turbine₊out)).u; label="ieeeg1.hp_turbine.out (HP)", color=Cycled(4), linewidth=2, linestyle=:dash)
    lines!(ax5, ref.time, ref[!, Symbol("iEEEG1.imSimpleLag1.y")]; label="iEEEG1.imSimpleLag1.y (IP)", color=Cycled(5), linewidth=2, alpha=0.7)
    lines!(ax5, ts, sol(ts, idxs=VIndex(:GEN1, :ieeeg1₊ip_turbine₊out)).u; label="ieeeg1.ip_turbine.out (IP)", color=Cycled(5), linewidth=2, linestyle=:dash)
    axislegend(ax5)

    # Plot 6: Governor - Valve Position
    ax6 = Axis(fig[3,2]; xlabel="Time [s]", ylabel="Valve Position [pu]", title="Governor: Valve Position")
    lines!(ax6, ref.time, ref[!, Symbol("iEEEG1.limIntegrator.y")]; label="iEEEG1.limIntegrator.y", color=Cycled(6), linewidth=2, alpha=0.7)
    lines!(ax6, ts, sol(ts, idxs=VIndex(:GEN1, :ieeeg1₊valve_integrator₊out)).u; label="ieeeg1.valve_integrator.out", color=Cycled(6), linewidth=2, linestyle=:dash)
    axislegend(ax6)

    # Plot 7: Exciter - Field Voltage Output
    ax7 = Axis(fig[4,1]; xlabel="Time [s]", ylabel="EFD [pu]", title="Exciter: Field Voltage Output")
    lines!(ax7, ref.time, ref[!, Symbol("iEEET1.EFD")]; label="iEEET1.EFD", color=Cycled(7), linewidth=2, alpha=0.7)
    lines!(ax7, ts, sol(ts, idxs=VIndex(:GEN1, :ieeet1₊exciter₊EFD)).u; label="ieeet1.exciter.EFD", color=Cycled(7), linewidth=2, linestyle=:dash)
    axislegend(ax7)

    # Plot 8: Governor - Lead-Lag Output
    ax8 = Axis(fig[4,2]; xlabel="Time [s]", ylabel="Lead-Lag Output [pu]", title="Governor: Lead-Lag Compensator")
    lines!(ax8, ref.time, ref[!, Symbol("iEEEG1.imLeadLag.y")]; label="iEEEG1.imLeadLag.y", color=Cycled(8), linewidth=2, alpha=0.7)
    lines!(ax8, ts, sol(ts, idxs=VIndex(:GEN1, :ieeeg1₊leadlag₊out)).u; label="ieeeg1.leadlag.out", color=Cycled(8), linewidth=2, linestyle=:dash)
    axislegend(ax8)

    # Plot 9: Exciter - Voltage Transducer
    ax9 = Axis(fig[5,1]; xlabel="Time [s]", ylabel="[pu]", title="Exciter: Voltage Transducer")
    lines!(ax9, ref.time, ref[!, Symbol("iEEET1.TransducerDelay.y")]; label="iEEET1.TransducerDelay.y", color=Cycled(9), linewidth=2, alpha=0.7)
    lines!(ax9, ts, sol(ts, idxs=VIndex(:GEN1, :ieeet1₊transducer₊out)).u; label="ieeet1.transducer.out", color=Cycled(9), linewidth=2, linestyle=:dash)
    axislegend(ax9)

    # Plot 10: Exciter - Amplifier Output
    ax10 = Axis(fig[5,2]; xlabel="Time [s]", ylabel="[pu]", title="Exciter: Amplifier Output")
    lines!(ax10, ref.time, ref[!, Symbol("iEEET1.simpleLagLim.y")]; label="iEEET1.simpleLagLim.y", color=Cycled(10), linewidth=2, alpha=0.7)
    lines!(ax10, ts, sol(ts, idxs=VIndex(:GEN1, :ieeet1₊amplifier₊out)).u; label="ieeet1.amplifier.out", color=Cycled(10), linewidth=2, linestyle=:dash)
    axislegend(ax10)

    # Plot 11: Exciter - Derivative Feedback
    ax11 = Axis(fig[6,1]; xlabel="Time [s]", ylabel="[pu]", title="Exciter: Derivative Feedback")
    lines!(ax11, ref.time, ref[!, Symbol("iEEET1.derivativeLag.y")]; label="iEEET1.derivativeLag.y", color=Cycled(11), linewidth=2, alpha=0.7)
    lines!(ax11, ts, sol(ts, idxs=VIndex(:GEN1, :ieeet1₊derivative_lag₊out)).u; label="ieeet1.derivative_lag.out", color=Cycled(11), linewidth=2, linestyle=:dash)
    axislegend(ax11)

    # Plot 12: Generator - Terminal Voltage
    ax12 = Axis(fig[6,2]; xlabel="Time [s]", ylabel="Vt [pu]", title="Generator: Terminal Voltage")
    lines!(ax12, ref.time, ref[!, Symbol("gENROU.Vt")]; label="gENROU.Vt", color=Cycled(12), linewidth=2, alpha=0.7)
    lines!(ax12, ts, sol(ts, idxs=VIndex(:GEN1, :genrou₊Vt)).u; label="genrou.Vt", color=Cycled(12), linewidth=2, linestyle=:dash)
    axislegend(ax12)

    fig
end
=#
