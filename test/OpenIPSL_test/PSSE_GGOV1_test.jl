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
    joinpath(pkgdir(PowerDynamics),"test","OpenIPSL_test","GGOV1","modelica_results.csv.gz"),
    DataFrame;
    drop=(i,name) -> contains(string(name), "nrows="),
    silencewarnings=true
)

BUS = let
    # System parameters
    S_b = 100e6
    M_b = 100e6

    Xppd = 0.2
    Xppq = 0.2
    Xpp = 0.2
    Xl = 0.12
    Tpd0 = 5
    Tppd0 = 0.50000E-01
    Tppq0 = 0.1
    H = 4.0000
    D = 0
    Xd = 1.41
    Xq = 1.3500
    Xpd = 0.3
    S10 = 0.1
    S12 = 0.5
    Xpq = 0.6
    Tpq0 = 0.7
    R_a = 0 # from defaults

    # Power flow from OpenIPSL ESST4B test
    v_0 = 1.0
    angle_0 = 0.070620673811798

    # Create GENROU machine with EFD input enabled for exciter control
    @named genrou = PSSE_GENROU(;
        S_b, H, M_b, D,
        Xd, Xq, Xpd, Xpq, Xppd, Xppq, Xl,
        Tpd0, Tpq0, Tppd0, Tppq0,
        S10, S12, R_a,
        pmech_input=true, efd_input=false  # EFD controlled by ESST4B
    )

    @named ggov1 = PSSE_GGOV1(
        R = 0.,
        T_pelec = 1,
        maxerr = 0.05,
        minerr = -0.05,
        Kpgov = 10,
        Kigov = 2,
        Kdgov = 0,
        Tdgov = 1,
        Vmax = 1,
        Vmin = 0.15,
        Tact = 0.5,
        Kturb = 1.5,
        Wfnl = 0.2,
        Tb = 0.1,
        Tc = 0,
        # Teng = 0.5,
        Teng = 0,
        Tfload = 3,
        Kpload = 2,
        Kiload = 0.67,
        Ldref = 1,
        Dm = 0,
        Ropen = 0.1,
        Rclose = -0.1,
        Kimw = 0,
        Aset = 0.1,
        Ka = 10,
        Ta = 0.1,
        # Trate = 0,     # NOT IMPLEMENTED IN OPENIPSL
        db = 0,
        Tsa = 4,
        Tsb = 5,
        # Rup = 99,    # NOT IMPLEMENTED IN OPENIPSL
        # Rdown = -99,
        DELT = 0.0001,
        Flag = 0,
        Rselect = 0
    )

    con = [
        connect(ggov1.PELEC_in, genrou.PELEC_out)
        connect(ggov1.PMECH_out, genrou.PMECH_in)
        connect(genrou.SPEED_out, ggov1.SPEED_in)
    ]

    # Create bus model with proper connections
    busmodel = MTKBus([genrou, ggov1], con; name=:GEN1)
    bm = compile_bus(busmodel, pf=pfSlack(V=v_0, δ=angle_0))

    add_initconstraint!(bm, @initconstraint(:ggov1₊Pref * :ggov1₊R - :ggov1₊Pmwset))
    bm
end

dump_initial_state(BUS)

s0 = OpenIPSL_SMIB(BUS; just_init=true)

sol = OpenIPSL_SMIB(BUS; tol=1, nwtol=1);

## Validation tests for GENROU machine variables (core 5 variables)
@test ref_rms_error(sol, ref, VIndex(:GEN1, :genrou₊w), "gENROU.w") < 2e-5    # Actual: 1.32e-5
@test ref_rms_error(sol, ref, VIndex(:GEN1, :genrou₊delta), "gENROU.delta") < 1e-3   # Actual: 9.34e-4
@test ref_rms_error(sol, ref, VIndex(:GEN1, :genrou₊P), "gENROU.P") < 1e-3    # Actual: 7.16e-4
@test ref_rms_error(sol, ref, VIndex(:GEN1, :genrou₊Q), "gENROU.Q") < 1e-3    # Actual: 9.08e-4
@test ref_rms_error(sol, ref, VIndex(:GEN1, :genrou₊Vt), "gENROU.Vt") < 5e-4  # Actual: 4.23e-4


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
