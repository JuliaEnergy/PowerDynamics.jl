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
using Test

#=
This model has several problems:
- does not support Rup/Rdown (same as OpenIPSL)
- does not support delay (unless OpenIPSL)
- there seems to be bugs in the OpenIPSL Turbine implementation:
  - SPEED+1 is doen in an external block but once again in the dm select and speed flag
  - the dm select bloc is connected wrongly, it allways passes to
  - since DM=0 and FLAG=0 this does not show up in the recorded data

The initialization is quite tricky to. In this file we do a completly manual initialization
of the guess values.Crucially, the parameter
- Pref is set to 0, i don't quite know what its supposed to do
  but there is a ninitialisation ambiguity between Pref and the integrator behind PMWSet-Pel
- the parameter LDRef was left free for initialization, otherwise we cannot guarantee
  that the KILOAD integrator is in steady state
- the state of the load integrator KILoad is "free" and musst be chosen manually
=#

@warn "Carfull, the GGOV1 reference data has known bugs. Also, the OpenIPSL was automaticially \
       modified while generating the reference data to REMOVE THE DELAY BLOCK, since we cannot \
       model this in PowerDynamics (yet). The delay played a MAJOR role in the original data!"
ref = CSV.read(
    joinpath(pkgdir(PowerDynamics),"test","OpenIPSL_test","GGOV1","modelica_results_modified.csv.gz"),
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

    @named ggov1 = PSSE_GGOV1_EXPERIMENTAL(
        Pref=0,
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
        # Ldref = 1,
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

    # load/temp limiter initialization
    # XXX This is a choice! It won't affect the result so we can chose whatever
    # we like here, however it defines how far away from the "limiting" we start
    # with our simulation
    set_default!(bm, :ggov1₊load_limiter₊integral_state, 1)

    guessformulas = @guessformula begin
        # base equations
        Pe0 = - 1 * (:busbar₊u_r*:busbar₊i_r + :busbar₊u_i*:busbar₊i_i)
        Pmech0 = Pe0
        :ggov1₊Pmwset = Pe0

        turbine_output = if :ggov1₊Dm ≥ 0
            Pmech0 + 1 * :ggov1₊Dm
        else
            Pmech0
        end
        :ggov1₊turbine₊turbine_dynamics₊internal = turbine_output
        fuel_flow = turbine_output/ :ggov1₊Kturb + :ggov1₊Wfnl
        :ggov1₊turbine₊valve_integrator = fuel_flow

        TEXM = fuel_flow # ifels of DM is 1 either way
        :ggov1₊turbine₊temp_leadlag₊internal = TEXM
        :ggov1₊turbine₊temp_lag₊out = TEXM

        # pid govenor initialization
        rsel_contribution = if :ggov1₊_Rselect_static == 0
            0
        elseif :ggov1₊_Rselect_static
            :ggov1₊R * Pe0
        else
            :ggov1₊R * fuel_flow
        end
        Pref_plus_int = -rsel_contribution
        :ggov1₊pid_governor₊power_controller₊out = Pref_plus_int - :ggov1₊Pref
        :ggov1₊pid_governor₊power_transducer₊out = Pe0 # power measurment
        :ggov1₊pid_governor₊pid_integral_state = fuel_flow # integral state
        # :ggov1₊pid_governor₊speed_derivative₊internal = 0

        # accel limiter initialization
        # the derivateive is zero which is fine
        # however, the accelerator intoruced an algebraic loop, we need guess its output
        :ggov1₊accel_limiter₊FSR = fuel_flow # same as govenor output

        :ggov1₊Ldref = (TEXM - :ggov1₊Wfnl) * :ggov1₊Kturb
    end
    add_guessformula!(bm, guessformulas)

    bm
end
sol = OpenIPSL_SMIB(BUS);

## Validation tests for GENROU machine variables
@test ref_rms_error(sol, ref, VIndex(:GEN1, :genrou₊w), "gENROU.w") < 1e-6
@test ref_rms_error(sol, ref, VIndex(:GEN1, :genrou₊delta), "gENROU.delta") < 8e-5
@test ref_rms_error(sol, ref, VIndex(:GEN1, :genrou₊P), "gENROU.P") < 5e-5
@test ref_rms_error(sol, ref, VIndex(:GEN1, :genrou₊Q), "gENROU.Q") < 2e-5
@test ref_rms_error(sol, ref, VIndex(:GEN1, :genrou₊Vt), "gENROU.Vt") < 7e-6

## Validation tests for GGOV1 governor
@test ref_rms_error(sol, ref, VIndex(:GEN1, :ggov1₊PMECH_out₊u), "gGOV1.PMECH") < 5e-5
@test ref_rms_error(sol, ref, VIndex(:GEN1, :ggov1₊load_limiter₊FSRT), "gGOV1.gGOV1_Temp.FSRT") < 1e-5
@test ref_rms_error(sol, ref, VIndex(:GEN1, :ggov1₊accel_limiter₊FSRA), "gGOV1.gGOV1_Accel.FSRA") < 3e-5
@test ref_rms_error(sol, ref, VIndex(:GEN1, :ggov1₊pid_governor₊FSRN), "gGOV1.gGOV1_Power.FSRN") < 3e-5
@test ref_rms_error(sol, ref, VIndex(:GEN1, :ggov1₊turbine₊VSTROKE), "gGOV1.gGOV1_Turb.VSTROKE") < 3e-5

# This is broken becaus we left LDREF free thus achieving a true steady state
# in OpenIPSL, they fix LDREF=1 so the KILOAD integrator is constantly integrating up
@test_broken ref_rms_error(sol, ref, VIndex(:GEN1, :ggov1₊turbine₊TEXM), "gGOV1.gGOV1_Turb.TEXM") < 1e-3

if isdefined(Main, :EXPORT_FIGURES) && Main.EXPORT_FIGURES
    fig = let
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

        # Plot 3: Angular frequency ω
        ax3 = Axis(fig[2,1]; xlabel="Time [s]", ylabel="ω [pu]", title="Generator: Angular Frequency")
        lines!(ax3, ref.time, ref[!, Symbol("gENROU.w")]; label="OpenIPSL", color=Cycled(2), linewidth=2, alpha=0.7)
        lines!(ax3, ts, sol(ts, idxs=VIndex(:GEN1, :genrou₊w)).u; label="PowerDynamics", color=Cycled(2), linewidth=2, linestyle=:dash)
        axislegend(ax3)

        # Plot 4: Active power P
        ax4 = Axis(fig[2,2]; xlabel="Time [s]", ylabel="P [pu]", title="Generator: Active Power")
        lines!(ax4, ref.time, ref[!, Symbol("gENROU.P")]; label="OpenIPSL", color=Cycled(2), linewidth=2, alpha=0.7)
        lines!(ax4, ts, sol(ts, idxs=VIndex(:GEN1, :genrou₊P)).u; label="PowerDynamics", color=Cycled(2), linewidth=2, linestyle=:dash)
        axislegend(ax4)

        # Plot 5: Reactive power Q
        ax5 = Axis(fig[3,1]; xlabel="Time [s]", ylabel="Q [pu]", title="Generator: Reactive Power")
        lines!(ax5, ref.time, ref[!, Symbol("gENROU.Q")]; label="OpenIPSL", color=Cycled(3), linewidth=2, alpha=0.7)
        lines!(ax5, ts, sol(ts, idxs=VIndex(:GEN1, :genrou₊Q)).u; label="PowerDynamics", color=Cycled(3), linewidth=2, linestyle=:dash)
        axislegend(ax5)

        # Plot 6: Turbine mechanical power (GGOV1 output)
        ax6 = Axis(fig[3,2]; xlabel="Time [s]", ylabel="PMECH [pu]", title="GGOV1: Mechanical Power Output")
        lines!(ax6, ref.time, ref[!, Symbol("gGOV1.PMECH")]; label="OpenIPSL", color=Cycled(4), linewidth=2, alpha=0.7)
        lines!(ax6, ts, sol(ts, idxs=VIndex(:GEN1, :ggov1₊PMECH_out₊u)).u; label="PowerDynamics", color=Cycled(4), linewidth=2, linestyle=:dash)
        axislegend(ax6)

        # Plot 7: Exhaust temperature (GGOV1 TEXM)
        ax7 = Axis(fig[4,1]; xlabel="Time [s]", ylabel="TEXM [pu]", title="GGOV1: Exhaust Temperature")
        lines!(ax7, ref.time, ref[!, Symbol("gGOV1.gGOV1_Turb.TEXM")]; label="OpenIPSL", color=Cycled(5), linewidth=2, alpha=0.7)
        lines!(ax7, ts, sol(ts, idxs=VIndex(:GEN1, :ggov1₊turbine₊TEXM)).u; label="PowerDynamics", color=Cycled(5), linewidth=2, linestyle=:dash)
        axislegend(ax7)

        fig
    end
    save(joinpath(pkgdir(PowerDynamics),"docs","src","assets","OpenIPSL_valid","GGOV1.png"), fig)
end
