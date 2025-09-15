using PowerDynamics
using PowerDynamics.Library
using PowerDynamicsTesting
using ModelingToolkit
using NetworkDynamics
using Graphs
using OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqNonlinearSolve
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkit: ModelingToolkit as MTK
using OrderedCollections
using DiffEqCallbacks
using CairoMakie
using Test

@info "Start Library tests"

@testset "PiLine" begin
    @named branchA = PiLine(X=0.022522, R=0.01)
    @named branchB = PiLine(X=0.04189, R=0.02)
    line = compile_line(MTKLine(branchA, branchB));
    toi = line_between_slacks(line)
    # isinteractive() && plottoi(toi)
    @reftest "PiLine_1" toi

    @named branchA = PiLine(X=0.022522, R=0.01)
    line = compile_line(MTKLine(branchA));
    toi = line_between_slacks(line)
    # isinteractive() && plottoi(toi)
    @reftest "PiLine_2" toi
end

@testset "Swing bus" begin
    @named swing = Swing(Pm=1, D=0.1, M=0.005, θ=0, ω=1, V=1)
    bus = compile_bus(MTKBus(swing));
    toi = bus_on_slack(bus)
    # isinteractive() && plottoi(toi)
    @reftest "SwingBus_1" toi

    # swing bus with load
    @named swing = Swing(Pm=1, D=0.1, M=0.005, θ=0, ω=1, V=1)
    @named pqload = PQLoad(Pset=-0.5, Qset=-0.2)
    bm = MTKBus(swing, pqload)
    @test length(full_equations(simplify_mtkbus(bm))) == 2
    bus = compile_bus(bm)
    toi = bus_on_slack(bus)
    toi["active power"]["electric power of swing"] = VIndex(2,:swing₊Pel)
    toi["active power"]["electric power of load"] = VIndex(2,:pqload₊P)
    toi["reactive power"]["electric power of load"] = VIndex(2,:pqload₊Q)
    # isinteractive() && plottoi(toi)
    @reftest "Swing_and_load" toi
end

@testset "SauerPai Generator" begin
    @mtkmodel GenBus begin
        @components begin
            machine = PowerDynamics.Library.SauerPaiMachine(;
                vf_input=false,
                τ_m_input=false,
                S_b=100,
                V_b=18,
                ω_b=2π*60,
                X_d=0.146, X′_d=0.0608, X″_d=0.06,
                X_q=0.1, X′_q=0.0969, X″_q=0.06,
                R_s=0.000124,
                X_ls=0.01460,
                T′_d0=8.96,
                T″_d0=0.01,
                T′_q0=0.31,
                T″_q0=0.01,
                H=23.64,
            )
            busbar = BusBar()
        end
        @equations begin
            connect(machine.terminal, busbar.terminal)
        end
    end
    @named mtkbus = GenBus()
    bus = compile_bus(mtkbus)

    set_voltage!(bus; mag=1.017, arg=0.0295)
    set_current!(bus; P=0.716, Q=0.3025)
    initialize_component!(bus)

    toi = bus_on_slack(bus; tmax=600, toilength=10_000)
    # isinteractive() && plottoi(toi)

    @reftest "SauerPai" toi
end

@testset "SauerPai Generator with AVR and GOV" begin
    A, B = Library.solve_ceilf(3.3 => 0.6602, 4.5 => 4.2662)

    @mtkmodel GenBus begin
        @components begin
            machine = PowerDynamics.Library.SauerPaiMachine(;
                vf_input=true,
                τ_m_input=true,
                S_b=100,
                V_b=18,
                ω_b=2π*60,
                X_d=0.146, X′_d=0.0608, X″_d=0.06,
                X_q=0.1, X′_q=0.0969, X″_q=0.06,
                R_s=0.000124,
                X_ls=0.01460,
                T′_d0=8.96,
                T″_d0=0.01,
                T′_q0=0.31,
                T″_q0=0.01,
                H=23.64,
            )
            avr = AVRTypeI(
                vr_min=-5,
                vr_max=5,
                Ka=20,
                Ta=0.2,
                Kf=0.063,
                Tf=0.35,
                Ke=1,
                Te=0.314,
                A, B,
                tmeas_lag=false)
            gov = TGOV1(
                R=0.05,
                T1=0.05,
                T2=2.1,
                T3=7.0,
                DT=0,
                V_max=5,
                V_min=-5)
            busbar = BusBar()
        end
        @equations begin
            connect(machine.terminal, busbar.terminal)
            connect(machine.v_mag_out, avr.v_mag)
            connect(avr.vf, machine.vf_in)
            connect(gov.τ_m, machine.τ_m_in)
            connect(machine.ωout, gov.ω_meas)
        end
    end
    @named mtkbus = GenBus()

    bus = compile_bus(mtkbus)
    set_voltage!(bus; mag=1.017, arg=0.0295)
    set_current!(bus; P=0.716, Q=0.3025)
    initialize_component!(bus)

    toi = bus_on_slack(bus; tmax=600, toilength=10_000)
    # isinteractive() && plottoi(toi)

    @reftest "SauerPai_with_AVR_GOV" toi
end

@testset "test loads" begin
    @named load = PQLoad(Pset=-0.5, Qset=-0.5)
    bus = compile_bus(MTKBus(load));
    toi = bus_on_slack(bus)
    # isinteractive() && plottoi(toi)
    @reftest "PQLoad_1" toi

    @named load = VoltageDependentLoad(Pset=-0.5, Qset=-0.5, Vn=1, αP=1, αQ=1)
    bus = compile_bus(MTKBus(load));
    toi = bus_on_slack(bus)
    # isinteractive() && plottoi(toi)
    @reftest "VoltageDependentLoad_1" toi

    Sload = -0.5-im*0.5
    Vset = 1.0
    Y = -conj(Sload)/Vset^2
    @named load = ConstantYLoad(B=imag(Y), G=real(Y))
    bus = compile_bus(MTKBus(load));
    toi = bus_on_slack(bus)
    # isinteractive() && plottoi(toi)
    @reftest "ConstantYLoad" toi
end

@testset "test zip load" begin
    # constant power
    @named load = ZIPLoad(Pset=-1, Qset=-1, Vset=1,
                          KpZ=0, KpI=0, KpC=1,
                          KqZ=0, KqI=0, KqC=1)
    bus = compile_bus(MTKBus(load))
    toi = bus_on_slack(bus)
    # isinteractive() && plottoi(toi)
    @reftest "ZIPLoad_constant_power" toi

    # constant current
    @named load = ZIPLoad(Pset=-1, Qset=-1, Vset=1,
                          KpZ=0, KpI=1, KpC=0,
                          KqZ=0, KqI=1, KqC=0)
    bus = compile_bus(MTKBus(load))
    toi = bus_on_slack(bus)
    toi["current"] = OrderedDict(
        "current magnitude at bus" => VIndex(2,:busbar₊i_mag),
        "real current" => VIndex(2,:load₊terminal₊i_r),
        "imaginary current" => VIndex(2,:load₊terminal₊i_i),
    )
    # isinteractive() && plottoi(toi)
    @reftest "ZIPLoad_constant_current" toi

    # constant Z
    @named load = ZIPLoad(Pset=-1, Qset=-1, Vset=1,
                          KpZ=1, KpI=0, KpC=0,
                          KqZ=1, KqI=0, KqC=0)
    bus = compile_bus(MTKBus(load))
    toi = bus_on_slack(bus)
    toi["current"] = OrderedDict(
        "current magnitude at bus" => VIndex(2,:busbar₊i_mag),
        "real current" => VIndex(2,:load₊terminal₊i_r),
        "imaginary current" => VIndex(2,:load₊terminal₊i_i),
    )
    # isinteractive() && plottoi(toi) # neither i nor p is constant
    @reftest "ZIPLoad_constant_Z" toi
end

@testset "Classical machine" begin
    @mtkmodel GenBus begin
        @components begin
            machine = Library.ClassicalMachine(;
                τ_m_input=false,
                S_b=100,
                V_b=18,
                ω_b=2π*60,
                X′_d=0.0608,
                R_s=0.000124,
                H=23.64,
            )
            busbar = BusBar()
        end
        @equations begin
            connect(machine.terminal, busbar.terminal)
        end
    end
    @named mtkbus = GenBus()

    bus = compile_bus(mtkbus)
    set_voltage!(bus; mag=1.017, arg=0.0295)
    set_current!(bus; P=0.716, Q=0.3025)
    initialize_component!(bus)

    toi = bus_on_slack(bus; tmax=600, toilength=10_000)
    # isinteractive() && plottoi(toi)
    @reftest "ClassicalMachine" toi
end

@testset "AVR model" begin
    E1 = 3.5461
    E2 = 4.7281
    Se1 = 0.08
    Se2 = 0.26
    se_quad = x -> Library.quadratic_ceiling(x, E1, E2, Se1, Se2)
    Ae, Be = Library.solve_ceilf(E1=>Se1, E2=>Se2)
    se_exp  = x -> Ae* exp(Be*x)

    if isinteractive()
        let fig = Figure()
            ax = Axis(fig[1,1])
            xs = range(0, 6; length=100)
            lines!(ax, xs, se_quad.(xs); label="quad")
            lines!(ax, xs, se_exp.(xs); label="exp")
            axislegend(ax)
            scatter!(ax, [(E1, Se1), (E2, Se2)])
            fig
        end
    end
end
