using Test
using PowerDynamics
using PowerDynamics.Library
using PowerDynamics: autoconnections, CompositeInjector
using ModelingToolkitBase
using ModelingToolkitBase: t_nounits as t, D_nounits as Dt
using SciCompDSL: @mtkmodel
using ModelingToolkitStandardLibrary.Blocks

@info "Start modeling_tools tests"

# Helper to test connect equations while string(::Equation) is broken upstream.
# string(con.lhs) returns "connect()" instead of "connect(a, b)", so we search
# the full string for each endpoint individually.
# TODO: replace with direct string equality once upstream bug is fixed.
function is_connection(con, a, b)
    s = string(con)
    contains(s, a) && contains(s, b)
end

@testset "CompositConenctor and autoconnections function" begin
    @testset "basic name matching with suffix stripping" begin
        @mtkmodel TestOutput begin
            @components begin
                signal_meas = RealOutput()
            end
            @variables begin
                y(t) = 1.0
            end
            @equations begin
                signal_meas.u ~ y
            end
        end

        @mtkmodel TestInput begin
            @components begin
                signal_in = RealInput()
            end
            @variables begin
                x(t)
            end
            @equations begin
                x ~ signal_in.u
            end
        end

        @named out_sys = TestOutput()
        @named in_sys = TestInput()

        connections = autoconnections([out_sys, in_sys])

        @test length(connections) == 1
        # Sentinel: captures current broken string repr of connect equations.
        # If this fails, string() was fixed upstream — replace is_connection() with something better
        @test string(connections[1]) == "connect(Model out_sys.signal_meas:\nEquations (1):\n  1 connecting: see equations(expand_connections(out_sys.signal_meas))\nUnknowns (1): see unknowns(out_sys.signal_meas)\n  u(t): Inner variable in RealOutput signal_meas, Model in_sys.signal_in:\nEquations (1):\n  1 connecting: see equations(expand_connections(in_sys.signal_in))\nUnknowns (1): see unknowns(in_sys.signal_in)\n  u(t): Inner variable in RealInput signal_in)"
        @test is_connection(connections[1], "out_sys.signal_meas", "in_sys.signal_in")
    end

    @testset "real components autoconnections" begin
        @named machine = SauerPaiMachine(;
            vf_input=true, τ_m_input=true, S_b=100, V_b=1, Sn=100, Vn=1, ω_b=2π*60, R_s=0.000125, T″_d0=0.01,
            T″_q0=0.01, X_ls=0.01460, X_d=0.1460, X′_d=0.0608, X″_d=0.06, X_q=0.1000, X′_q=0.0969,
            X″_q=0.06, T′_d0=8.96, T′_q0=0.310, H=23.64
        )

        @named avr = AVRTypeI(
            vr_min=-5, vr_max=5, Ka=20, Ta=0.2, Kf=0.063, Tf=0.35, Ke=1, Te=0.314,
            E1=3.3, Se1=0.6602, E2=4.5, Se2=4.2662, tmeas_lag=false
        )

        @named gov = TGOV1(R=0.05, T1=0.05, T2=2.1, T3=7.0, DT=0, V_max=5, V_min=-5)

        connections = autoconnections([machine, avr, gov])

        @test length(connections) == 4

        # Test that expected connections are present
        @test any(is_connection.(connections, "machine.ωout", "gov.ω_meas"))
        @test any(is_connection.(connections, "gov.τ_m", "machine.τ_m_in"))
        @test any(is_connection.(connections, "machine.v_mag_out", "avr.v_mag"))
        @test any(is_connection.(connections, "avr.vf", "machine.vf_in"))
    end

    @testset "error handling" begin
        @mtkmodel UnmatchedSystem begin
            @components begin
                terminal = Terminal()
                weird_input = RealInput()
            end
        end

        @named machine = SauerPaiMachine(;
            vf_input=true, τ_m_input=true, S_b=100, V_b=1, Sn=100, Vn=1, ω_b=2π*60, R_s=0.000125,
            T″_d0=0.01, T″_q0=0.01, H=23.64
        )

        @named unmatched = UnmatchedSystem()

        @test_throws ErrorException autoconnections([machine, unmatched])

        @mtkmodel Matched1 begin
            @components begin
                τ_m_out = RealOutput()
                vf = RealOutput()
            end
        end
        @mtkmodel Matched2 begin
            @components begin
                τ_m_meas = RealOutput()
            end
        end
        @named matched1 = Matched1()
        @named matched2 = Matched2()
        @test_throws ErrorException autoconnections([machine, matched1, matched2])
    end

    @testset "CompositeInjector with autoconnections" begin
        @named machine = SauerPaiMachine(;
            vf_input=true, τ_m_input=true, S_b=100, V_b=1, Sn=100, Vn=1, ω_b=2π*60, R_s=0.000125, T″_d0=0.01,
            T″_q0=0.01, X_ls=0.01460, X_d=0.1460, X′_d=0.0608, X″_d=0.06, X_q=0.1000, X′_q=0.0969,
            X″_q=0.06, T′_d0=8.96, T′_q0=0.310, H=23.64
        )

        @named avr = AVRTypeI(
            vr_min=-5, vr_max=5, Ka=20, Ta=0.2, Kf=0.063, Tf=0.35, Ke=1, Te=0.314,
            E1=3.3, Se1=0.6602, E2=4.5, Se2=4.2662, tmeas_lag=false
        )

        @named gov = TGOV1(R=0.05, T1=0.05, T2=2.1, T3=7.0, DT=0, V_max=5, V_min=-5)

        composite = CompositeInjector([machine, avr, gov])

        @test composite isa System
        @test isinjectormodel(composite)

        # Test that it can be used in MTKBus
        mtkbus = MTKBus(composite)
        @test mtkbus isa System
        compile_bus(mtkbus) # can compile?
    end
end

