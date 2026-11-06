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
            vf_input=true, τ_m_input=true, R_s=0.000125, T″_d0=0.01,
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
            vf_input=true, τ_m_input=true, R_s=0.000125,
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
            vf_input=true, τ_m_input=true, R_s=0.000125, T″_d0=0.01,
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

@testset "compile_bus/compile_line wrap bare injectors/branches" begin
    @named swing = Swing(Pm=1.0)
    @named piline = PiLine(R=0.01, X=0.1)

    # a bare injector/branch is wrapped in MTKBus/MTKLine, giving the same model
    # as the explicit two-step spelling — name included
    @test compile_bus(swing).psym == compile_bus(MTKBus(swing)).psym
    @test compile_bus(swing).name == compile_bus(MTKBus(swing)).name == :bus
    @test compile_line(piline).psym == compile_line(MTKLine(piline)).psym
    @test compile_line(piline).name == compile_line(MTKLine(piline)).name == :line
    # an explicit name still wins, in both spellings
    @test compile_bus(swing; name=:mybus).name == :mybus
    @test compile_line(piline; name=:myline).name == :myline

    # kwargs still apply to the wrapped model
    @test get_default(compile_bus(swing; Vbase=110), :busbar₊Vbase) == 110
    @test get_default_from(compile_bus(swing; current_source=true), :busbar₊Vbase) ==
        (:hub, :busbar₊Vbase)
    # without an explicit `Vbase`, the weak fallback is *materialized* into a real default —
    # `default_from` reads defaults in a pre-pass and cannot see a formula, so a bus that
    # kept the fallback would silently fail to feed its line ends
    @test get_default(compile_bus(swing), :busbar₊Vbase) == get_Vbase()
    @test !has_initformula(compile_bus(swing))
    # ... and neither spelling of the voltage base may be given twice
    @test_throws ArgumentError compile_bus(swing; Vbase=110, var"busbar₊Vbase"=220)

    # a system which satisfies neither interface is still rejected
    @named neither = RealOutput()
    @test_throws ArgumentError compile_bus(neither)
    @test_throws ArgumentError compile_line(neither)
end

@testset "MTKBus with extra equations keeps their bus level variables" begin
    @named load = PQLoad(Pset=-1, Qset=0)
    @variables extra(t)
    @parameters kp=2.0
    bus = MTKBus([load], [extra ~ kp * load.terminal.u_r]; name=:bus)
    @test any(isequal(extra), ModelingToolkitBase.get_unknowns(bus))
    @test any(isequal(kp), ModelingToolkitBase.get_ps(bus))
    @test :kp ∈ compile_bus(bus).psym
end


import NetworkDynamics
@testset "global bases with the model cache" begin
    # a reused (cached) model must still pick up the global bases that are current at build time
    cacheext = Base.get_extension(NetworkDynamics, :NetworkDynamicsMTKExt)
    cache_hits() = sum(c -> sum(e -> e.hits, values(c.entries); init=0), cacheext.ALL_CACHES)
    function two_bus()
        gen = compile_bus(MTKBus(SauerPaiMachine(; name=:machine, vf_input=false, τ_m_input=false));
            vidx=1, pf=pfSlack(V=1))
        load = compile_bus(MTKBus(ZIPLoad(; name=:load)); vidx=2, pf=pfPQ(P=-1, Q=-0.1))
        line = compile_line(MTKLine(PiLine(; name=:piline)); src=1, dst=2)
        Network([gen, load], [line])
    end
    function with_bases(f, Sbase, fbase)
        set_Sbase!(Sbase); set_fbase!(fbase)
        try f() finally set_Sbase!(); set_fbase!() end
    end
    base_syms(c) = filter(s -> occursin("base", string(s)), vcat(psym(c), obssym(c)))
    function base_values(nw)
        s = NWState(nw)
        vcat([s[VIndex(i, sym)] for i in 1:length(nw.im.vertexm) for sym in base_syms(nw[VIndex(i)])],
             [s[EIndex(i, sym)] for i in 1:length(nw.im.edgem) for sym in base_syms(nw[EIndex(i)])])
    end

    nw_ref = NetworkDynamics.with_mtk_cache(:off) do
        with_bases(two_bus, 200, 60)
    end
    NetworkDynamics.with_mtk_cache(:on) do
        nw_default = two_bus()
        hits = cache_hits()
        nw_changed = with_bases(two_bus, 200, 60)
        @test cache_hits() > hits # second build was served from the cache

        @test NWState(nw_default)[VIndex(1, :machine₊Sbase)] == 100
        @test NWState(nw_changed)[VIndex(1, :machine₊Sbase)] == 200
        @test NWState(nw_changed)[VIndex(1, :machine₊ωbase)] ≈ 2π*60
        @test isequal(pflat(NWParameter(nw_changed)), pflat(NWParameter(nw_ref)))
        @test isequal(base_values(nw_changed), base_values(nw_ref))
        @test !isequal(base_values(nw_default), base_values(nw_ref))
    end
end

# cached components whose bodies capture global settings, through submodels they build
NetworkDynamics.@cached_component function WrappedLoadBus(; name)
    MTKBus(ZIPLoad(; name=:load); name)
end
NetworkDynamics.@cached_component function WrappedLimInt(; name)
    LimIntegrator(; name, K=1, outMin=-1, outMax=1)
end
@testset "global settings are part of the component cache key" begin
    eqstrings(sys) = string.(ModelingToolkitBase.equations(sys))
    rhs_ref = with_saturation_config(() -> WrappedLimInt(; name=:lim, cache=:off), :rhs_hard)
    NetworkDynamics.with_mtk_cache(:on) do
        NetworkDynamics.wipe_mtk_cache!()
        a = compile_bus(WrappedLoadBus(; name=:a))
        set_Sbase!(200); set_Vbase!(220)
        b = try compile_bus(WrappedLoadBus(; name=:b)) finally set_Sbase!(); set_Vbase!() end
        @test get_default(a, :systembase₊Sbase) == 100
        @test get_default(b, :systembase₊Sbase) == 200
        @test get_default(b, :busbar₊Vbase) == 220

        cb = WrappedLimInt(; name=:lim)
        rhs = with_saturation_config(() -> WrappedLimInt(; name=:lim), :rhs_hard)
        @test eqstrings(rhs) == eqstrings(rhs_ref)
        @test eqstrings(rhs) != eqstrings(cb)
    end
end
