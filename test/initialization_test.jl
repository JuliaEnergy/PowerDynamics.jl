using PowerDynamics
using NetworkDynamics
using PowerDynamics: @capture, postwalk
using Graphs
using ModelingToolkitBase
using ModelingToolkitBase: t_nounits as t, D_nounits as Dt
using ModelingToolkitStandardLibrary.Blocks: RealInput, RealOutput
using PowerDynamics.Library
using NetworkDynamics: set_mtk_defaults
PowerDynamics.load_pdtesting()
using Main.PowerDynamicsTesting

@info "Start Initialization tests"

@testset "Initialization constraint construction" begin
    ic1 = @pfinitconstraint :x + :y + @pf(:z)
    ic2 = PFInitConstraint([:x, :y], [:z], 1) do out, u, pfu
        out[1] = u[:x] + u[:y] + pfu[:z]
    end

    out1 = [0.0]
    out2 = [0.0]
    u = rand(2)
    pfu = rand(1)
    ic1(out1, u, pfu)
    ic2(out2, u, pfu)
    @test out1 == out2

    ic1 = @pfinitconstraint begin
        :x + :y + @pf :z
        :z^2 - @pf :x
    end
    ic2 = PFInitConstraint([:x, :y, :z], [:x, :z], 2) do out, u, upf
        out[1] = u[:x] + u[:y] + upf[:z]
        out[2] = u[:z]^2 - upf[:x]
    end
    out1 = [0.0, 0.0]
    out2 = [0.0, 0.0]
    u = rand(3)
    pfu = rand(2)
    ic1(out1, u, pfu)
    ic2(out2, u, reverse(pfu))
    @test out1 == out2

    # test with captured runtime variables
    ics = Any[]
    for _scale in [1,2,3]
        ic = @pfinitconstraint :x * _scale - @pf(:z)
        push!(ics, ic)
    end
    for (ic3, _scale) in zip(ics, [1,2,3])
        out3 = [0.0]
        u3 = rand(1)
        pfu3 = rand(1)
        ic3(out3, u3, pfu3)
        @test out3[1] ≈ u3[1] * _scale - pfu3[1]
        @test !occursin("Expr(:escape", repr(ic3))
    end
end

@testset "PFInitFormula construction" begin
    # Test single formula
    if1 = @pfinitformula :Vset = sqrt(:u_r^2 + :u_i^2)
    if2 = PFInitFormula([:Vset], [:u_r, :u_i], Symbol[]) do out, u, pfu
        out[:Vset] = sqrt(u[:u_r]^2 + u[:u_i]^2)
    end

    out1 = [0.0]  # For :Vset
    out2 = [0.0]  # For :Vset
    u = [3.0, 4.0]  # For :u_r, :u_i
    pfu = Float64[]  # Empty for no pf variables

    if1(out1, u, pfu)
    if2(out2, u, pfu)
    @test out1[1] ≈ out2[1] ≈ 5.0

    # Test formula with power flow variable
    if3 = @pfinitformula :Pset = :u_r * :i_r + :u_i * :i_i + @pf(:Pload)
    if4 = PFInitFormula([:Pset], [:u_r, :i_r, :u_i, :i_i], [:Pload]) do out, u, pfu
        out[:Pset] = u[:u_r] * u[:i_r] + u[:u_i] * u[:i_i] + pfu[:Pload]
    end

    out3 = [0.0]  # For :Pset
    out4 = [0.0]  # For :Pset
    u = [1.0, 3.0, 2.0, 4.0]  # For :u_r, :i_r, :u_i, :i_i
    pfu = [5.0]  # For :Pload

    if3(out3, u, pfu)
    if4(out4, u, pfu)
    @test out3[1] ≈ out4[1] ≈ 16.0  # 1*3 + 2*4 + 5 = 16

    # Test multiple formulas in begin/end block
    if5 = @pfinitformula begin
        :Vset = sqrt(:u_r^2 + :u_i^2)
        :Pset = :u_r * :i_r + :u_i * :i_i + @pf(:Pload)
    end

    out5 = [0.0, 0.0]  # For :Vset, :Pset
    u = [1.0, 2.0, 3.0, 4.0]  # For :u_r, :u_i, :i_r, :i_i
    pfu = [5.0]  # For :Pload

    if5(out5, u, pfu)
    @test out5[1] ≈ sqrt(5.0)  # sqrt(1^2 + 2^2) = sqrt(5)
    @test out5[2] ≈ 16.0       # 1*3 + 2*4 + 5 = 16

    # test with captured runtime variables
    _offset = 2.0
    if6 = @pfinitformula :Vset = sqrt(:u_r^2 + :u_i^2) + _offset
    out6 = [0.0]
    u6 = [3.0, 4.0]
    pfu6 = Float64[]
    if6(out6, u6, pfu6)
    @test out6[1] ≈ 5.0 + _offset
    @test !occursin("Expr(:escape", repr(if6))

    _scale = 0.5
    if7 = @pfinitformula :Pset = :u_r * :i_r + _scale * @pf(:Pload)
    out7 = [0.0]
    u7 = [1.0, 3.0]  # u_r, i_r
    pfu7 = [4.0]     # Pload
    if7(out7, u7, pfu7)
    @test out7[1] ≈ 1.0 * 3.0 + _scale * 4.0
    @test !occursin("Expr(:escape", repr(if7))
end

@testset "Test end to end initialziation" begin
    @testset "test happy path" begin
        nw = PowerDynamicsTesting.load_ieee9bus()
        s0_nonmut = initialize_from_pf(nw)
        s0_nonmut_meta = NWState(nw)
        s0_mut = initialize_from_pf!(nw)
        s0_mut_meta = NWState(nw)

        equal_states(a, b) = uflat(a)==uflat(b) && pflat(a) == pflat(b)
        @test equal_states(s0_nonmut, s0_mut)
        @test equal_states(s0_mut, s0_mut_meta)
        @test !equal_states(s0_nonmut, s0_nonmut_meta)
    end

    @testset "test add of identical formula/constraint" begin
        nw = PowerDynamicsTesting.load_ieee9bus()
        eidx = EIndex(4)
        em = nw[eidx]

        pfif = @pfinitformula :pibranch₊active = @pf(:pbranch₊active)
        @test add_pfinitformula!(nw, eidx, pfif)
        @test !add_pfinitformula!(em, pfif)

        pfic = @pfinitconstraint :pibranch₊active - @pf(:pbranch₊active)
        @test add_pfinitconstraint!(em, pfic)
        @test !add_pfinitconstraint!(em, pfic)

        f1 = copy_pf_parameters(em)
        f2 = copy_pf_parameters(em)
        @test f1 == f2
    end

    @testset "test pfinitformula" begin
        nw = PowerDynamicsTesting.load_ieee9bus()
        pfnw = powerflow_model(nw)
        pfs0 = NWState(pfnw)
        pfs0.p.e[9, :pibranch₊active] = 0 # deactivate a line
        pfs = solve_powerflow(pfnw; pfs0)
        @test pfs.p.e[9, :pibranch₊active] == 0

        # cannot initialize because not steadystate
        @test_throws ComponentInitError initialize_from_pf(nw; pfs=pfs)

        s0 = initialize_from_pf(nw; pfs=pfs, nwtol=Inf, tol=Inf)
        # without formula/constraint, the parameter is not initialized
        @test s0.p.e[9, :pibranch₊active] != pfs0.p.e[9, :pibranch₊active]

        # remove the default parameter, initialization can recover
        # the active state from the interface variables
        for i in 1:ne(nw)
            em = nw[EIndex(i)]
            param = only(filter(s -> contains(string(s), "active"), psym(em)))
            delete_default!(em, param)
            set_guess!(em, param, 1.0)
        end
        s1 = initialize_from_pf(nw; pfs=pfs, pfnw=nothing, pfs0=nothing, subverbose=EIndex(9));
        @test s1.p.e[9, :pibranch₊active] ≈ pfs0.p.e[9, :pibranch₊active] atol=1e-10

        # alternative, we can use the copy formula to copy them over
        for i in 1:ne(nw)
            em = nw[EIndex(i)]
            # set default again
            param = only(filter(s -> contains(string(s), "active"), psym(em)))
            set_default!(em, param, 1.0)

            # add copy formula
            form = copy_pf_parameters(em)
            add_pfinitformula!(em, form)
        end

        s2 = initialize_from_pf(nw; pfs=pfs, subverbose=EIndex(9));
        @test s2.p.e[9, :pibranch₊active] == pfs0.p.e[9, :pibranch₊active]
    end

    @testset "pf model with missing defaults" begin
        nw = PowerDynamicsTesting.load_ieee9bus()
        pfnw = powerflow_model(nw)
        delete_default!(pfnw[VIndex(3)], :pv₊P)
        delete_default!(pfnw[VIndex(3)], :busbar₊u_i)

        @test_throws NetworkInitError solve_powerflow(pfnw)
        pfs0 = NWState(pfnw)
        pfs0.v[3, :busbar₊u_i] = 0
        pfs0.p.v[3, :pv₊P] = 0.85
        solve_powerflow(nw; pfs0) # runs now
    end

    @testset "Test handling of keyword overrides in initialize_from_pf" begin
        nw = PowerDynamicsTesting.load_ieee9bus()
        delete_default!(nw, VIndex(1, :generator₊gov₊T1))

        @test_throws ComponentInitError initialize_from_pf(nw) # missing guess
        initialize_from_pf(nw; default_overrides=Dict(VIndex(1,:generator₊gov₊T1)=>0.05))
        s0 = initialize_from_pf(nw; guess_overrides=Dict(VIndex(1,:generator₊gov₊T1)=>0.04))
        @test s0[VIndex(1,:generator₊gov₊T1)] ≈ 0.05 atol=1e-2

        @test_throws ArgumentError initialize_from_pf(nw; additional_initconstraint=:foo)
        @test_throws ArgumentError initialize_from_pf(nw; additional_initformula=:foo)
    end
end

# End-to-end backward initialization through nested components: every model carries its
# own init recipe (`initf`), the AVR pins its PI block's output observable via `set_initf`,
# and composition chains everything into one DAG which determines every free variable from
# the powerflow interface — no nonlinear solve. Mirrors docs/tutorials/nested_component_init.jl.
@testset "nested backward init: initf + set_initf pin" begin
    @component function nbi_pi_block(; name, defaults...)
        @parameters K_p K_i
        @variables begin
            err(t)
            y(t)
        end
        @variables x(t), [guess=0, initf = (y - K_p*err)/K_i]
        sys = System([Dt(x) ~ err, y ~ K_p*err + K_i*x], t; name)
        sys = set_mtk_defaults(sys, defaults)
    end
    @component function nbi_avr(; name, defaults...)
        @named v_mag_in = RealInput()
        @named vf_out = RealOutput()
        @variables begin
            v_meas(t), [guess=1, initf = v_mag_in.u]
            v_f(t), [guess=1]
        end
        @parameters begin
            T_m; T_e; K_e
            v_ref, [guess=1, initf = v_meas]
        end
        @named pi = nbi_pi_block()
        eqs = [
            T_m*Dt(v_meas) ~ v_mag_in.u - v_meas,
            pi.err ~ v_ref - v_meas,
            T_e*Dt(v_f) ~ pi.y - K_e*v_f,
            vf_out.u ~ v_f,
        ]
        sys = System(eqs, t; name, systems=[v_mag_in, vf_out, pi])
        sys = set_initf(sys, pi.y => K_e*v_f)
        sys = set_mtk_defaults(sys, defaults)
    end
    @component function nbi_gov(; name, defaults...)
        @named ω_in = RealInput()
        @named Pm_out = RealOutput()
        @variables begin
            P_m(t), [guess=1]
            P_in(t)
        end
        @parameters begin
            R; T_g
            P_ref, [guess=1, initf = P_m - (1 - ω_in.u)/R]
        end
        eqs = [
            P_in ~ P_ref + (1 - ω_in.u)/R,
            T_g*Dt(P_m) ~ P_in - P_m,
            Pm_out.u ~ P_m,
        ]
        sys = System(eqs, t; name, systems=[ω_in, Pm_out])
        sys = set_mtk_defaults(sys, defaults)
    end
    @component function nbi_machine(; name, defaults...)
        @named terminal = Terminal()
        @named vf_in = RealInput()
        @named Pm_in = RealInput()
        @named v_mag_out = RealOutput()
        @named ω_out = RealOutput()
        @parameters R_s X′_d H
        ## frequency base is taken structurally from the bus's `systembase`
        @parameters ωbase, [bound_to = :systembase₊ωbase]
        @variables begin
            ω(t)=1
            V_d(t); V_q(t); I_d(t); I_q(t); τ_e(t); v_mag(t); e_r(t); e_i(t)
        end
        @variables begin
            δ(t), [guess=0, initf = atan(e_i, e_r)]
            v_f(t), [guess=1, initf = sqrt(e_r^2 + e_i^2)]
            P_m(t), [guess=1, initf = ω*τ_e]
        end
        T_to_loc(α)  = [ sin(α) -cos(α); cos(α)  sin(α)]
        T_to_glob(α) = [ sin(α)  cos(α); -cos(α) sin(α)]
        eqs = [
            v_f ~ vf_in.u
            P_m ~ Pm_in.u
            [terminal.u_r, terminal.u_i] .~ T_to_glob(δ)*[V_d, V_q]
            [I_d, I_q] .~ T_to_loc(δ)*[terminal.i_r, terminal.i_i]
            Dt(δ) ~ ωbase*(ω - 1)
            2*H*Dt(ω) ~ P_m/ω - τ_e
            τ_e ~ (V_q + R_s*I_q)*I_q + (V_d + R_s*I_d)*I_d
            0 ~ V_q + R_s*I_q + X′_d*I_d - v_f
            0 ~ V_d + R_s*I_d - X′_d*I_q
            e_r ~ terminal.u_r + R_s*terminal.i_r - X′_d*terminal.i_i
            e_i ~ terminal.u_i + R_s*terminal.i_i + X′_d*terminal.i_r
            v_mag ~ sqrt(V_d^2 + V_q^2)
            v_mag_out.u ~ v_mag
            ω_out.u ~ ω
        ]
        sys = System(eqs, t; name, systems=[terminal, vf_in, Pm_in, v_mag_out, ω_out])
        sys = set_mtk_defaults(sys, defaults)
    end

    @named machine = nbi_machine(; R_s=0.01, X′_d=0.3, H=5)
    @named avr = nbi_avr(; T_m=0.02, T_e=0.5, K_e=1.0, pi₊K_p=20, pi₊K_i=5)
    @named gov = nbi_gov(; R=0.05, T_g=0.5)
    gen = CompositeInjector([machine, avr, gov]; name=:gen)
    genbus = compile_bus(MTKBus(gen; name=:genbus); vidx=2, pf=pfPV(V=1, P=1))

    # the AVR's set_initf pinned the PI output observable
    @test NetworkDynamics.pinned_obssyms(genbus) == Set([:gen₊avr₊pi₊y])

    slackbus = compile_bus(PowerDynamics.VariableFrequencySlack(name=:slack); vidx=1, pf=pfSlack(V=1))
    line = compile_line(MTKLine(PiLine(; name=:piline, R=0.01, X=0.1, B_src=0, B_dst=0)); src=1, dst=2)
    nw = Network([slackbus, genbus], line)

    s0 = initialize_from_pf!(nw; verbose=false)

    # re-run the (non-mutating) component init on the seeded bus to inspect the log:
    # the whole DAG resolves — pin marked, everything computed, nothing solved
    io = IOBuffer()
    initialize_component(nw[VIndex(2)]; verbose=true, io)
    out = String(take!(io))
    @test occursin("(pinned observable)", out)
    @test occursin("No free variables!", out)

    # hand-computed reference: v_meas = V_pf = 1 at the PV bus, so err = 0 and the PI
    # state is exactly y/K_i = K_e*v_f/K_i
    @test s0[VIndex(2, :gen₊avr₊v_meas)] ≈ 1.0
    @test s0[VIndex(2, :gen₊avr₊v_ref)] ≈ 1.0
    v_f = s0[VIndex(2, :gen₊avr₊v_f)]
    @test s0[VIndex(2, :gen₊avr₊pi₊x)] ≈ 1.0 * v_f / 5
    # governor holds the machine's electrical power at nominal speed
    @test s0[VIndex(2, :gen₊gov₊P_ref)] ≈ s0[VIndex(2, :gen₊gov₊P_m)]

    # flat start: the RHS at the initialized state is zero
    du = similar(uflat(s0))
    nw(du, uflat(s0), pflat(s0), 0.0)
    @test maximum(abs, du) < 1e-10
end
