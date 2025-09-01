using PowerDynamics
using NetworkDynamics
using PowerDynamics: @capture, postwalk
using Graphs

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
end

(isinteractive() && @__MODULE__()==Main ? includet : include)("testsystems.jl")
@testset "Test end to end initialziation" begin
    @testset "test happy path" begin
        nw = TestSystems.load_ieee9bus()
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
        nw = TestSystems.load_ieee9bus()
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
        nw = TestSystems.load_ieee9bus()
        pfnw = powerflow_model(nw)
        pfs0 = NWState(pfnw)
        pfs0.p.e[9, :pibranch₊active] = 0 # deactivate a line
        pfs = solve_powerflow(pfnw; pfs0)
        @test pfs.p.e[9, :pibranch₊active] == 0

        # cannot initialize because not steadystate
        @test_throws ErrorException initialize_from_pf(nw; pfs=pfs)

        s0 = initialize_from_pf(nw; pfs=pfs, nwtol=5, tol=5)
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
        nw = TestSystems.load_ieee9bus()
        pfnw = powerflow_model(nw)
        delete_default!(pfnw[VIndex(3)], :pv₊P)
        delete_default!(pfnw[VIndex(3)], :busbar₊u_i)

        @test_throws ArgumentError solve_powerflow(pfnw)
        pfs0 = NWState(pfnw)
        pfs0.v[3, :busbar₊u_i] = 0
        pfs0.p.v[3, :pv₊P] = 0.85
        solve_powerflow(nw; pfs0) # runs now
    end
end
