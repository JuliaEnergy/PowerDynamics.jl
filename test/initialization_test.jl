using OpPoDyn
using NetworkDynamics
using OpPoDyn: @capture, postwalk

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
        em = nw[EIndex(4)]

        pfif = @pfinitformula :pibranch₊active = @pf(:pbranch₊active)
        @test add_pfinitformula!(em, pfif)
        @test !add_pfinitformula!(em, pfif)

        pfic = @pfinitconstraint :pibranch₊active - @pf(:pbranch₊active)
        @test add_pfinitconstraint!(em, pfic)
        @test !add_pfinitconstraint!(em, pfic)
    end

    @testset "test pfinitformula" begin
        nw = TestSystems.load_ieee9bus()

        em = nw[EIndex(1)]



        pfnw = powerflow_model(nw)

    end
end
