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
