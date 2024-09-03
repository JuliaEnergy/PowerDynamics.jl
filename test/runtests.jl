using PowerDynamicsPrototype
using ModelingToolkit
using ModelingToolkit: t_nounits as t
using Test

@testset "PowerDynamicsPrototype.jl" begin
    @testset "pin parameters" begin
        @mtkmodel Inner begin
            @parameters begin
                a = 0
                b = 0
            end
            @variables begin
                i(t)
            end
            @equations begin
                i ~ a + b
            end
        end
        @mtkmodel Outer begin
            @components begin
                inner = Inner()
            end
            @parameters begin
                a = 0
                b = 0
            end
            @variables begin
                o(t)
            end
            @equations begin
                o ~ inner.i * (a + b)
            end
        end
        @named outer = Outer()
        full_equations(outer)
        sys1 = pin_parameters(outer, :a => 1)
        sys2 = pin_parameters(outer, :outer₊b => 1)
        full_equations(sys1)
        full_equations(sys2)
        sys1 = pin_parameters(outer, outer.inner.a => 1)
        sys2 = pin_parameters(outer, :inner₊a => 1)
        @test full_equations(sys1) == full_equations(sys2)
    end
end
