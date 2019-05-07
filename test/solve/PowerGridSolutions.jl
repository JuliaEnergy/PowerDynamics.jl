using PowerDynBase
import PowerDynBase.symbolsof
import PowerDynBase.dimension
using NetworkDynamics
using LightGraphs
using Test

struct Foo
end

symbolsof(f::Foo) = [:a, :b, :c]
dimension(f::Foo) = 3

function tests()
        @testset "PowerGridSolution should find proper variable_index" begin
                idx = PowerDynBase.variable_index([Foo()], 1, :c)
                @test idx == 3
                idx = PowerDynBase.variable_index([Foo(),Foo()], 2, :b)
                @test idx == 5
        end
end

tests()
