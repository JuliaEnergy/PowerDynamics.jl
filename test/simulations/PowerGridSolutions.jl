using Test: @test, @testset
using PowerDynBase: variable_index

import PowerDynBase.symbolsof
import PowerDynBase.dimension

struct Foo
end

symbolsof(f::Foo) = [:a, :b, :c]
dimension(f::Foo) = 3

@testset "PowerGridSolution should find proper variable_index" begin
        idx = variable_index([Foo()], 1, :c)
        @test idx == 3
        idx = variable_index([Foo(),Foo()], 2, :b)
        @test idx == 5
end
