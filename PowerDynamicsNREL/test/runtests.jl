using Test
using SafeTestsets

@safetestset "MetaGenerator Tests" begin include("MetaGenerator_test.jl") end
@safetestset "Bus Tests" begin include("Bus_test.jl") end
