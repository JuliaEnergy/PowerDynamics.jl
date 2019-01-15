
using Test
using Crayons

testlist = [
    ("nodedynamics.jl", "Single Node Tests"),
    ("dynamicnodemacro.jl", "Dynamic Node Macro Tests"),
    ("complexview.jl", "Complex View Tests"),
    ("griddynamics.jl", "Grid Construction Tests"),
    ("outputanderrors.jl", "Output and Error Tests"),
    ("states.jl", "States Tests"),
]
@testset "All Tests" begin
    @testset "$desc" for (file, desc) in testlist
        t = @elapsed include(file)
        println(Crayon(foreground = :green, bold = true), "$desc:", Crayon(reset = true), " $t s")
    end
end
