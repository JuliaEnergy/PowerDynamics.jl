
using Test
using Crayons

testlist = [
    ("nodes/PQAlgebraic.jl", "PQAlgebraic Node Tests"),
    ("nodes/PVAlgebraic.jl", "PVAlgebraic Node Tests"),
    ("nodes/SlackAlgebraic.jl", "SlackAlgebraic Node Tests"),
    ("nodes/SwingEq.jl", "SwingEq Node Tests"),
    ("nodes/SwingEqLVS.jl", "SwingEqLVS Node Tests"),
    ("nodes/FourthOrderEquation.jl", "FourthOrderEq Node Tests"),
    ("nodes/VoltageSourceInverterMinimal.jl", "VoltageSourceInverterMinimal Node Tests")
    #("dynamicnodemacro.jl", "Dynamic Node Macro Tests"),
    #("outputanderrors.jl", "Output and Error Tests"),
    #("states.jl", "States Tests"),
]
@testset "All Tests" begin
    @testset "$desc" for (file, desc) in testlist
        t = @elapsed include(file)
        println(Crayon(foreground = :green, bold = true), "$desc:", Crayon(reset = true), " $t s")
    end
end
