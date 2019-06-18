
using Test
using Crayons

testlist = [
    ("nodes/PQAlgebraic.jl", "PQAlgebraic Node Tests"),
    ("nodes/PVAlgebraic.jl", "PVAlgebraic Node Tests"),
    ("nodes/SlackAlgebraic.jl", "SlackAlgebraic Node Tests"),
    ("nodes/SwingEq.jl", "SwingEq Node Tests"),
    ("nodes/SwingEqLVS.jl", "SwingEqLVS Node Tests"),
    ("nodes/FourthOrderEq.jl", "FourthOrderEq Node Tests"),
    ("nodes/VoltageSourceInverterMinimal.jl", "VoltageSourceInverterMinimal Node Tests"),
    ("nodes/VoltageSourceInverterVoltagePT1.jl", "VoltageSourceInverterVoltagePT1 Node Tests"),
    ("nodes/CurrentSourceInverterMinimal.jl", "CurrentSourceInverterMinimal Node Tests"),
    ("nodes/ExponentialRecoveryLoad.jl", "ExponentialRecoveryLoad Node Tests"),
    ("nodes/FourthOrderEqGovernorExciterAVR.jl", "FourthOrderEqGovernorExciterAVR Node Tests"),
    ("lines/LineMacro.jl", "LineMacro tests"),
    ("common/PowerGrid.jl", "PowerGrid Tests"),
    ("common/states.jl", "States Tests"),
]
@testset "All Tests" begin
    @testset "$desc" for (file, desc) in testlist
        t = @elapsed include(file)
        println(Crayon(foreground = :green, bold = true), "$desc:", Crayon(reset = true), " $t s")
    end
end
