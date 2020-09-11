using Test
using Crayons

testlist = [
    ("operationpoint/test_example.jl", "PowerModels operationpoint Tests"),
    ("common/States.jl", "States Tests"),
    ("common/PowerGrid.jl", "PowerGrid Tests"),
    ("common/PowerGridSolutions.jl", "PowerGridSolutions Tests"),
    ("nodes/PQAlgebraic.jl", "PQAlgebraic Node Tests"), 
    ("nodes/VoltageDependentLoad.jl", "VoltageDependentLoad Node Tests"),
    ("nodes/PVAlgebraic.jl", "PVAlgebraic Node Tests"),
    ("nodes/SlackAlgebraic.jl", "SlackAlgebraic Node Tests"),
    ("nodes/SwingEq.jl", "SwingEq Node Tests"),
    ("nodes/SwingEqLVS.jl", "SwingEqLVS Node Tests"),
    ("nodes/FourthOrderEq.jl", "FourthOrderEq Node Tests"),
    ("nodes/VoltageSourceInverterMinimal.jl", "VoltageSourceInverterMinimal Node Tests"),
    ("nodes/VoltageSourceInverterVoltagePT1.jl", "VoltageSourceInverterVoltagePT1 Node Tests"),
    ("nodes/CurrentSourceInverterMinimal.jl", "CurrentSourceInverterMinimal Node Tests"),
    ("nodes/experimental/CompositeNode.jl", "CompositeNode Node Tests"),
    ("nodes/ExponentialRecoveryLoad.jl", "ExponentialRecoveryLoad Node Tests"),
    ("nodes/FourthOrderEqGovernorExciterAVR.jl", "FourthOrderEqGovernorExciterAVR Node Tests"),
    ("nodes/experimental/RLCLoad.jl", "RLCLoad Node Tests"),
    ("nodes/experimental/PVInverterWithFrequencyControl.jl","PVInverter Tests"),
    ("nodes/experimental/WindTurbineGenType4.jl","WindGenType 4 Tests"),
    ("nodes/experimental/WindTurbineGenType4_RotorControl.jl","WindGenType 4 with Rotor Control Tests"),
    ("nodes/experimental/CurtailedPowerPlantWithInertia.jl","CurtailedPowerPlantWithInertia Tests"),
    ("lines/StaticLine.jl", "StaticLine tests"),
    ("lines/PiModelLine.jl", "PiModelLine tests"),
    ("lines/RLLine.jl", "RLLine tests"),
    ("lines/Transformer.jl", "Transformer tests"),
    ("faults/AbstractPerturbation.jl", "AbstractPerturbation Tests"),
    ("faults/ChangeInitialConditions.jl", "ChangeInitialConditions Tests"),
    ("faults/NodeParameterChange.jl", "NodeParameterChange Tests"),
    ("faults/LineFailure.jl", "LineFailure Tests"),
    ("faults/PowerPerturbation.jl", "PowerPerturbation Simulation Tests"),
    ("faults/NodeShortCircuit.jl", "NodeShortCircuit Simulation Tests"),
    ("operationpoint/operationpoint.jl", "operationpoint Tests"),
    ("parsers/JsonParser.jl", "JsonParser Tests")
]
@testset "All Tests" begin
    @testset "$desc" for (file, desc) in testlist
        t = @elapsed include(file)
        println(Crayon(foreground = :green, bold = true), "$desc:", Crayon(reset = true), " $t s")
    end
end
