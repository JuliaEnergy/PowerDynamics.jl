# (C) 2018 Potsdam Institute for Climate Impact Research, authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)

module PowerDynamics

using Markdown # for the @doc

using DataFrames

include("common/Errors.jl")
include("common/Helpers.jl")
include("common/PowerGrid.jl")
include("common/States.jl")
include("common/PowerGridSolutions.jl")

include("parsers/Format.jl")
include("parsers/JsonParser.jl")

# all possible node dynamics
include("nodes/controller/PIControl.jl")
include("nodes/AbstractNode.jl")
include("nodes/NodeMacro.jl")
include("nodes/PQAlgebraic.jl")
include("nodes/VoltageDependentLoad.jl")
include("nodes/PVAlgebraic.jl")
include("nodes/SlackAlgebraic.jl")
include("nodes/SwingEquation.jl")
include("nodes/FourthOrderEq.jl")
include("nodes/FourthOrderEqGovernorExciterAVR.jl")
include("nodes/VoltageSourceInverterMinimal.jl")
include("nodes/VoltageSourceInverterVoltagePT1.jl")
include("nodes/CurrentSourceInverterMinimal.jl")
include("nodes/ExponentialRecoveryLoad.jl")
include("nodes/FourthOrderEqExciterIEEEDC1A.jl")
include("nodes/FourthOrderEqGovernorIEEEG1.jl")
include("nodes/ExciterSaturtionEq.jl")
include("nodes/experimental/RLCLoad.jl")
include("nodes/experimental/PVInverterWithFrequencyControl.jl")
include("nodes/experimental/LinearGenerator.jl")
include("nodes/experimental/LinearPTO.jl")
include("nodes/experimental/HydraulicPTO.jl")
include("nodes/experimental/DirectDrivePTO.jl")
include("nodes/experimental/WindTurbineGenType4.jl")
include("nodes/experimental/WindTurbineGenType4_RotorControl.jl")
include("nodes/experimental/CurtailedPowerPlantWithInertia.jl")
include("nodes/experimental/CompositeNode.jl")
include("nodes/experimental/FluctuationNode.jl")
include("nodes/experimental/NormalForm.jl")

# requirements for the IONodes
include("IONodes/utils.jl")
include("IONodes/IONode.jl")
include("IONodes/BusNode.jl")
include("IONodes/IOComponents.jl")
include("IONodes/GFI_MTK.jl")
include("IONodes/MTK_Load.jl")

module ModularInverter
include("IONodes/MIComponents.jl")
include("IONodes/ModularInverter.jl")
include("IONodes/LTI.jl")
end

# all line types

include("lines/AbstractLine.jl")
include("lines/PiModel.jl")
include("lines/LineMacro.jl")
include("lines/StaticLine.jl")
include("lines/PiModelLine.jl")
include("lines/Transformer.jl")
include("lines/RLLine.jl")

include("operationpoint/operationpoint.jl")
include("operationpoint/find_valid_initial_condition.jl")
include("operationpoint/power_flow.jl")

include("faults/ChangeInitialConditions.jl")
include("faults/AbstractPerturbation.jl")
include("faults/LineFailure.jl")
include("faults/NodeParameterChange.jl")
include("faults/PowerPerturbation.jl")
include("faults/NodeShortCircuit.jl")

export AbstractNode

# export of the main types and functions
export PowerDynamicsError,NodeDynamicsError,StateError,GridSolutionError,OperationPointError
export no_internal_masses
export @DynamicNode, showdefinition
export construct_vertex

export State

export @Line, StaticLine
export construct_edge

export convert, promote_rule # only so the autodocs work properly

export find_operationpoint

export PowerGrid
export PowerGridSolution
export rhs
export systemsize
export symbolsof
export total_current
export LinearGenerator
export HydraulicPTO

end
