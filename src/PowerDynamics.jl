# (C) 2018 Potsdam Institute for Climate Impact Research, authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)

module PowerDynamics

using Markdown # for the @doc

include("common/Errors.jl")
include("common/Helpers.jl")
include("common/PowerGrid.jl")
include("common/States.jl")

include("parsers/Format.jl")
include("parsers/JsonParser.jl")

# all possible node dynamics
include("nodes/AbstractNode.jl")
include("nodes/NodeMacro.jl")
include("nodes/PQAlgebraic.jl")
include("nodes/PVAlgebraic.jl")
include("nodes/SlackAlgebraic.jl")
include("nodes/SwingEquation.jl")
include("nodes/FourthOrderEq.jl")
include("nodes/FourthOrderEqGovernorExciterAVR.jl")
include("nodes/VoltageSourceInverterMinimal.jl")
include("nodes/VoltageSourceInverterVoltagePT1.jl")
include("nodes/CurrentSourceInverterMinimal.jl")
include("nodes/ExponentialRecoveryLoad.jl")
include("nodes/experimental/RLCLoad.jl")


# all line types
include("lines/AbstractLine.jl")
include("lines/LineMacro.jl")
include("lines/StaticLine.jl")
include("lines/PiModelLine.jl")
include("lines/PiModel.jl")
include("lines/Transformer.jl")

include("operationpoint/operationpoint.jl")

include("simulations/PowerGridSolutions.jl")
include("simulations/simulations.jl")
include("simulations/NodeShortCircuit.jl")


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

end
