# (C) 2018 Potsdam Institute for Climate Impact Research, authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)

module PowerDynBase

using Markdown # for the @doc

# errors
include("Errors.jl")

include("Helpers.jl")

include("NodeDynamicsBase.jl")

include("PowerGrid.jl")
include("parsers/csv_parser.jl")

# all possible node dynamics
include("DynamicNodeMacro.jl")
include("NodeDynamics/PQAlgebraic.jl")
include("NodeDynamics/PVAlgebraic.jl")
include("NodeDynamics/SlackAlgebraic.jl")
include("NodeDynamics/SwingEquation.jl")
include("NodeDynamics/FourthOrderEquation.jl")
include("NodeDynamics/FourthOrderEquationGovernorExciterAVR.jl")
include("NodeDynamics/VoltageSourceInverterMinimal.jl")
include("NodeDynamics/VoltageSourceInverterVoltagePT1.jl")
include("NodeDynamics/CurrentSourceInverterMinimal.jl")
include("NodeDynamics/ExponentialRecovery.jl")

# all line types
include("Lines/LineMacro.jl")
include("Lines/StaticLine.jl")
include("Lines/PiModelLine.jl")

include("States.jl")

include("operationpoint/operationpoint.jl")

include("simulations/PowerGridSolutions.jl")
include("simulations/simulations.jl")


# export of the main types and functions
export GridDynamicsError, NodeDynamicsError, MissingParameterError, StateError
export no_internal_masses
export @DynamicNode, showdefinition
export construct_node_dynamics

export State

export @Line, StaticLine
export construct_edge

export convert, promote_rule # only so the autodocs work properly

export find_operationpoint

export read_network_from_csv
export PowerGrid
export ode_function
export systemsize
export symbolsof

end
