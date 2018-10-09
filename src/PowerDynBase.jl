# (C) 2018 Potsdam Institute for Climate Impact Research, authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)

module PowerDynBase

using Lazy: @>>
using MacroTools
using Markdown # for the @doc
using Parameters: @with_kw

#####################################################################
# helper
begin
    approxzero(val) = isapprox(val, 0; atol=100*eps(Float64), rtol=0)

    macro def(name, definition)
        return quote
            macro $(esc(name))()
                esc($(Expr(:quote, definition)))
            end
        end
    end

    _removelinereferences(b) = b
    function _removelinereferences(b::Expr)
        if ~(b.head == :block)
            return b
        end
        b.args = [ el for el in b.args if ~(el isa LineNumberNode) ]
        return b
    end
    removelinereferences(q::Expr) = MacroTools.postwalk( _removelinereferences , q)
    rlr = removelinereferences
end

#####################################################################



# errors
include("Errors.jl")

# complex view (interpreting part of an array of real values as an array with complex values)
include("complexview.jl")

# DEVariables
include("DEVariables.jl")

# base definitions for all node dynamics
include("NodeParametersBase.jl")
include("NodeDynamicsBase.jl")

# all possible node dynamics
include("DynamicNodeMacro.jl")
include("NodeDynamics/PQAlgebraic.jl")
include("NodeDynamics/PVAlgebraic.jl")
include("NodeDynamics/SlackAlgebraic.jl")
include("NodeDynamics/SwingEquation.jl")
include("NodeDynamics/FourthOrderEquation.jl")


# the structures building the grid dynamics from the node dynamics
include("NetwirkRHSs.jl")
include("GridDynamics.jl")

# States (kind of an interface)
include("States.jl")

# export of the main types and functions
export GridDynamicsError, NodeDynamicsError, MissingParameterError, StateError
export OrdinaryNodeDynamics, OrdinaryNodeDynamicsWithMass
export no_internal_masses, no_internal_differentials
export @DynamicNode, AbstractNodeParameters, showdefinition
export construct_node_dynamics
export NetworkRHS, GridDynamics
export SystemSize, Nodes, AdmittanceLaplacian, masses, differentials
export State
export internalsymbolsof, internaldsymbolsof, internaloutsymbolsof, parametersof

export convert, promote_rule # only so the autodocs work properly

end # module DPSABase
