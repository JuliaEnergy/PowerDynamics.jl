module Library

using ArgCheck: @argcheck
using ModelingToolkit: ModelingToolkit
using ModelingToolkit: @connector, @mtkmodel, @variables, @parameters
using ModelingToolkit: @named, @unpack, ODESystem, System, Equation, Num, unwrap
using ModelingToolkit: connect, simplify, getname, unknowns, parameters, iscomplete, rename, defaults
using ModelingToolkit: get_name, get_eqs, get_observed, get_ctrls, get_defaults, get_schedule,
                       get_connector_type, get_gui_metadata, get_preface, get_initializesystem,
                       get_continuous_events, get_discrete_events, get_parameter_dependencies, get_iv,
                       get_discrete_subsystems, get_solved_unknowns, get_systems, get_tspan, get_guesses
using ModelingToolkit: t_nounits as t, D_nounits as Dt
using Symbolics: Symbolics, Symbolic, iscall, fixpoint_sub
using ModelingToolkitStandardLibrary.Blocks: RealInput, RealOutput
using NonlinearSolve: NonlinearProblem
using SciMLBase: SciMLBase, solve

export Terminal

"""
    Terminal

A ModelingToolkit connector for electrical terminals in power system components.

Represents an electrical connection point with complex voltage and current in dq coordinates.
The terminal defines the interface between power system components like buses, lines, and machines.

# Variables
- `u_r(t)`: d-axis voltage component
- `u_i(t)`: q-axis voltage component
- `i_r(t)`: d-axis current component (flow variable)
- `i_i(t)`: q-axis current component (flow variable)

# Notes
Current variables are defined as flow variables, meaning they sum to zero at connection points
according to Kirchhoff's current law.

See also: [`BusBar`](@ref), [`LineEnd`](@ref)
"""
@connector Terminal begin
    u_r(t), [description="d-voltage"]
    u_i(t), [description="q-voltage"]
    i_r(t), [guess=0, description="d-current", connect=Flow]
    i_i(t), [guess=0, description="q-current", connect=Flow]
end

Ibase(S, V) = S/V
Zbase(S, V) = V^2/S
Ybase(S, V) = S/V^2

export BusBar, MTKBus, SlackAlgebraic, SlackDifferential, CompositeInjector
include("Bus.jl")

export LineEnd, MTKLine
include("Lines.jl")

export isinjectormodel, isbusmodel, isbranchmodel, islinemodel
include("Interfaces.jl")

include("lib_utils.jl")

####
#### Machine Models
####
export SauerPaiMachine
include("Machines/SauerPaiMachine.jl")

export DynawoMachine
include("Machines/DynawoMachine.jl")

export Swing
include("Machines/Swing.jl")

export IPSLPSATOrder4
include("Machines/IPSLPSAT.jl")

export ClassicalMachine
include("Machines/ClassicalMachine.jl")

export ClassicalMachine_powerfactory
include("Machines/ClassicalMachine_powerfactory.jl")

#export StandardModel_pf
#include("Machines/StandardModel_pf.jl")

export StandardModel_pf_testneu
include("Machines/StandardModel_pf_testneu.jl")

export StandardModel_pf
include("Machines/StandardModel_pf.jl")

#export MarconatoModel
#include("Machines/MarconatoModel.jl")

####
#### Control Models
####
export AVRFixed, AVRTypeI
include("Controls/AVRs.jl")

export GovFixed, TurbineGovTypeI, TGOV1
include("Controls/Govs.jl")

####
#### Load Models
####
export PQLoad, VoltageDependentLoad, ConstantYLoad, ZIPLoad
include("Loads/PQLoad.jl")

####
#### Line Models
####
export DynawoPiLine
include("Branches/DynawoPiLine.jl")
export PiLine
include("Branches/PiLine.jl")
export PiLine_fault
include("Branches/PiLine_fault.jl")

export DynawoFixedRatioTransformer
include("Transformers/DynawoFixedRatioTransformer.jl")

####
#### Fault models
####
export RXGroundFault
include("Faults/Faults.jl")

####
#### Powerflow models
####
using NetworkDynamics: ComponentModel, has_metadata, get_metadata
include("powerflow_models.jl")

end
