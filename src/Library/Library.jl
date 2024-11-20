module Library

using ArgCheck: @argcheck
using ModelingToolkit: ModelingToolkit
using ModelingToolkit: @connector, @mtkmodel, @variables, @parameters
using ModelingToolkit: @named, @unpack, ODESystem, Equation, Num, unwrap
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
@connector Terminal begin
    u_r(t), [description="d-voltage"]
    u_i(t), [description="q-voltage"]
    i_r(t), [description="d-current", connect=Flow]
    i_i(t), [description="q-current", connect=Flow]
end

@mtkmodel SystemBase begin
    @parameters begin
        SnRef = 100, [description="System base"]
        fNom = 50, [description="AC system frequency"]
        ωNom = 2 * π * fNom, [description="System angular frequency"]
        ωRef0Pu = 1, [description="Reference for system angular frequency (pu base ωNom)"]
        ω0Pu = 1, [description="System angular frequency (pu base ωNom)"]
   end
end

@mtkmodel PUBase begin
    @parameters begin
        S, [description="Base power in MVA"]
        V, [description="Base voltage in kV"]
        ω, [description="System angular frequency in rad/s"]
        I=S/V, [description="Base current in kA"]
        Z=V^2/S, [description="Base impedance in Ω"]
        Y=S/V^2, [description="Base admittance in S"]
    end
end
Ibase(S, V) = S/V
Zbase(S, V) = V^2/S
Ybase(S, V) = S/V^2

export BusBar, MTKBus, SlackAlgebraic, SlackDifferential
include("Bus.jl")

export LineEnd, MTKLine
include("Lines.jl")

export iscomponentmodel, isbusmodel, isbranchmodel, islinemodel
include("Interfaces.jl")

export pin_parameters
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

####
#### Control Models
####
export AVRTypeI
include("Controls/AVRs.jl")

export TurbineGovTypeI, TGOV1
include("Controls/Govs.jl")

####
#### Load Models
####
export PQLoad, VoltageDependentLoad, ConstantYLoad
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
