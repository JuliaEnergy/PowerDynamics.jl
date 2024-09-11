module Library

using ArgCheck: @argcheck
using ModelingToolkit: ModelingToolkit
using ModelingToolkit: @connector, @mtkmodel, @variables, @parameters
using ModelingToolkit: @named, @unpack, ODESystem, Equation, Num, unwrap
using ModelingToolkit: connect, simplify, getname, unknowns, parameters, iscomplete, rename, get_name
using ModelingToolkit: t_nounits as t, D_nounits as Dt
using Symbolics: Symbolics, Symbolic, iscall, fixpoint_sub
using ModelingToolkitStandardLibrary.Blocks: RealInput

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

export BusBar, BusModel, SlackAlgebraic, SlackDifferential
include("Bus.jl")

export LineEnd, LineModel
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

####
#### Line Models
####
export DynawoPiLine
include("Lines/DynawoPiLine.jl")

export DynawoFixedRatioTransformer
include("Transformers/DynawoFixedRatioTransformer.jl")

end
