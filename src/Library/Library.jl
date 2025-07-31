module Library

using ArgCheck: @argcheck
using ..PowerDynamics: Terminal, BusBase, Ibase
using ModelingToolkit: ModelingToolkit, @named, @mtkmodel, @variables, @parameters, simplify,
                       t_nounits as t, D_nounits as Dt, ODESystem
using ModelingToolkit: @unpack, Equation, Num, System # needed for @mtkmodel?
using ModelingToolkitStandardLibrary.Blocks: RealInput, RealOutput
using NonlinearSolve: NonlinearProblem
using SciMLBase: SciMLBase, solve

@mtkmodel SystemBase begin
    @parameters begin
        SnRef = 100, [description="System base"]
        fNom = 50, [description="AC system frequency"]
        ωNom = 2 * π * fNom, [description="System angular frequency"]
        ωRef0Pu = 1, [description="Reference for system angular frequency (pu base ωNom)"]
        ω0Pu = 1, [description="System angular frequency (pu base ωNom)"]
   end
end

export SlackAlgebraic, SlackDifferential

@mtkmodel SlackAlgebraic begin
    @components begin
        busbar = BusBase()
    end
    @parameters begin
        u_set_r=1, [description="bus d-voltage setpoint"]
        u_set_i=0, [description="bus q-voltage setpoint"]
    end
    @equations begin
        busbar.u_r ~ u_set_r
        busbar.u_i ~ u_set_i
    end
end

@mtkmodel SlackDifferential begin
    @parameters begin
        u_init_r=1, [description="bus d-voltage initial value"]
        u_init_i=0, [description="bus q-voltage initial value"]
    end
    @components begin
        busbar = BusBase(;u_r=u_init_r, u_i=u_init_i)
    end
    @equations begin
        Dt(busbar.u_r) ~ 0
        Dt(busbar.u_i) ~ 0
    end
end

####
#### Machine Models
####
export SauerPaiMachine
include("Machines/SauerPaiMachine.jl")

export Swing
include("Machines/Swing.jl")

export ClassicalMachine
include("Machines/ClassicalMachine.jl")

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
export PiLine
include("Branches/PiLine.jl")
export PiLine_fault
include("Branches/PiLine_fault.jl")


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
