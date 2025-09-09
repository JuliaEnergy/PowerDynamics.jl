module Library

using ArgCheck: @argcheck
using ..PowerDynamics: Terminal, BusBase, Ibase
using ModelingToolkit: ModelingToolkit, @named, simplify, t_nounits as t, D_nounits as Dt
# needed for @mtkmodel
using ModelingToolkit: @mtkmodel, @variables, @parameters, @unpack, Num, System, Equation, connect
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

export SlackAlgebraic, SlackDifferential, VariableFrequencySlack

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

@mtkmodel VariableFrequencySlack begin
    @components begin
        busbar = BusBase()
    end
    @parameters begin
        V, [guess=1, description="bus voltage magnitude"]
        ω = 1, [description="slack frequency in pu (base ωNom)"]
        ω_b=2π*50, [description="System base frequency in rad/s"]
    end
    @variables begin
        δ(t), [guess=0, description="voltage angle"]
    end
    @equations begin
        Dt(δ) ~ ω_b*(ω - 1)
        busbar.u_r ~ V * cos(δ)
        busbar.u_i ~ V * sin(δ)
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
#### Renewables
####
export IdealDroopInverter
include("Renewables/IdealDroopInverter.jl")

####
#### Powerflow models
####
include("powerflow_models.jl")

end
