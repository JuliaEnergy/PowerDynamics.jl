module Library

using ArgCheck: @argcheck
using ..PowerDynamics: PowerDynamics, Terminal, BusBase, Ibase
using NetworkDynamics: NetworkDynamics, ComponentCondition, ComponentAffect,
                       VertexModel, VIndex, EIndex, NWState,
                       VectorContinuousComponentCallback, DiscreteComponentCallback, ComponentPostprocessing
using ModelingToolkit: ModelingToolkit, @named, simplify, t_nounits as t, D_nounits as Dt,
                       @component
# needed for @mtkmodel
using ModelingToolkit: @mtkmodel, @variables, @parameters, @unpack, Num, System, Equation, connect, setmetadata
using ModelingToolkitStandardLibrary.Blocks: RealInput, RealOutput
using NonlinearSolve: NonlinearProblem
using SciMLBase: SciMLBase, solve
using Symbolics: Symbolics
using LinearAlgebra: LinearAlgebra
using SparseConnectivityTracer: GradientTracer
using ScopedValues: ScopedValue

export CallbackVerbose, set_callback_verbosity!
export SaturationConfiguration, SaturationConfig, set_saturation_config!
"""
    const CallbackVerbose = ScopedValue(true)

Toggle verbosity of callbacks during simulation. Use [`set_callback_verbosity!`](@ref) to change
globally or use
```julia
with(CallbackVerbose => false) do
    # your code here
end
```
to set temporarily via ScopedValue mechanism.
"""
const CallbackVerbose = ScopedValue(true)
"""
    set_callback_verbosity!(v::Bool)

Sets [`CallbackVerbose`](@ref) to `v`.
"""
function set_callback_verbosity!(v::Bool)
    CallbackVerbose[] = v
end


"""
    simplify_barrier(x) = x

Symbolicially registers a function that acts as a barrier to simplification.
It does nothing nummericially but is opaque to structural simplify / mtkcompile.
Can be used to prevent unwanted simplifications which might lead to devision by zero.
"""
simplify_barrier(x) = x
Symbolics.@register_symbolic simplify_barrier(x)

"""
    @no_simplify a ~ a + b

Macro to prevent simplification of an equation during mtkcompile.
Transforms the equation to an explicit constraint opaque to structural simplification:

    0 ~ simplify_barrier(rhs - lhs)

where `lhs ~ rhs` is the original equation.
"""
macro no_simplify(ex)
    if ex isa Expr && ex.head == :call && ex.args[1] == :~
        lhs, rhs = ex.args[2], ex.args[3]
        return :(0 ~ $(simplify_barrier)($(esc(rhs)) - $(esc(lhs))))
    else
        throw(ArgumentError("@no_simplify can only be used on equations. Can't handle $ex"))
    end
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

####
#### Slack Models
####
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

export SimpleLag, SimpleLead, LeadLag, Derivative, SimpleGain
export SimpleLagLim, LimIntegrator
export DeadZone
export QUAD_SE, EXP_SE
export ss_to_mtkmodel, siso_tf_to_ss
include("building_blocks.jl")

include("Machines/PSSE_BaseMachine.jl")

# Synchronous Machine Models
export PSSE_GENCLS
include("Machines/PSSE_GENCLS.jl")

export PSSE_GENROU, PSSE_GENROE
include("Machines/PSSE_GENROUND.jl")

export PSSE_GENSAL, PSSE_GENSAE
include("Machines/PSSE_GENSALIENT.jl")

export SauerPaiMachine
include("Machines/SauerPaiMachine.jl")

export Swing
include("Machines/Swing.jl")

export ClassicalMachine
include("Machines/ClassicalMachine.jl")

####
#### Control Systems
####

# Exciters & AVRs
export PSSE_EXST1
include("Controls/EX/PSSE_EXST1.jl")

export PSSE_ESST4B
include("Controls/EX/PSSE_ESST4B.jl")

export PSSE_ESST1A
include("Controls/EX/PSSE_ESST1A.jl")

export PSSE_SCRX
include("Controls/EX/PSSE_SCRX.jl")

export PSSE_IEEET1
include("Controls/EX/PSSE_IEEET1.jl")

export AVRFixed, AVRTypeI
include("Controls/AVRs.jl")

# Governors and Turbines
export PSSE_IEEEG1
include("Controls/GOV/PSSE_IEEEG1.jl")

export PSSE_HYGOV
include("Controls/GOV/PSSE_HYGOV.jl")

export GovFixed, TurbineGovTypeI, TGOV1
include("Controls/Govs.jl")

export PSSE_GGOV1_EXPERIMENTAL
include("Controls/GOV/PSSE_GGOV1.jl")

# Power System Stabilizers (PSS)
export PSSE_IEEEST
include("Controls/PSS/PSSE_IEEEST.jl")

####
#### Load Models
####

export PQLoad, VoltageDependentLoad, ConstantYLoad, ZIPLoad, ConstantCurrentLoad
include("Loads/StaticLoads.jl")

# Static Load Models
export PSSE_Load
include("Loads/PSSE_Load.jl")

####
#### Line Models
####

# Transmission Line Models
export PiLine
include("Branches/PiLine.jl")

export PiLine_fault
include("Branches/PiLine_fault.jl")

export Breaker
include("Branches/Breaker.jl")

####
#### Fault Models
####

# Ground Fault Models
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
