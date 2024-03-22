#=
```@meta
CurrentModule = PowerDynamics.ModularInverter
```
# Modular Inverter Framework
This section describes the modular inverter framework developed in the MARiE-Project funded by the
*German Federal Ministry for Economic Affairs and Climate Action*.

## Structure

The inverter model is split in two parts: the inner loop and the outer loop.

The *inner loop* either acts as a voltage source
```
           +---------+
u_d_ref -->|         |
u_q_ref -->| Voltage |--> u_d
i_d_dis -->| Source  |--> u_q
i_q_dis -->|         |
           +---------+
```
which gets a `dq` voltage reference and tries to achive this reference despite the `dq` current disturbance.

Alternatively, the inner loop can mimic a current source, which tries to realsize a given current setpoint despite a disturbance voltage at the output terminals:
```
           +---------+
i_d_ref -->|         |
i_q_ref -->| Current |--> i_d
u_d_dis -->| Source  |--> i_q
u_q_dis -->|         |
           +---------+
```

The *outer loop* provides references to the inner loop blocks.
As an input, outerloops have acces to the measured currents and voltages at the point of common coupling.
They either eject a dq-reference for the inner voltage source or the inner current soruce:

```
            +---------+
i_d_meas -->|         |
i_q_meas -->|  Outer  |--> u_d_ref or i_d_ref
u_d_meas -->|  Loop   |--> u_q_ref or i_q_ref
u_q_meas -->|         |
            +---------+
```

`PowerDynamics.jl` provides predefined outer loop as well as inner loop models as `IOBlocks` provide by `BlockSystems.jl`, which can be used to create `IONodes` (see [Custom-Nodes-using-BlockSystems.jl](@ref)).

## Example
In order to construct a modular inverter we need to import the `ModularInverter` module to the namespace:
=#
using PowerDynamics
using BlockSystems
using PowerDynamics.ModularInverter
nothing #hide
#=
First, we need to define the inner loop as an `IOBlock` object.
You can either use your own or make use of the [Predefined Inner Loops](@ref).
=#
BlockSystems.WARN[] = false # hide
inner = CascadedVoltageControl()
#=
As you can see, the input/output structure of that innerloop adheres to the volten source innerloop convention.

Next you define the outer loop, either create you own based on the interface or pick from the [Predefined Outer Loops](@ref).
=#
outer = DroopControl()
#=
There are still a lot of open parameters, we need to assign specific values to create the `IONode`(@ref) later.
We can do so by using the keyword aguments of the  constructor:
=#
outer = DroopControl(;
    ω_ref = 2π*50, # frequency reference
    V_ref = 1,     # voltage reference (pu)
    P_ref = 1,     # active power reference (pu)
    Q_ref = 0,     # reactive power reference (pu)
    K_P   = 0.4,   # active droop gain
    K_Q   = 0.004, # reactive droop gain
    τ_P   = 0.1,   # active power filter constant
    τ_Q   = 0.1,   # reactive power filter constant
)

#=
The full inverter model is created by combining outer and inner model:
=#
inv = Inverter(inner, outer)
#=
Which satisfies the interface for `IONodes`. To create the node
=#
node = IONode(inv);

#=
Which is a valid node definition which can be used in any `PowerDynamics.jl` context.
=#

#=
## Model Reduction

`PowerDynamics.jl` also implements automated model order reduction of linear models using the balance residualization technique.
Internaly, we make use of the method implementaion in `ControlSystems.jl`. For that, a linear block (such as the cascaded inner loops) can be decomposed in to its `A`, `B`, `C` and `D` matrix for regular LTI representation.

For example, we could creat a 10th order representation of the previously defined inner loop by calling
=#
reduced_inner = balresid(inner, 10; reconstruct=false)
#=
To obtain the reduced order model we can construct the inverter again:
=#
reduced_inverter = Inverter(reduced_inner, outer)
# ... and the reduced node
reduced_node = IONode(reduced_inverter);

#=
```@docs
balresid
```
=#

#=
## Predefined Inner Loops
```@docs
CascadedVoltageControl
CascadedVoltageControlCS
CascadedCurrentControl
PT1Source
IdealSource
IdealCSource
```

## Predefined Outer Loops
```@docs
FixedVoltage
DroopControl
Synchronverter
PLLCurrent
ConstantPower
FiltConstantPower
```

## Inverter Construction
```@docs
Inverter
InverterCS
```

=#
