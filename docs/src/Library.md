# Component Library

This page documents all the pre-built components available in the PowerDynamics.jl library.

!!! warning "PowerDynamics.Library Under Active Development"
    **The PowerDynamics.Library component library is currently excluded from semantic versioning and is under heavy development.**

    While PowerDynamics itself follows semantic versioning, the Library submodule's API is highly unstable and variable names, function signatures, and model interfaces may change frequently without notice. If you are using specific models from PowerDynamics.Library in their current state, we strongly recommend copying them to your own source code to avoid breaking changes in future updates.

## Building Blocks

The following building blocks can be used to construct custom control systems and machine models.

### Basic Blocks
```@docs
SimpleLag
SimpleLead
LeadLag
Derivative
SimpleGain
SimpleLagLim
LimIntegrator
DeadZone
```

### Utility Functions
```@docs
ss_to_mtkmodel
siso_tf_to_ss
```

### Saturation Functions
```@docs
QUAD_SE
EXP_SE
```

## Slack Models

### Algebraic and Differential Slack
```@docs
SlackAlgebraic
SlackDifferential
VariableFrequencySlack
```

## Machine Models

### Synchronous Machine Models
```@docs
PSSE_GENCLS
PSSE_GENROU
PSSE_GENROE
PSSE_GENSAL
PSSE_GENSAE
SauerPaiMachine
Swing
ClassicalMachine
```

## Control Systems

### Exciters & AVRs
```@docs
PSSE_EXST1
PSSE_ESST4B
PSSE_ESST1A
PSSE_SCRX
PSSE_IEEET1
AVRFixed
AVRTypeI
```

### Governors and Turbines
```@docs
PSSE_IEEEG1
PSSE_HYGOV
GovFixed
TurbineGovTypeI
TGOV1
PSSE_GGOV1_EXPERIMENTAL
```

### Power System Stabilizers (PSS)
```@docs
PSSE_IEEEST
```

## Load Models

### Static Load Models
```@docs
PQLoad
VoltageDependentLoad
ConstantYLoad
ZIPLoad
ConstantCurrentLoad
PSSE_Load
```

## Line Models

### Transmission Line Models
```@docs
PiLine
PiLine_fault
Breaker
DynamicRLBranch
```

## Shunt Models

### Static and Dynamic Shunts
```@docs
StaticShunt
DynamicRCShunt
```

## Fault Models

### Ground Fault Models
```@docs
RXGroundFault
```
