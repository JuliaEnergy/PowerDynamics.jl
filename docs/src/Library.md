# Component Library

This page documents all the pre-built components available in the PowerDynamics.jl library.

!!! warning "PowerDynamics.Library Under Active Development"
    **The PowerDynamics.Library component library is currently excluded from semantic versioning and is under heavy development.**

    While PowerDynamics itself follows semantic versioning, the Library submodule's API is highly unstable and variable names, function signatures, and model interfaces may change frequently without notice. If you are using specific models from PowerDynamics.Library in their current state, we strongly recommend copying them to your own source code to avoid breaking changes in future updates.

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
SauerPaiMachine
Swing
ClassicalMachine
```

## Control Systems

### Automatic Voltage Regulators (AVRs)
```@docs
AVRFixed
AVRTypeI
```

### Governors and Turbines
```@docs
GovFixed
TurbineGovTypeI
TGOV1
```

## Load Models

### Static Load Models
```@docs
PQLoad
VoltageDependentLoad
ConstantYLoad
ZIPLoad
```

## Line Models

### Transmission Line Models
```@docs
PiLine
PiLine_fault
```

## Fault Models

### Ground Fault Models
```@docs
RXGroundFault
```
