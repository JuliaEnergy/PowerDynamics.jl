```@meta
CurrentModule = PowerDynamics
```

# PowerDynamics

PowerDynamics.jl is a Julia package for modeling and simulating power grid dynamics. It provides a comprehensive framework for analyzing electrical power systems, including synchronous machines, loads, lines, and various control elements. The package is built on top of [NetworkDynamics.jl](https://github.com/JuliaDynamics/NetworkDynamics.jl) and offers both predefined component models and the flexibility to create custom power system components.

!!! warning "PowerDynamics.Library Under Active Development"
    **The PowerDynamics.Library component library is currently excluded from semantic versioning and is under heavy development.**

    While PowerDynamics itself follows semantic versioning, the Library submodule's API is highly unstable and variable names, function signatures, and model interfaces may change frequently without notice. If you are using specific models from PowerDynamics.Library in their current state, we strongly recommend copying them to your own source code to avoid breaking changes in future updates.

## Getting Started

- **[Modeling Concepts](@ref)** - Learn the fundamental concepts behind PowerDynamics modeling
- **[Component Library](@ref)** - Explore the available power system component models
- **[Powergrid Initialization](@ref)** - Understand how to properly initialize power system simulations
- **[API Reference](@ref API)** - Complete function and type documentation

It is also highly recommend to out check the docs on
[NetworkDynamics.jl](https://github.com/JuliaDynamics/NetworkDynamics.jl)
as those explain lots of the underlying functionality and concepts

## Tutorials
- **[Custom Components](@ref custom-bus)** - Shows how to implement Milano's classical synchronous machine model with a power system stabilizer (PSS)
- **[Custom Transmission Lines](@ref custom-line)** - Demonstrates creating a PI-branch transmission line model with overcurrent protection that can trip during faults

## Examples
- **[IEEE 9-Bus System](@ref ieee9bus)** - Simulates the complete 9-bus IEEE test system with synchronous generators and dynamic load changes
- **[IEEE 39-Bus System Part 1](@ref ieee39-part1)** - Shows how to build the 39-bus New England test system from custom CSV data files with proper component modeling
- **[IEEE 39-Bus System Part 2](@ref ieee39-part2)** - Demonstrates the detailed initialization process for the 39-bus system including power flow and initialization of dynamic models
- **[IEEE 39-Bus System Part 3](@ref ieee39-part3)** - Runs dynamic simulation of the 39-bus system with a short circuit disturbance and fault clearing
- **[IEEE 39-Bus System Part 4](@ref ieee39-part4)** - Implements a custom droop-controlled inverter model and performs parameter optimization using sensitivity analysis
- **[EMT Toy Model Example](@ref emt-toymodel)** - Demonstrates very basic EMT modeling using dynamic shunt capacitor and RL transmission line components in rotating dq coordinates

## Reproducibility

```@raw html
<details><summary>Direct dependencies used for this documentation:</summary>
```

```@example
using Pkg #hide
Pkg.status() #hide
```

```@raw html
</details>
```

```@raw html
<details><summary>Julia Version:</summary>
```

```@example
using InteractiveUtils #hide
versioninfo() #hide
```

```@raw html
</details>
```

```@raw html
<details><summary>Full Manifest:</summary>
```

```@example
using Pkg #hide
Pkg.status(; mode = PKGMODE_MANIFEST) #hide
```

```@raw html
</details>
```

## Funding
Development of this project was in part funded by the *German Federal Ministry for Economic Affairs and Climate Action* as part of the *OpPoDyn*-Project (Project ID 01258425/1, 2024-2027).

```@raw html
<img src="assets/bmwk_logo_en.svg" width="300"/>
```
