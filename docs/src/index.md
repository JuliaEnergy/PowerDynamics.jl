```@meta
CurrentModule = PowerDynamics
```

# PowerDynamics

PowerDynamics.jl is a Julia package for modeling and simulating power grid dynamics. It provides a comprehensive framework for analyzing electrical power systems, including synchronous machines, loads, lines, and various control elements. The package is built on top of [NetworkDynamics.jl](https://github.com/JuliaDynamics/NetworkDynamics.jl) and offers both predefined component models and the flexibility to create custom power system components.

## Getting Started

- **[Modeling Concepts](@ref)** - Learn the fundamental concepts behind PowerDynamics modeling
- **[Component Library](@ref)** - Explore the available power system component models
- **[Powergrid Initialization](@ref)** - Understand how to properly initialize power system simulations
- **[API Reference](@ref API)** - Complete function and type documentation

It is also highly recommend to out check the docs on
[NetworkDynamics.jl](https://github.com/JuliaDynamics/NetworkDynamics.jl)
as those explain lots of the underlying functionality and concepts

## Tutorials
- **[Custom Components](@ref custom-bus)** - Creating your own bus models

## Examples
- **[IEEE 9-Bus System](@ref ieee9bus)** - Basic power system simulation
- **[IEEE 39-Bus System Part 1](@ref ieee39-part1)** - Modeling of a larger System 
- **[IEEE 39-Bus System Part 2](@ref ieee39-part2)** - Advanced Initialization Tutorial 
- **[IEEE 39-Bus System Part 3](@ref ieee39-part3)** - Dynamic Simulation and Parameter Tuning 

!!! warning "PowerDynamics.Library Under Active Development"
    **The PowerDynamics.Library component library is currently excluded from semantic versioning and is under heavy development.**

    While PowerDynamics itself follows semantic versioning, the Library submodule's API is highly unstable and variable names, function signatures, and model interfaces may change frequently without notice. If you are using specific models from PowerDynamics.Library in their current state, we strongly recommend copying them to your own source code to avoid breaking changes in future updates.

## Reproducibility

```@raw html
<details><summary>Direct dependencies used for this documentation:</summary>
```

```@example
using Pkg # hide
Pkg.status() # hide
```

```@raw html
</details>
```

```@raw html
<details><summary>Julia Version:</summary>
```

```@example
using InteractiveUtils # hide
versioninfo() # hide
```

```@raw html
</details>
```

```@raw html
<details><summary>Full Manifest:</summary>
```

```@example
using Pkg # hide
Pkg.status(; mode = PKGMODE_MANIFEST) # hide
```

```@raw html
</details>
```

## Funding
Development of this project was in part funded by the *German Federal Ministry for Economic Affairs and Climate Action* as part of the *OpPoDyn*-Project (Project ID 01258425/1, 2024-2027).

```@raw html
<img src="assets/bmwk_logo_en.svg" width="300"/>
```
