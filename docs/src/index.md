```@meta
CurrentModule = PowerDynamics
```

# PowerDynamics

PowerDynamics.jl is a Julia package for modeling and simulating power grid dynamics. It provides a comprehensive framework for analyzing electrical power systems, including synchronous machines, loads, lines, and various control elements. The package is built on top of [NetworkDynamics.jl](https://github.com/PIK-ICoNe/NetworkDynamics.jl) and offers both predefined component models and the flexibility to create custom power system components.

## Getting Started

- **[Modeling Concepts](@ref)** - Learn the fundamental concepts behind PowerDynamics modeling
- **[Component Library](@ref)** - Explore the available power system component models
- **[Initialization](@ref)** - Understand how to properly initialize power system simulations
- **[API Reference](@ref)** - Complete function and type documentation

It is also highly recommend to out check the docs on
[NetworkDynamics.jl](https://github.com/PIK-ICoNe/NetworkDynamics.jl)
as those explain lots of the underlying functionality and concepts

## Tutorials
- **[Custom Components](@ref)** - Creating your own bus models

## Examples
- **[IEEE 9-Bus System](@ref)** - Basic power system simulation
- **[IEEE 39-Bus System Part 1](@ref)** - Large-scale system modeling
- **[IEEE 39-Bus System Part 2](@ref)** - Advanced analysis techniques

!!! warning "PowerDynamics.Library Under Active Development"
    **The PowerDynamics.Library component library is currently excluded from semantic versioning and is under heavy development.**

    While PowerDynamics itself follows semantic versioning, the Library submodule's API is highly unstable and variable names, function signatures, and model interfaces may change frequently without notice. If you are using specific models from PowerDynamics.Library in their current state, we strongly recommend copying them to your own source code to avoid breaking changes in future updates.

## Funding
Development of this project was in part funded by the *German Federal Ministry for Economic Affairs and Climate Action* as part of the *OpPoDyn*-Project (Project ID 01258425/1, 2024-2027).

```@raw html
<img src="assets/bmwk_logo_en.svg" width="300"/>
```
