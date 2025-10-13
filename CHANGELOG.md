# PowerDynamices.jl Changelog

## Version 4.3.0 Changelog
- [#232](https://github.com/JuliaDynamics/PowerDynamics.jl/pull/232): Added lots of new models based on the great OpenIPSL Library.
  Those models are validated against OpenIPSL by reproducing their component test-harness and comparing trajectories of internal
  states.

## Version 4.2.0 Changelog
- [#230](https://github.com/JuliaDynamics/PowerDynamics.jl/pull/230): 
  - deprecate `Bus(...)` → `compile_bus(...)` and `Line(...)` → `compile_line(...)`
  - remove `PowerDynamicsTesting` as separate package and just load it as module (less env hassle)
  - add new `asciiart` code style for documentation

## Version 4.1.0 Changelog
- [#221](https://github.com/JuliaDynamics/PowerDynamics.jl/pull/221) update for ModelingToolkit.jl v10 compatibility:
  - Update minimum ModelingToolkit.jl requirement from to v10
  - Update minimum NetworkDynamics.jl requirement from v0.10.3 to v0.10.4
  - Remove internal `pin_parameters` function - was never part of public API
  - Rename all `ODESystem` → `System` throughout codebase (follows MTK v10 API)
  - Replace `structural_simplify` with `mtkcompile` for model compilation
  - Replace custom `_to_zero()` hack with `implicit_output()` from NetworkDynamics
  - Update custom metadata system to use `CustomMetadata` wrapper for safer metadata handling

## Version 4.0.0 - Major Breaking Release
In Q2 2024, we began a complete rewrite of PowerDynamics.jl, bringing it much closer in alignment with the modern SciML stack. This rewrite heavily leverages ModelingToolkit.jl for equation-based models and includes a vastly modernized version of our backend, NetworkDynamics.jl.

The 4.0.0 update incorporates all of these changes. We consider the modeling concepts and simulation tools stable enough for release. The library, however, is marked as experimental for now and may change in upcoming minor versions — but you can copy the model definitions into your own code if you rely on an “old” model.

The new version has improved significantly in terms of modeling, initialization, and solution analysis. However, some models and tools previously available are not yet available. If you want to continue using the (unmaintained) old version, stick with PowerDynamics@v3.

--------------------------
## Version 2.4.0

* ![enhancement](https://img.shields.io/badge/PD-enhancement-%23a2eeef.svg) [added Architecture.md](https://github.com/JuliaEnergy/PowerDynamics.jl/issues/52)
* ![enhancement](https://img.shields.io/badge/PD-enhancement-%23a2eeef.svg) [using github actions for CI]()

## Version 2.3.3

* ![enhancement](https://img.shields.io/badge/PD-enhancement-%23a2eeef.svg) [refactored the fault design](https://github.com/JuliaEnergy/PowerDynamics.jl/issues/87)
* ![enhancement](https://img.shields.io/badge/PD-enhancement-%23a2eeef.svg) [added OrderDicts for nodes and lines definition](https://github.com/JuliaEnergy/PowerDynamics.jl/issues/86)

## Version 2.3.0

* ![enhancement](https://img.shields.io/badge/PD-enhancement-%23a2eeef.svg) [added a voltage dependent load model](https://github.com/JuliaEnergy/PowerDynamics.jl/pull/109)
* ![enhancement](https://img.shields.io/badge/PD-enhancement-%23a2eeef.svg) [added the feature to find operationpoints using SteadyStateProblem](https://github.com/JuliaEnergy/PowerDynamics.jl/pull/97)
* ![enhancement](https://img.shields.io/badge/PD-enhancement-%23a2eeef.svg) [added new implementation of NodeShortCircuit](https://github.com/JuliaEnergy/PowerDynBase.jl/pull/93)
* ![enhancement](https://img.shields.io/badge/PD-enhancement-%23a2eeef.svg) [added dynamic RL Line](https://github.com/JuliaEnergy/PowerDynamics.jl/pull/96)
* ![enhancement](https://img.shields.io/badge/PD-enhancement-%23a2eeef.svg) [added composition of node types](https://github.com/JuliaEnergy/PowerDynamics.jl/pull/75)

## Version 0.8

* ![bugfix](https://img.shields.io/badge/PD-bugfix-%23d73a4a.svg) [fixed state conversion](https://github.com/JuliaEnergy/PowerDynBase.jl/pull/62)
* ![bugfix](https://img.shields.io/badge/PD-bugfix-%23d73a4a.svg) [Voltage measurement of the exciter in grid reference frame](https://github.com/JuliaEnergy/PowerDynBase.jl/pull/63)
* ![bugfix](https://img.shields.io/badge/PD-bugfix-%23d73a4a.svg) [removed pi-model function from Transformer.jl and added export statement](https://github.com/JuliaEnergy/PowerDynamics.jl/pull/26)

## Version 0.7

* ![enhancement](https://img.shields.io/badge/PD-enhancement-%23a2eeef.svg) [added CHANGELOG.md and check whether CHANGELOG.md has been modified to ci/travis and proper splitting of ci on travis](https://github.com/JuliaEnergy/PowerDynBase.jl/pull/36)
* ![bugfix](https://img.shields.io/badge/PD-bugfix-%23d73a4a.svg) [symbolsof is now defined on the class level via the @DynamicNode Macro](https://github.com/JuliaEnergy/PowerDynBase.jl/pull/35)
* ![bugfix](https://img.shields.io/badge/PD-bugfix-%23d73a4a.svg) & ![enhancement](https://img.shields.io/badge/PD-enhancement-%23a2eeef.svg) [add Julia 1.1. to travis/ci and fixed wrong coverage reporting (thus)](https://github.com/JuliaEnergy/PowerDynBase.jl/pull/38)
* ![enhancement](https://img.shields.io/badge/PD-enhancement-%23a2eeef.svg) [cleaning .travis.yml](https://github.com/JuliaEnergy/PowerDynBase.jl/pull/39)
* ![enhancement](https://img.shields.io/badge/PD-enhancement-%23a2eeef.svg) [fixing node docs](https://github.com/JuliaEnergy/PowerDynBase.jl/pull/46)
* ![miscellaneous](https://img.shields.io/badge/PD-miscellaneous-lightgrey.svg) [add Sabine Auer to AUTHORS](https://github.com/JuliaEnergy/PowerDynBase.jl/pull/45)
* ![bugfix](https://img.shields.io/badge/PD-bugfix-%23d73a4a.svg) [temporarily pin MbedTLS version to 0.6.6 to fix travis ci issue](https://github.com/JuliaEnergy/PowerDynBase.jl/pull/49)
* ![bugfix](https://img.shields.io/badge/PD-bugfix-%23d73a4a.svg) [adding and package update before the tests in travis to ensure the latest packges are used](https://github.com/JuliaEnergy/PowerDynBase.jl/pull/50)
* ![enhancement](https://img.shields.io/badge/PD-enhancement-%23a2eeef.svg) [add current source inverter to NodeDynamics](https://github.com/JuliaEnergy/PowerDynBase.jl/pull/52)
* ![enhancement](https://img.shields.io/badge/PD-enhancement-%23a2eeef.svg) [Add Exponential Recovery Load Model](https://github.com/JuliaEnergy/PowerDynBase.jl/pull/54)
* ![enhancement](https://img.shields.io/badge/PD-enhancement-%23a2eeef.svg) [Fourth Order Equation with Exciter, AVR and Governor](https://github.com/JuliaEnergy/PowerDynBase.jl/pull/53)
* ![enhancement](https://img.shields.io/badge/PD-enhancement-%23a2eeef.svg) [`make` has now an `open-docs` target that builds the docs and then opens it](https://github.com/JuliaEnergy/PowerDynBase.jl/pull/55)
* ![bugfix](https://img.shields.io/badge/PD-bugfix-%23d73a4a.svg) [removing accidentally added build file `.DS_Store`](https://github.com/JuliaEnergy/PowerDynBase.jl/pull/49)
* ![bugfix](https://img.shields.io/badge/PD-bugfix-%23d73a4a.svg) [adding include statement for PiModel.jl](https://github.com/JuliaEnergy/PowerDynamics.jl/pull/29)
* ![bugfix](https://img.shields.io/badge/PD-bugfix-%23d73a4a.svg) [enabling simulations without slack bus](https://github.com/JuliaEnergy/PowerDynamics.jl/pull/34)
* ![bugfix](https://img.shields.io/badge/PD-bugfix-%23d73a4a.svg) [docs fixed, inertia constant in SM model corrected](https://github.com/JuliaEnergy/PowerDynamics.jl/pull/37)
