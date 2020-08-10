# PowerDynBase.jl Changelog

## Version 2.3.0

* ![enhancement](https://img.shields.io/badge/PD-enhancement-%23a2eeef.svg) [added a voltage dependent load model]()
* ![enhancement](https://img.shields.io/badge/PD-enhancement-%23a2eeef.svg) [added the feature to find operationpoints using SteadyStateProblem]()
* ![enhancement](https://img.shields.io/badge/PD-enhancement-%23a2eeef.svg) [added new implementation of NodeShortCircuit](https://github.com/JuliaEnergy/PowerDynBase.jl/pull/93)
* ![enhancement](https://img.shields.io/badge/PD-enhancement-%23a2eeef.svg) [added dynamic RL Line]()



## Version 0.8

* ![enhancement](https://img.shields.io/badge/PD-enhancement-%23a2eeef.svg) [added dynamic RL Line]()
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
