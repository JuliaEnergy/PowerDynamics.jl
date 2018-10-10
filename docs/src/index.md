
# PowerDynamics.jl - Dynamic Power System Analysis in Julia

This package provides all the tools you need to create a dynamic power grid model
and analyze it.

## Installation

For now install all the packages using git directly.
**Please be aware of the order and make sure you have access rights.**
They will be registered with the official package managers upon publishing.

TBD: add new installation instructions

## Usage

Generally, we distinguish three types of user for `PowerDynamics.jl`:
- [Grid Modeler](@ref)
- [Grid Component Developer](@ref)
- [`PowerDynamics.jl` Developer](@ref)

### Grid Modeler

**Your Goal** is to use `PowerDynamics.jl` to model your grid of preference. You don't
want to implement new types of nodes.

We recommend you to choose your favorite example from [`PowerDynamicsExamples`](https://gitlab.com/JuliaEnergy/PowerDynamicsExamples),
read [Node Types](@ref) and try to understand it. That should give you the kickstart you need. If you
have any questions, [contact us](@ref Contact).

### Grid Component Developer

**Your Goal** is to use `PowerDynamics.jl` to develop types of nodes, e.g. new control schemes for inverters or
new descriptions of synchronous machines.

After going through the introduction for a [Grid Modeler](@ref), we recommend that you read
through [Node Dynamics Types](@ref) and [Custom Node Types](@ref) and try to implement
a new node type for an example grid. With that, you should have all the tools you need.
If you have any questions, [contact us](@ref Contact).

### `PowerDynamics.jl` Developer

**Your Goal** is to extend `PowerDynamics.jl` with new fundamental functionalities.

After going throught the introduction for a [Grid Modeler](@ref) and a [Grid Component Developer](@ref),
read through the code where hopefully all of this documentation will helpful for you.
Afterwards, it's probably best to open an issue explainng the idea you want to implement
and we can discuss how you can transform your idea into a pull request.
