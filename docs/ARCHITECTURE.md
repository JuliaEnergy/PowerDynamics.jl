# Architecture Overview

In the following the architecture of PowerDynamics with its components is explained in more detail.

## Components of PowerDynamics

The main components of PowerDynamics are 
- nodes
- lines
- PowerGrid
- operationpoint
- faults
- PowerGridSolutions

The data flow between some of the components is illustrated in the figure below:

<img src="figures/PowerDynamics_Architecturemd.svg">

### nodes
The [node component](https://github.com/JuliaEnergy/PowerDynamics.jl/tree/master/src/nodes) includes a library of different standard generators, loads and inverters. The definition of new nodes is simplified with the node macro `@DynamicNode`.

### lines
The [line component](https://github.com/JuliaEnergy/PowerDynamics.jl/tree/master/src/lines) includes a standard library of different line and transformer types, e.g. a simple admittance line (`StaticLine`), the `PiModelLine` and a simple transfomer model based on the PiModel.

### PowerGrid
The [PowerGrid component](https://github.com/JuliaEnergy/PowerDynamics.jl/blob/master/src/common/PowerGrid.jl) is built from nodes and lines. It contains all information about the graph and is used to derive the right-hand-side function with the [NetworkDynamics.jl](https://github.com/JuliaEnergy/PowerDynamics.jl/blob/master/src/common/PowerGrid.jl)-library. The powergrid can also be parsed from a json-file with read_powergrid.

### operationpoint
The [operationpoint](https://github.com/JuliaEnergy/PowerDynamics.jl/tree/master/src/operationpoint) represents the steady-state solution of the dynamic power system. There are different methods for finding the operation point: rootfind, nlsolve and dynamic. The function `find_operationpoint` returns the operationpoint which is of Type State.

### faults
All simulation in PowerDynamics assume some kind of [fault scenario](https://github.com/JuliaEnergy/PowerDynamics.jl/tree/master/src/faults). It can be either a change in initial conditions (different from the operation point of the system) at the beginning of a simulation. Or it is a perturbation of node or line paramers for a certain time span `tspan_fault` that alter the PowerGrid. These can be e.g. a sudden loss of load, `PowerPerturbation`, or a line failure, `LineFailure`. The latter faults are all derived from `AbstractPerturbation`, and all share the same `simulate` function. 

### PowerGridSolution
After a simulation in PowerDynamics the solution returned is of type [PowerGridSolutions](https://github.com/JuliaEnergy/PowerDynamics.jl/blob/master/src/common/PowerGridSolutions.jl). It allows to access the solution output of the ODE solver in a user friendly manner.