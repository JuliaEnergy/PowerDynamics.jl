# Custom Node Types

There are two ways to define custom node types:
- equation based modeling using [`ModelingToolkit.jl`](https://github.com/SciML/ModelingToolkit.jl) and [`BlockSystems.jl`](https://github.com/hexaeder/BlockSystems.jl)
- using the `@DynamicNode` macro

## Custom Nodes using `BlockSystems.jl`
`PowerDynmaics.jl` provides a node constructor to build nodes from `BlockSystem.jl` objects based on voltage setting blocks:
```@docs
IONode(blk::IOBlock, parameters::Dict)
```
Check out the `examples/BlockSystems` folder for more examples.

There is another constructor to combine several current injecting blocks to a single node:
```@docs
BusNode
BlockPara
```

### Component library
`PowerDynamics` has a submodule `IOComponents` which contains several predefined building blocks which may be used to build complex nodes.
```@autodocs
Modules = [PowerDynamics.IOComponents]
```

## Custom Nodes using the Node Macro 
To define your own Node Types, use the [`PowerDynamics.@DynamicNode`](@ref) macro. The new node type will be a subtype of [`PowerDynamics.AbstractNode`](@ref).

```@docs
@DynamicNode
```

```@docs
PowerDynamics.MassMatrix
```
