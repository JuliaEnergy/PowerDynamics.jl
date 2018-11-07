# Dynamics Types

## Node Dynamics

In this section, the implemented general types of node dynamics (which are all subtypes of [`PowerDynBase.AbstractNodeDynamics`](@ref))
with the corresponding helper functions and constants are introduced.
The documentation is done for each type below.
The main types are:
- [`PowerDynBase.OrdinaryNodeDynamics`](@ref)
- [`PowerDynBase.OrdinaryNodeDynamicsWithMass`](@ref)

```@docs
PowerDynBase.AbstractNodeDynamics
PowerDynBase.OrdinaryNodeDynamics
PowerDynBase.OrdinaryNodeDynamicsWithMass
```

### Helper Functions

```
to be done
```

```@docs
PowerDynBase.internalsymbolsof
```

## Grid Dynamics

Analogously, for each of the node dynamics types exists a corresponding
grid dynamics type, that represents the dynamics of a whole power grid model.
They are all subtypes of [`PowerDynBase.GridDynamics`](@ref)
These are:
- [`PowerDynBase.OrdinaryGridDynamics`](@ref)
- [`PowerDynBase.OrdinaryGridDynamicsWithMass`](@ref)
- [`PowerDynBase.AlgebraicGridDynamics`](@ref)

```@docs
PowerDynBase.GridDynamics
PowerDynBase.OrdinaryGridDynamics
PowerDynBase.OrdinaryGridDynamicsWithMass
PowerDynBase.AlgebraicGridDynamics
```

### NetworkRHS

The logic of building a full power grid model (i.e. a subtype of [`PowerDynBase.GridDynamics`](@ref)) is
encoded in [`PowerDynBase.NetworkRHS`](@ref). Currently, the docs are a bit thin here.

```@docs
PowerDynBase.NetworkRHS
```

### Helper Functions

```
to be done
```
