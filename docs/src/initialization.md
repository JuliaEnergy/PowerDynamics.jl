# Powergrid Initialization

Initialization of power grid simulations requires a multi-step approach that combines steady-state power flow analysis with dynamic component initialization. OpPoDyn.jl provides a structured framework for this process, building on the initialization capabilities of NetworkDynamics.jl.

## Overview

Power grid initialization involves finding valid initial conditions that satisfy both:
1. **Power flow constraints**: Electrical power balance equations (steady-state)
2. **Dynamic constraints**: Initial conditions for dynamic components (generators, controllers, etc.)

This is typically achieved through a two-step process:
1. Solve the power flow problem to determine steady-state electrical conditions
2. Initialize dynamic components using the power flow solution as boundary conditions

## Multi-Step Initialization Process

```julia
nw = get_dynamic_network(...)

# extract powerflow model       # ⎫                 ⎫
pfnw = powerflow_model(nw)      # │                 │
# initial guess for powerflow   # ⎬ solve_powerflow │
pfs0 = NWState(pfnw)            # │                 │  
# find fixpoint for pf model    # │                 │
pfs = find_fixpoint(pfnw, pfs0) # ⎭                 ⎬ initialize_from_pf[!]
# extract interface (u/i values)#                   │
interf = interface_values(pfs)  #                   │
# initialize around powerflow   #                   │
initialize_componentwise[!](    #                   │
    nw;                         #                   │
    default_overrides = interf  #                   │
)                               #                   ⎭
```
This low-level step-wise interface allows users full control and complete management of the initialization process. However, OpPoDyn.jl also provides higher-level wrapper functions [`solve_powerflow`](@ref) and [`initialize_from_pf`](@ref) that combine these steps for common use cases.

### Step 1: Power Flow Model Extraction

The first step creates a simplified, algebraic representation of the power grid that captures the essential power flow relationships:

```julia
pfnw = powerflow_model(nw)
```

This function extracts the power flow network from the full dynamic network model, creating a steady-state representation. The power flow network itself is also a `Network` in the NetworkDyanamics.jl sense.

The `powerflow_model` function determines the appropriate power flow representation for each dynamic Node and LineModel by:
checking if `:pfmodel` metadata is set, which points to a different component model specifically designed for power flow analysis

If the `:pfmodel` is **not** set, it assumes that the same model is used for both power flow and dynamic simulation. This is the case for purly static models such as PiLines or PQ-Loads.

### Step 2: Power Flow Solution

The power flow problem is solved using NetworkDynamics.jl's [`find_fixpoint`](@exref) function,
which internally uses NonlinearSolve.jl:

```julia
pfs0 = NWState(pfnw)            # Initial guess for power flow state
pfs = find_fixpoint(pfnw, pfs0) # Solve power flow equations
```

This step finds the steady-state solution where:
- Active and reactive power are balanced at each bus
- Generation and load are in equilibrium

### Step 3: Interface Value Extraction

The power flow solution provides boundary conditions for dynamic component initialization:

```julia
interf = interface_values(pfs)
```

This extracts voltage magnitudes, voltage angles, and current flows at each network node, which serve as interface constraints for the dynamic components.

### Step 4: Component-wise Dynamic Initialization

Finally, each dynamic component is initialized individually using the power flow solution as boundary conditions:

```julia
initialize_componentwise!(nw; default_overrides = interf)
```

This step leverages NetworkDynamics.jl's component-wise initialization to determine free internal states and parameters (such as rotor angles or controller setpoints), such that the
**steady state** of the overall network matches the flows from the power flow solution (i.e. all currents and voltages match).

## Advanced Component Initialization

In some cases, the standard initialization process may not be sufficient. For example, when component initialization constraints cannot be expressed solely in terms of **interface variables** (voltages and currents), but need access to other variables from the complete power flow solution.

### Power Flow Dependent Initialization Constraints

OpPoDyn.jl provides [`PFInitConstraint`](@ref) for these advanced scenarios. Unlike regular [`InitConstraint`](@extref)s from NetworkDynamics.jl, PFInitConstraints can access any variable from the solved power flow state, not just interface variables.
I.e. you get access to states, parameters and observables from the power flow model of the same component.

#### When to Use PFInitConstraints

Use `PFInitConstraints` when your component's initialization requires:
- Access to internal variables of other components (e.g., generator power setpoints)
- Constraints involving power flow quantities that aren't interface variables
- Complex initialization logic that depends on the complete network state

#### Syntax and Usage

The [`@pfinitconstraint`](@ref) macro provides convenient syntax for defining these constraints:

```julia
# Single constraint accessing both component and power flow variables
constraint = @pfinitconstraint :dynamicload₊P - @pf(:PQ₊Pset)

# Multiple constraints in a single block
constraints = @pfinitconstraint begin
    :pibranch₊X - @pf(:pibranch₊X) # "copy" parameters from pf
    :P_gen - @pf(:P_load)          # Power balance constraint
    :AVR₊Vset - :busbar₊u_mag      # init controller setpoints
end

# Attach to a component
set_pfinitconstraint!(my_generator, constraints)
```

#### Key Syntax Elements

- `:symbol` - Access component's own variables during initialization
- `@pf :symbol` - Access variables from the solved power flow state
- Multiple constraints can be combined in a `begin...end` block

### Integration with Initialization Process

PFInitConstraints are automatically handled during [`initialize_from_pf[!]`](@ref):

1. **Power flow solution**: The power flow equations are solved first
2. **Constraint specialization**: All `PFInitConstraints`` are converted to regular `InitConstraints`` by "specializing" them with the power flow solution (i.e. the `@pf(:x)` blocks are replaced by the actual values)
3. **Component initialization**: The specialized constraints are passed to NetworkDynamics.jl's component initialization

This process is transparent to the user - simply define your `PFInitConstraints` and use `initialize_from_pf[!]` as usual.
