# Powergrid Initialization

Initialization of power grid simulations requires a multi-step approach that combines steady-state power flow analysis with dynamic component initialization. OpPoDyn.jl provides a structured framework for this process, building on the initialization capabilities of NetworkDynamics.jl.

For general background on NetworkDynamics initialization concepts, see the [NetworkDynamics Initialization Guide](@extref Initialization).

## Overview

Power grid initialization involves finding valid initial conditions that satisfy both:
1. **Power flow constraints**: Electrical power balance equations (steady-state)
2. **Dynamic constraints**: Initial conditions for dynamic components (generators, controllers, etc.)

This is typically achieved through a two-step process:
1. Solve the power flow problem to determine steady-state electrical conditions
2. Initialize dynamic components using the power flow solution as boundary conditions

This follows the [Two-Step Initialization Pattern](@extref Component-wise-Network-Initialization) described in NetworkDynamics.jl, specialized for power grid applications.

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

Note: This workflow above is slightly simplified, see [Integration with Initialization Process](@ref) below for the full set of commands.

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

The power flow problem is solved using NetworkDynamics.jl's [`find_fixpoint`](@extref NetworkDynamics.find_fixpoint) function,
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

This step leverages NetworkDynamics.jl's [component-wise initialization](@extref Component-wise-Network-Initialization) to determine free internal states and parameters (such as rotor angles or controller setpoints), such that the
**steady state** of the overall network matches the flows from the power flow solution (i.e. all currents and voltages match).

For details on how component initialization works, see the [Single Component Initialization](@extref Single-Component-Initialization) section in NetworkDynamics.jl.

## Advanced Component Initialization

In some cases, the standard initialization process may not be sufficient. For example, when component initialization constraints cannot be expressed solely in terms of **interface variables** (voltages and currents), but need access to other variables from the complete power flow solution.

NetworkDynamics.jl provides general [InitFormulas and InitConstraints](@extref  Advanced-Component-Initialization:-Formulas-and-Constraints) for advanced component initialization. OpPoDyn.jl extends these concepts with power flow-aware variants that can access the complete power flow solution.

### PFInitConstraints vs PFInitFormulas

| Method | Purpose | Usage |
|--------|---------|-------|
| [`PFInitConstraint`](@ref) | Add constraint equations that must be satisfied | When you need to enforce specific relationships between variables |
| [`PFInitFormula`](@ref) | Set default initial values directly | When you need to initialize variables based on power flow results |

Both methods can access any variable from the solved power flow state, not just interface variables. You get access to states, parameters and observables from the power flow model of the same component.

**Key difference**: Constraints **increase the number of equations** that must be satisfied during initialization, while formulas **reduce the number of free variables** by setting additional default values.

These are power flow-aware extensions of NetworkDynamics.jl's standard [`InitConstraint`](@extref NetworkDynamics.InitConstraint) and [`InitFormula`](@extref NetworkDynamics.InitFormula) mechanisms.

### Power Flow Dependent Initialization Constraints

[`PFInitConstraint`](@ref) adds constraint equations to the initialization problem. Unlike regular [`InitConstraint`](@extref NetworkDynamics.InitConstraint)s from NetworkDynamics.jl, PFInitConstraints can access power flow variables.

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

### Power Flow Dependent Initialization Formulas

[`PFInitFormula`](@ref) sets default initial values for variables using both component and power flow variables. Unlike constraints, formulas directly assign values without adding equations to solve.

The [`@pfinitformula`](@ref) macro provides convenient syntax:

```julia
# Single formula - set variable from component variables
@pfinitformula :Vset = sqrt(:u_r^2 + :u_i^2)

# Formula using power flow variables
@pfinitformula :Pset = @pf(:generator_power)

# Multiple formulas in a block
@pfinitformula begin
    :Vset = sqrt(:u_r^2 + :u_i^2)
    :Pset = @pf(:generator_power)
end

# Attach to a component
set_pfinitformula!(my_generator, formulas)
```

### Integration with Initialization Process

Both PFInitConstraints and PFInitFormulas are automatically handled during [`initialize_from_pf[!]`](@ref initialize_from_pf):

1. **Power flow solution**: The power flow equations are solved first
2. **Specialization**: All `PFInitConstraints` and `PFInitFormulas` are converted to regular `InitConstraints` and `InitFormulas` by "specializing" them with the power flow solution (i.e. the `@pf(:x)` blocks are replaced by the actual values)
3. **Component initialization**: The specialized constraints and formulas are passed to NetworkDynamics.jl's component initialization

This process is transparent to the user - simply define your power flow dependent initialization methods and use `initialize_from_pf[!]` as usual.

The underlying mechanism follows NetworkDynamics.jl's [component initialization pipeline](@extref Single-Component-Initialization), with the power flow solution providing additional context for constraint and formula evaluation.

The extended initializaion workflow (automaticially done in `initialize_from_pf[!]`) looks like this:

```julia
nw = get_dynamic_network(...)
pfnw = powerflow_model(nw)
pfs0 = NWState(pfnw)
pfs = find_fixpoint(pfnw, pfs0)
interf = interface_values(pfs)

# specialize the constaints an formulas and pass then down
pfconstraints = specialize_pfinitconstraints(nw, pfs)
pfformulas    = specialize_pfinitformulas(nw, pfs)
initialize_componentwise[!](
    nw;
    default_overrides = interf,
    additional_initconstraints = pfconstraints,
    additional_initformulas = pfformulas,
)
```
