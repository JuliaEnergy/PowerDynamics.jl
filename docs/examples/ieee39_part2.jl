#=
# [IEEE39 Bus Tutorial - Part II: Initialization](@id ieee39-part2)

This tutorial can be downloaded as a normal Julia script [here](@__NAME__.jl). #md

This is the second part of a four-part tutorial series for the IEEE 39-bus test system:

- **Part I: Model Creation** - Build the network structure with buses, lines, and components
- **Part II: Initialization** (this tutorial) - Perform power flow calculations and dynamic initialization
- **Part III: Dynamic Simulation** - Run time-domain simulations and analyze system behavior
- **Part IV: Advanced Modeling & Parameter Optimization** - Create custom components and optimize system parameters

The goal of this tutorial is to get an understanding of the initialization process in PowerDynamics.jl.

For comprehensive documentation on initialization, see:
- [NetworkDynamics.jl initialization docs](@extref NetworkDynamics initialization-guide)
- [PowerDynamics.jl initialization docs](@ref "Advanced Component Initialization")

!!! tip "Quick Start"
    If you're looking for the practical initialization approach without diving into implementation details,
    jump directly to the [Initialize all Components](@ref initialize-all-components) section at the end.

This tutorial goes deep into the initialization internals for educational purposes. In practice, 
you'll typically use the high-level functions shown at the end rather than the detailed step-by-step 
process demonstrated here.

As a prerequisite, we load part I of the tutorial, which contains the network model:
=#

using PowerDynamics
EXAMPLEDIR = joinpath(pkgdir(PowerDynamics), "docs", "examples")
include(joinpath(EXAMPLEDIR, "ieee39_part1.jl"))
nw # nw object now available

#=
The initialization process in PowerDynamics.jl is a two-step process: first we
solve the power flow, then we use the power flow results to initialize the
individual network components.

There are shortcut functions to do this (as shown [later](@ref initialize-all-components)), 
but we will go through the steps in detail for educational purposes.

## Power Flow
To solve the power flow, we first need to get the power flow model.
We can use the function [`powerflow_model`](@ref).
=#
pfnw = powerflow_model(nw)
#=
The power flow model is a `Network` object like the original network.
It is built from the original network, by calling `powerflow_model` on the
individual components. For example, we have this rather complex dynamic model at bus 30:
=#
nw[VIndex(30)]
#=
From the printout, we can see that a power flow model `:pvbus` is attached to the generator model.

We can extract the attached PV power flow model by calling `powerflow_model`
=#
powerflow_model(nw[VIndex(30)])
#=
The function `powerflow_model` checks if there is a power flow model attached (it checks the `:pfmodel` [metadata](@extref NetworkDynamics Metadata)).

Per component, the `powerflow_model` function will do the following:
1. If the model has the `:pfmodel` metadata set (see [Metadata](@extref)), it will return the `VertexModel` stored in the metadata. In [Part I](@ref ieee39-part1) we set the `:pfmodel` metadata using the `pf` keyword to the [`compile_bus`].
2. If the model **does not have** the `:pfmodel` metadata set, PowerDynamics will check if the model itself is a valid power flow model. If so, it'll just use the dynamic model as the power flow model.

**What is a valid power flow model?**:
A valid power flow model is a model that has no internal dynamics,
i.e. it either has `dim(model) == 0` OR it has a zero mass matrix (i.e. only constraints).

For example, the PiLine models are completely static:
=#
nw[EIndex(1)]
#=
This model has no internal states and no `:pfmodel` metadata.
=#
@assert ispfmodel(nw[EIndex(1)]) # is pf model itself
powerflow_model(nw[EIndex(1)]) === nw[EIndex(1)]

#=
As a result, when we call `powerflow_model(nw)`, we get a **completely static** network, i.e. only constraints, no dynamics.
=#
all(iszero, pfnw.mass_matrix)

#=
The fully static network has the form
```math
\dot{x} = 0 = f_{\mathrm{nw}}(x, p, t)
```
Where `x` are the network states (mainly voltages $u_r$ and $u_i$ at the buses) and `p` are all the parameters such as $P$, $V$ and $Q$ values for bus models and line parameters for branch models.

To solve this root-finding problem, we first need to find an initial guess.
The initial guess is prefilled with all the default values from the components.
In our case, all parameters and states have default values attached, so the
initial state is fully determined.
=#
pfs0 = NWState(pfnw)
#=
With the default state, we can call [`find_fixpoint`](@extref `NetworkDynamics.find_fixpoint`), which keeps $p$ constant and tries to find an $x$ such that
the root-finding problem stated above is fulfilled:
=#
pfs = find_fixpoint(pfnw, pfs0)
#=
As a result, we get a `NWState` object which contains the **full state** for the power flow model.

Since power flow model and dynamic model share the same topology and the same *network interface* (i.e. nodes create voltages, edges create currents), we can extract the interface values from the power flow model state and apply them to the dynamic model:
=#
interf = interface_values(pfs)
#=
The interface values give us the **inputs** and **outputs** for every component in the network, i.e. for all buses we get values for
- `busbar₊i_r` and `busbar₊i_i`: current input
- `busbar₊u_r` and `busbar₊u_i`: voltage output

For all branches we get:
- `src₊u_r` and `src₊u_i`: source side voltage input
- `dst₊u_r` and `dst₊u_i`: destination side voltage input
- `src₊i_r` and `src₊i_i`: source side current output
- `dst₊i_r` and `dst₊i_i`: destination side current output

With those interface values fixed, we can go over to the second step of the initialization process: the initialization of the dynamic component.

## Initialization of Dynamic Components

Initialization of a bus model means, we want to find a **steady state** of the dynamics given the interface values from the power flow model.

Therefore, the equations for the bus models (recall from [Modeling Concepts](@ref)) become
```math
\begin{aligned}
M_{\mathrm v}\,\frac{\mathrm{d}}{\mathrm{d}t}x_{\mathrm v} = \color{red}{0} &= f^{\mathrm v}\left(x^{\mathrm v}, \color{red}{\sum_k\begin{bmatrix}i^k_r\\ i^k_i\end{bmatrix}}, p_{\mathrm v}, t\right)\\
\color{red}{\begin{bmatrix}u_r\\ u_i\end{bmatrix}} &= g^{\mathrm v}(x^\mathrm{v},p_{\mathrm v}, t)
\end{aligned}
```
where red symbols are fixed by either the power flow solution or our steady state condition (i.e. $\dot{x}=0$).
This leaves us with a system of $N=\mathrm{dim}(x) + \mathrm{dim}(u)$ equations -- we can solve for $N$ unknowns.

!!! details "Edge Initialization Equations"
    ```math
    \begin{aligned}
    M_{\mathrm e}\,\frac{\mathrm{d}}{\mathrm{d}t}x_{\mathrm e} = \color{red}{0} &= f_{\mathrm e}\left(x_{\mathrm e}, \color{red}{\begin{bmatrix} u_r^\mathrm{src}\\u_i^\mathrm{src}\end{bmatrix}}, \color{red}{\begin{bmatrix} u_r^\mathrm{dst}\\u_i^\mathrm{dst}\end{bmatrix}},p_\mathrm{e}, t\right)\\
    \color{red}{\begin{bmatrix}i_r^\mathrm{src}\\i_i^\mathrm{src}\end{bmatrix}} &= g^\mathrm{src}_{\mathrm e}\left(x_{\mathrm e},\color{red}{ \begin{bmatrix} u_r^\mathrm{src}\\u_i^\mathrm{src}\end{bmatrix}}, \color{red}{\begin{bmatrix} u_r^\mathrm{dst}\\u_i^\mathrm{dst}\end{bmatrix}}, p_\mathrm{e}, t\right)\\
    \color{red}{\begin{bmatrix}i_r^\mathrm{dst}\\i_i^\mathrm{dst}\end{bmatrix}} &= g^\mathrm{dst}_{\mathrm e}\left(x_{\mathrm e}, \color{red}{\begin{bmatrix} u_r^\mathrm{src}\\u_i^\mathrm{src}\end{bmatrix}}, \color{red}{\begin{bmatrix} u_r^\mathrm{dst}\\u_i^\mathrm{dst}\end{bmatrix}}, p_\mathrm{e}, t\right)\\
    \end{aligned}
    ```

Notably, the unknowns can come from either the set of states or the set of parameters.
We divide the set of symbols into two sets:
- **fixed** symbols have a `default` metadata set. They are considered fixed in the solution of the nonlinear system.
- **free** symbols only have a `guess` metadata set. They are considered free in the solution of the nonlinear system.

Let's take a look at the bus model at bus 30:
=#
gen = nw[VIndex(30)]
#=
We have a system with 15 States and 42 parameters. Of those, most are fixed, i.e. have a `default` metadata set.
In the VertexModel printout, defaults are shown with `=` while guesses are shown with `≈`.
We can use [`dump_initial_state`](@extref `NetworkDynamics.dump_initial_state`) to get an overview of the free and set states:
=#
dump_initial_state(gen; obs=false)
#=
Right now, we have 2 free parameters, 2 free inputs, 2 free outputs and 15 free states.
We can use [`initialize_component`](@extref `NetworkDynamics.initialize_component`) to find values for the free symbols:
=#
try #hide #md
initialize_component(gen)
catch e #hide #md
    @error e.msg #hide #md
end #hide #md
#=
Wait! The initialization failed! Why? Well, we need to apply additional defaults, so called `default_overwrites` for the inputs/outputs to make the system solvable.
=#
interf_v30 = Dict( # manually define interface values for demonstration
    :busbar₊u_r => 1.04573,
    :busbar₊u_i => -0.0609188,
    :busbar₊i_i => 1.53174,
    :busbar₊i_r => -2.30145,
)
initialize_component(gen; default_overrides=interf_v30)
nothing #hide #md
#=
!!! note "Mutating vs non-mutating Initialization"
    The non mutating [`initialize_component`](@extref `NetworkDynamics.initialize_component`) returns a dictionary
    containing a full initialized component state. That can be useful for certain purposes, often it is easier to work with the *mutating* version of initialization function, which will write the initialized values back to the component metadata (i.e setting the `:init` property for the symbols).

Let's call the mutating `initialize_component!` function to write the initialized values back to the component and inspect the initial state using `dump_initial_state`:
=#
initialize_component!(gen; default_overrides=interf_v30)
dump_initial_state(gen; obs=false)

#=
In the state dump we see how the initialization successfully set all the previously unknown values, including the control parameters `avr₊vref` and `gov₊p_ref`. We have our first initialized component!

However, in practice it's not always so easy.

## Handling Structurally Underconstrained Components
Recalling from [Part 1](@ref ieee39-part1), we have an uncontrolled machine
together with a load on bus 39.

```
            ╔════════════════════════════════╗
            ║ Unctr. Ma. Load Bus (compiled) ║
            ║  ┌────────────────────────┐    ║
  Network   ║  │MTKBus      ┌─────────┐ │    ║
 interface  ║  │          ┌─┤ Machine │ │    ║
  current ────→│ ┌──────┐ │ └─────────┘ │    ║
            ║  │ │BusBar├─o             │    ║
  voltage ←────│ └──────┘ │ ┌──────┐    │    ║
            ║  │          └─┤ Load │    │    ║
            ║  │            └──────┘    │    ║
            ║  └────────────────────────┘    ║
            ╚════════════════════════════════╝
```

If we try to initialize this component as before, we run into a problem:
=#
interf_v39 = Dict(
  :busbar₊u_r => 1.01419,
  :busbar₊u_i => -0.179795,
  :busbar₊i_i => -1.72223,
  :busbar₊i_r => 0.720135,
)
try #hide #md
initialize_component!(nw[VIndex(39)]; default_overrides=interf_v39)
catch e #hide #md
    @error e.msg #hide #md
end #hide #md

#=
Even though we set the interface values, the problem is still underconstrained!
Let's check the free symbols:
=#
println("free u: ", free_u(nw[VIndex(39)]))
println("free p: ", free_p(nw[VIndex(39)]))
nothing #hide #md
#=
We see 8 free states and 3 free parameters, however we only have 8 state + 2 output equations:
=#
nw[VIndex(39)] #hide #md

#=
Even though we have enough set parameters to initialize machine and load on its own, we cannot do it simultaneously.
Intuitively speaking, it's just not clear for the solver which of the two components provides how much power.

To solve this, we have essentially 3 methods:

### Method 1: Manual setting of defaults
The simplest solution is to manually set more defaults. For example, we know that we want to initialize the ZIP load around the initialization point, i.e. $V_\mathrm{set}$ should be the same as the bus voltage magnitude.
=#
vm_manual = copy(nw[VIndex(39)])
u_r = get_initial_state(vm_manual, :busbar₊u_r)
u_i = get_initial_state(vm_manual, :busbar₊u_i)
set_default!(vm_manual, :ZIPLoad₊Vset, sqrt(u_r^2 + u_i^2))
initialize_component!(vm_manual)
nothing #hide #md
#=
The initialization succeeded now!

!!! note "No more default_overrides"
    Note how we can skip the `default_overrides` keyword argument, since the first (failing) call of `initialize_component!` already "burned in" the default overrides! Mutating state is a powerful tool, but it needs care!

### Method 2: Adding an `init_formula`
The problem with the previous method is that it is quite manual.
In reality, we would never go through this very manual initialization process.
The fact that the model is structurally underconstrained is a property of the model and should therefore be handled by the model. To do so, NetworkDynamics.jl provides the [`InitFormula`](@extref `NetworkDynamics.InitFormula`) mechanism.

An `InitFormula` is a symbolic formula that is evaluated during the initialization process. It is attached to the VertexModel so it can be
evaluated automatically during the initialization process.
The "formula" we want to apply is simply the equation
```math
V_\mathrm{set} = \sqrt{u_r^2 + u_i^2}
```
=#
vm_formula = copy(nw[VIndex(39)])
formula = @initformula :ZIPLoad₊Vset = sqrt(:busbar₊u_r^2 + :busbar₊u_i^2)
set_initformula!(vm_formula, formula)
vm_formula #hide #md
#=
The printout shows 1 additional initialization equation was attached to the model.

The initialization works now:
=#
initialize_component!(vm_formula)
nothing #hide #md
#=
The init formula is applied early in the initialization process, essentially writing a new `default` for `ZIPLoad₊Vset` based on the other defaults.

This reduced the number of free variables to 10, thus the system was solvable.

### Method 3: Using an `InitConstraint`
Sometimes, your additional initialization needs are more complicated.
Similar to defining a *formula*, which is evaluated *before* the actual initialization, NetworkDynamics provides a mechanism for injecting additional
constraints into the initialization process.

In contrast to the formula, the constraint does not need to be explicitly solvable, as it defines a residual equation
```math
0 = c(x) = V_\mathrm{set} - \sqrt{u_r^2 + u_i^2}
```
=#
vm_constraint = copy(nw[VIndex(39)])
constraint = @initconstraint begin
  :ZIPLoad₊Vset - sqrt(:busbar₊u_r^2 + :busbar₊u_i^2)
end
set_initconstraint!(vm_constraint, constraint)
vm_constraint #hide #md
#=
With this added constraint, the initialization process is solvable again, since we now have 11 equations for the 11 free variables.
=#
initialize_component!(vm_constraint)
nothing #hide #md
#=
For this particular case, method (2) is the way to go. However there are cases where the constraint is more complex and cannot be expressed as a formula.

See NetworkDynamics docs on [Advanced Component Initialization: Formulas and Constraints](@extref NetworkDynamics) and the PowerDynamics specific extension
[Advanced Component Initialization](@ref) for more information on method 2 and 3.

## [Automatic Initialization of Full Network](@id initialize-all-components)

Let's return from our excursion into individual component initialization and focus on the whole network again.
As we've just seen, we have structurally underconstrained components in the network.
Let's define the init formulas for the two buses which have loads and machines:
=#
formula = @initformula :ZIPLoad₊Vset = sqrt(:busbar₊u_r^2 + :busbar₊u_i^2)
set_initformula!(nw[VIndex(31)], formula)
set_initformula!(nw[VIndex(39)], formula)
nothing #hide #md

#=
With that, the componentwise initialization of the whole network is possible:
=#
initialize_componentwise!(nw; default_overrides=interf)
nothing #hide #md
#=
Even shorter, we can just use [`initialize_from_pf!`](@ref) to do everything from
exporting the power flow model, finding the fixpoint and initializing all components:
=#
s0 = initialize_from_pf!(nw)
