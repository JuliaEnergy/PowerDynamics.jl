#=
# [IEEE39 Bus Tutorial - Part II: Initialization](@id ieee39-part2)

The goal of this tutorial is to get a understanding of the initialization process in PowerDynamics.jl

As a prequisit, we load part I of the tutorial, which contains the network model:
=#

include(joinpath(@__DIR__, "ieee39_part1.jl"))
nw # nw object now available

#=

The initialization process in PowerDynamics.jl is a two step process: first we
solve the power flow, then we use the power flow results to initialize the
indivicual network components.

There are shortcut functions to do this, but we will go through the steps in
detail for didactic reasons.

## Powerflow
To solve the power flow, we first need to get the power flow model.
We can use the funcion [`powerflow_model`](@ref).
=#
pfnw = powerflow_model(nw)
#=
The power flow model is a `Network` object like the original network.
It is build from the original network, by calling `powerflow_model` on the
individual components. For example, we have this rather complex dynamic model at bus 30:
=#
nw[VIndex(30)]
#=
We allready see from the printout, that we have a powerflow model :pvbus attached to the generator model.

We can extract the attached PV powerflow model by calling `powerflow_model`
=#
powerflow_model(nw[VIndex(30)])
#=
The function `powerflow_model` checks, if there is a powerflow model attached (it checks the `:pfmodel` [metadata](@extref NetworkDynamcis Metadata)).

Per component, the `powerflow_model` function will do the following:
1. If the model has the `:pfmodel` metadata set (see [Metadata](@extref)), it will return the `VertexModel` stored in the metadata. In [Part I](@ref ieee39-part1) we set the `:pfmodel` metadata using the `pf` keyword to the [`Bus`].
2. If the model **does not have** the `:pfmodel` metadata set, PowerDynamics will check, if the model itself is a valid powerflow model. If so, it'll just use the dynamic model as the powerflow model.

**What is a valid powerflow model?**:
A valid powerflow model is a model that has no internal dynamics,
i.e. it either has `dim(model) == 0` OR it has a zero massmatrix (i.e. only constraints).

For example, the PiLine models are completely static:
=#
nw[EIndex(1)]
#=
No internal states, not :pfmodel
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
\dot{x} &= 0 &= f_{\mathrm{nw}}(x, p, t)\\
```
Where `x` are the network states (mainly voltages $u_r$ and $u_i$ at the busses) and `p` are all the parameter such as $P$, $V$ and $Q$ values for bus models and line parameters for branch models.

To solve this rootfinding problem, we first need to find a initial guess.
The initial guess is prefilled with all the default values from the components.
In our case, all parameters and states have default values attached, so the
initial state is fully determined.
=#
pfs0 = NWState(pfnw)
#=
With the default state, we can call [`find_fixpoint`](@extref NetworkDynamics), which keeps $p$ constant and tries to find a $x$ such that
the rootfind problem state above is fulfilled:
=#
pfs = find_fixpoint(pfnw, pfs0)
#=
As a resul, we get a `NWState` object which contains the **full state** for the power flow model.

Since powerflow model and dynamic model share the same topology and the same interface (i.e. nodes create voltages, edges create currents), we can extract the interface values from the powerflow model state apply them to the dynamic model:
=#
interf = interface_values(pfs)
#=
The inteface values give us the **inputs** and **outputs** for every component in the network, i.e. for all buses we get values for
- `busbar₊i_r` and `busbar₊i_i`: current input
- `busbar₊u_r` and `busbar₊u_i`: voltage ioutput

For all branches we get:
- `src₊u_r` and `src₊u_i`: source side voltage input
- `dst₊u_r` and `dst₊u_i`: destination side voltage input
- `src₊i_r` and `src₊i_i`: source side current output
- `dst₊i_r` and `dst₊i_i`: destination side current output

With those interface values fixed, we can go over to the second step of the initialization process: the initialization of the dynamic component.

## Initialization of Dynamic Components

Initialization of a bus model means, we want to find a **steady state** of the dynamics given the interface values from the power flow model.

Therefor, the equations for the bus models (recall from [Modeling Concepts](@ref)) become
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
    \color{red}{\begin{bmatrix}i_r^\mathrm{src}\\i_r^\mathrm{dst}\end{bmatrix}} &= g^\mathrm{src}_{\mathrm e}\left(x_{\mathrm e},\color{red}{ \begin{bmatrix} u_r^\mathrm{src}\\u_i^\mathrm{src}\end{bmatrix}}, \color{red}{\begin{bmatrix} u_r^\mathrm{dst}\\u_i^\mathrm{dst}\end{bmatrix}}, p_\mathrm{e}, t\right)\\
    \color{red}{\begin{bmatrix}i_r^\mathrm{dst}\\i_r^\mathrm{dst}\end{bmatrix}} &= g^\mathrm{dst}_{\mathrm e}\left(x_{\mathrm e}, \color{red}{\begin{bmatrix} u_r^\mathrm{src}\\u_i^\mathrm{src}\end{bmatrix}}, \color{red}{\begin{bmatrix} u_r^\mathrm{dst}\\u_i^\mathrm{dst}\end{bmatrix}}, p_\mathrm{e}, t\right)\\
    \end{aligned}
    ```

Notably, the unknowns can come from either the set of states or the set of parameters.
We devide the set of symbols into two sets:
- **fixed** symbols have a `default` metadata set. They are considered fixed in the solution of the nonlinear system.
- **free** symbols only have a `guess` metadata set. They are considere free in the solution of the nonlinear system.

Lets take a look at the bus model at bus 30:
=#
gen = nw[VIndex(30)]
#=
We have a system with 15 States and 42 parameters. Of those, most are fixed, i.e. have a `default` metadata set.
In the VertexModel printout, defauls are shown with `=` while guesses are shown with `≈`.
We can use [`dump_initial_state`](@extref) to get an overview of the free and set states:
=#
dump_initial_state(gen; obs=false)
#=
Right now, we have 2 free parameters, 2 free inputs, 2 free outputs and 15 free states.
We can use [`initialize_component`](@extref NetworkDynamics) to find values for the free symbols:
=#
try #hide
initialize_component(gen)
catch e #hide
    @error e.msg #hide
end #hide
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
nothing #hide
#=
!!! note "Mutating vs non-mutating Initialization"
    The non mutating [`initialize_component`](@extref) returns a dictionary
    containing a full initialized component state. That can be useful for certain purposes, often it is easier to work with the *mutating* version of initialization function, which will write the initialized values back to the component metadata (i.e setting the `:init` property for the symbols).

Lets all the mutating `initialize_component!` function to write the initialized values back to the component and inspect the initial state using `dump_initial_state`:
=#
initialize_component!(gen; default_overrides=interf_v30)
dump_initial_state(gen; obs=false)

#=
In the state dump we see how the initialization sucessfully set all the previously unknown values, including the control parameters `avr₊vref` and `gov₊p_ref`. We have our first initialized component!

However, in practice its not allways so easy.

## Handling Structurally Underconstraint Components
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
try #hide
initialize_component!(nw[VIndex(39)]; default_overrides=interf_v39)
catch e #hide
    @error e.msg #hide
end #hide

#=
Even thugh the we set the interface values, the problem is still underconstraint!
Lets check the free symbols:
=#
println("free u: ", free_u(nw[VIndex(39)]))
println("free p: ", free_p(nw[VIndex(39)]))
nothing #hide
#=
We see 8 free states and 3 free parameters, however we only have 8 state + 2 output equations:
=#
nw[VIndex(39)] #hide

#=
Even though we have enough set parameters to initialize machine and load on its own, we cannot do it simultanously.
Intuitivly speaking, its just not clear for the solver which of the two components provides how much of power.

To solve this, we have essentially 3 methods:

### Method 1: Manual setting of defaults
The simples solution is to manually set more defaults. For example, we know that we want to initialize the ZIP load around the initialization point, i.e. $V_\mathrm{set}$ should be the same as the bus voltage magnitude.
=#
vm_manual = copy(nw[VIndex(39)])
u_r = get_initial_state(vm_manual, :busbar₊u_r)
u_i = get_initial_state(vm_manual, :busbar₊u_i)
set_default!(vm_manual, :ZIPLoad₊Vset, sqrt(u_r^2 + u_i^2))
initialize_component!(vm_manual)
nothing # hide
#=
The initialization suceeded now!

!!! note "No more default_overrides"
    Note how we can skip the `default_overrides` keyword argument, since the first (failing) call of `initialize_component!` allready "burned in" the default overrides! Mutating state is a powerful tool, but needs care!

### Method 2: Adding an `init_formula`
The problem with the previous method is, that it is quite manual.
In reality, we would never go thorug this very manual initialization process.
The fact, that the model is structurally underconstraint is a property of the model and should therefore be handled by the mode. To do so, NetworkDynamics.jl proveds the [`InitFormula`](@ref) mechanism.

An `InitFormula` is a symbolic formula that is evaluated during the initialization process. It is attached to the VertexModel so it can be
evaluated automaticially during the initialization process.
The "formula" we want to apply is simply the equation
```math
V_\mathrm{set} = \sqrt{u_r^2 + u_i^2}
```
=#
vm_formula = copy(nw[VIndex(39)])
formula = @initformula :ZIPLoad₊Vset = sqrt(:busbar₊u_r^2 + :busbar₊u_i^2)
set_initformula!(vm_formula, formula)
vm_formula # hide
#=
The printout shows 1 additional initialization equation was attached to the model.

The initialization works now:
=#
initialize_component!(vm_formula)
nothing # hide
#=
The init formula is applied early in the initialization process, essentially writing a new `default` for `ZIPLoad₊Vset` based on the other defaults.

This reduced the number of free variables to 10, thus the system was solvable.

### Method 3: Using a `InitConstraint`
Sometimes, you additional initialization needs are more complicated.
Similar to defining a *formula*, which is evaluated *before* the actual initialization, NetworkDynamics provids a mechanis for injecting additional
constraints into the initialization process.

In contrast to the formula, the constarint does not need to be explicitly sovlabe, as it defines a residual equation
```math
0 = c(x) = V_\mathrm{set} - \sqrt{u_r^2 + u_i^2}
```
=#
vm_constraint = copy(nw[VIndex(39)])
constraint = @initconstraint begin
  :ZIPLoad₊Vset - sqrt(:busbar₊u_r^2 + :busbar₊u_i^2)
end
set_initconstraint!(vm_constraint, constraint)
vm_constraint # hide
#=
With this added constaint, the initialization process is solvable again, since we now have 11 equations for the 11 free variables.
=#
initialize_component!(vm_constraint)
nothing #hide
#=
For this particular case, method (2) is the way to go. However there are cases where the constraint is more complex and cannot be expressed as a formula.

See NetworkDynamics docs on [Advanced Component Initialization: Formulas and Constraints](@extref NetworkDynamics) and the PowerDynamics specific extension
[Advanced Component Initialization](@ref) for more information on method 2 and 3.

## Initialize all Components
Lets return from our excurse into individual component initialization and focus on the whole network again.
As we've just seen, we have structurally underconstraint components in the network.
Lets define the init formulas for the two busses which have loads and machines:
=#
formula = @initformula :ZIPLoad₊Vset = sqrt(:busbar₊u_r^2 + :busbar₊u_i^2)
set_initformula!(nw[VIndex(31)], formula)
set_initformula!(nw[VIndex(39)], formula)
nothing #hide

#=
With that, the componentwise initialization of the whole network is possible:
=#
initialize_componentwise!(nw; default_overrides=interf)
nothing #hide
#=
Even shorte we can just use [`initialize_from_pf!`](@ref) to do everything from
exporting the power flow model, finding the fixpoint and initializing all components:
=#
s0 = initialize_from_pf!(nw)
