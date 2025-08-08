#=
# Tutorial on custom Line Models

In this tutorial we'll implement a custom line model:
- we start by defining a PI branch with optional fault admittance,
- we combine two pi branches to one MTKLine, to essentialy model a two-branch transmission line.

To make it more interesting, we add protection logic to the branches:
- each branch continously checks the current magnituged against a limit,
- if the current exceeds the limit, the branch is switched off after a delay time.

This script can be downloaded as a normal Julia script [here](@__NAME__.jl). #md
=#
using PowerDynamics
using ModelingToolkit
using ModelingToolkit: D_nounits as Dt, t_nounits as t
using NetworkDynamics
using OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqNonlinearSolve
using CairoMakie
using Graphs
#=
## Basic PI-Line Model

We start with defining a basic PI-Line model, which is similar to the one in `PiLine_fault.jl` as an
MTKModel. This model should fulfil the [Branch Interface](@ref), i.e. it needs to have two [`Terminal`](@ref), one called
`:src` the other called `:dst`:

```
      ┌───────────┐
(src) │           │ (dst)
  o←──┤  Branch   ├──→o
      │           │
      └───────────┘
```

The PiLine we want to describe looks like this.
We have:
- two terminals `:src` and `:dst` with their
    - voltages $V_\mathrm{src}$ and $V_\mathrm{dst}$,
    - currents $i_\mathrm{src}$ and $i_\mathrm{dst}$,
- two shunt admittances $Y_\mathrm{src}$ and $Y_\mathrm{dst}$,
- an impedance $Z$, which is split into two parts $Z_a$ and $Z_b$ by the fault position $p_\mathrm{fault}$.

```
              i_src  V₁   i_a   Vₘ   i_b   V₂  i_dst
     V_src o────←────o───Z_a─→──o───Z_b─→──o────→────o V_dst
              r_src  │          │          │   r_dst
                     ↓ i₁       ↓ i_f   i₂ ↓
                     ┴          ┴          ┴
Y_src = G_src+jB_src ┬          ┬ Y_f      ┬  Y_dst = G_dst+jB_dst
                     │          │          │
                     ⏚          ⏚          ⏚
                   (fault enabled by breaker)
```
The fault admittance $Y_f = G_f + jB_f$ can represent any fault impedance.

To model this, we we introduce the internal voltages $V_1$, $V_2$ and $V_\mathrm{m}$.
We consider the equations of the PI line in quasi-static-state. Therefore,
we can use complex variables to describe the voltages and the currents.
What we need in the end are equations for the currents at the terminals, i.e. $i_\mathrm{src}$ and $i_\mathrm{dst}$
as a function of all the parameters and the given node voltages. Lets start writing down the equations:

First, we "split" the impedance $Z$ into two parts $Z_a$ and $Z_b$:
```math
\begin{aligned}
Z_\mathrm{a} &= Z \, p_\mathrm{fault}\\
Z_\mathrm{b} &= Z \, (1-p_\mathrm{fault})
\end{aligned}
```
Next, we can define the internal voltages $V_1$ and $V_2$ in terms of the terminal
voltages and the transformation ratios:
```math
\begin{aligned}
V_1 &= r_\mathrm{src} \, V_\mathrm{src}\\
V_2 &= r_\mathrm{dst} \, V_\mathrm{dst}
\end{aligned}
```
Once we have the shunt voltages, we can directly calculate the shut currents
```math
\begin{aligned}
i_1 &= Y_\mathrm{src} \, V_1\\
i_2 &= Y_\mathrm{dst} \, V_2
\end{aligned}
```

To calculate the middle voltage $V_\mathrm{m}$, we need to consider the fault admittance $Y_f$.
The fault admittance is defined as:
```math
Y_f = G_f + jB_f
```
The effective fault admittance is controlled by the shortcircuit parameter:
```math
Y_{f,\text{eff}} = \mathrm{shortcircuit} \cdot Y_f
```
When the fault is active, we apply Kirchhoff's current law at the middle node:
$i_\mathrm{a} = i_\mathrm{b} + i_f$, which leads to the middle voltage:
```math
V_\mathrm{m} = \frac{V_1 \, (1-p_\mathrm{fault}) + V_2 \, p_\mathrm{fault}}{1 + Y_{f,\text{eff}} \, Z \, p_\mathrm{fault} \, (1-p_\mathrm{fault})}
```

Once we have the middle voltage defined, we can calculate the currents $i_\mathrm{a}$, $i_\mathrm{b}$, and $i_f$:
```math
\begin{aligned}
i_\mathrm{a} &= \frac{V_1 - V_\mathrm{m}}{Z_a}\\
i_\mathrm{b} &= \frac{V_\mathrm{m} - V_2}{Z_b}\\
i_f &= Y_{f,\text{eff}} \, V_\mathrm{m}
\end{aligned}
```
Finally, we can calculate the terminal currents using Kirchhoff law and the transformation ratios:
```math
\begin{aligned}
i_\mathrm{src} &= (-i_\mathrm{a} - i_1) \, r_\mathrm{src}\\
i_\mathrm{dst} &= (i_\mathrm{b} - i_2) \, r_\mathrm{dst}
\end{aligned}
```

## Excursion: Complex variables in MTK Models
!!! warning "Complex variables are not supported in MTK Models (at least not in PowerDynamics.jl)"
    In the end, all parameters and variables of NetworkDynamic models are real-valued, therfore, we cannot
    use complex parameters or states in our MTK Models.

However, there is a "hack" to prevent this issue. Lets say we want to model the compelx equation
```math
U = Z \cdot I
```
We could expand everything in real and imaginary parts and rewrite the equations.
However, we can also use ModelingToolkits capability to have **complex terms** even without having complex variables.
=#
@variables u_r u_i i_r i_i
@parameters R, X
Ic = i_r + im * i_i
#-
Z = R + im * X
#=
Here, `Ic` and `Z` are **not** a symbolic variables, they are julia variables which points to a complex term/expression.

Using Symboics/ModelingToolkit, we can also multiply complex terms:
=#
Uc = Z * Ic
#=
By applying `real` and `imag` to the complex term, we can extract the real and imaginary parts
to form separate equations for real and imaginary part:
=#
[
    u_r ~ real(Uc),
    u_i ~ imag(Uc)
]

#=
This trick can be used inside `@mtkmodel` as well, by just defining thos complex terms in an `begin...end` block.
=#

#=
## Implement the CustomPiBranch MTKModel

With the equations and the knowlege on how to use complex terms within MTK
Models the definition is relatively straigh forward:
=#
@mtkmodel CustomPiBranch begin
    @parameters begin
        R, [description="Resistance of branch in pu"]
        X, [description="Reactance of branch in pu"]
        G_src, [description="Conductance of src shunt"]
        B_src, [description="Susceptance of src shunt"]
        G_dst, [description="Conductance of dst shunt"]
        B_dst, [description="Susceptance of dst shunt"]
        r_src=1, [description="src end transformation ratio"]
        r_dst=1, [description="src end transformation ratio"]
        ## fault parameters
        pos=0.5, [description="Fault Position (from src, percent of the line)"]
        G_f=1, [description="Fault conductance in pu"]
        B_f=0, [description="Fault susceptance in pu"]
        shortcircuit=0, [description="shortcircuit on line"]
        ## parameter to "switch of" the line
        active=1, [description="Line active or switched off"]
    end
    @components begin
        src = Terminal()
        dst = Terminal()
    end
    begin
        ## define complex variables
        Z = R + im*X
        Ysrc = G_src + im*B_src
        Ydst = G_dst + im*B_dst
        Yf = G_f + im*B_f
        Vsrc = src.u_r + im*src.u_i
        Vdst = dst.u_r + im*dst.u_i
        ## define Z_a and Z_b in terms of Z
        Z_a = Z * pos
        Z_b = Z * (1-pos)
        ## define internal voltages using the
        V₁ = r_src * Vsrc
        V₂ = r_dst * Vdst
        ## currents through the shunt admittances
        i₁ = Ysrc * V₁
        i₂ = Ydst * V₂
        ## effective fault admittance (controlled by shortcircuit)
        Yf_eff = shortcircuit * Yf
        ## middle voltage with fault admittance effect
        V_m = (V₁*(1-pos) + V₂*pos) / (1 + Yf_eff * Z * pos * (1-pos))
        ## fault current to ground
        i_f = Yf_eff * V_m
        ## current through the two Z parts
        i_a = (V₁ - V_m) / Z_a
        i_b = (V_m - V₂) / Z_b
        ## terminal currents
        isrc = (-i_a - i₁)*r_src
        idst = (i_b - i₂)*r_dst
    end
    @equations begin
        src.i_r ~ active * real(isrc)
        src.i_i ~ active * imag(isrc)
        dst.i_r ~ active * real(idst)
        dst.i_i ~ active * imag(idst)
    end
end
nothing #hide #md
#=
Additionaly to the equations defined above, we multiply the currents by `active`.
This is equivalent of opening two ideal breakers on both ends of the line when `active=false`.

Lastly lets ensure that our model satisfies the [Branch Interface](@ref):
=#
@named pibranch = CustomPiBranch()
isbranchmodel(pibranch)

#=
## Extending the model for dynamic over-current Protection

In order to implement the overcurrent protection, we need to make a plan in terms of callbacks.
Callbacks are a neat [feature of DifferentialEquations.jl](https://docs.sciml.ai/DiffEqDocs/stable/features/callback_functions/),
which allow you to stop the solver under certain conditions and trigger a user-defined affect function to change the state of the system.
Their general capability is [extended in NetworkDynamics](@extref NetworkDynamics.Callbacks).

We want to implement the following behavior:
1. Continuously monitor the current magnitude and compare to the maximal current threshold.
2. If the maximum current is reached at time $t$, mark the line as to be switched off at time $t_\mathrm{cutoff} = t + \Delta t$.
3. Continuously monitor time of the simulation and switch off the line at $t_\mathrm{cutoff}$.

The way to implement this is by introducing 3 new parameters:
- `I_max`, the maximum current magnitude,
- `t_cutoff=Inf`, the time when the line should be switched off, *which defaults to infinity* and
- `t_delay`, the delay time after which the line should be switched off.

For robust overcurrent protection, we need to implement **multiple complementary callbacks**:
- A **continuous callback** that detects smooth threshold crossings using root-finding
- A **discrete callback** that catches instantaneous jumps above the threshold
- A **cutoff callback** that switches off the line at the scheduled time

This dual detection approach is necessary because discrete events (like short circuits) can cause
the current to jump above the threshold without crossing it smoothly, which continuous
callbacks might miss. Both overcurrent callbacks share the same affect function that schedules
the line cutoff, while the cutoff callback actually switches off the line.

!!! note
    NetworkDynamics currently does not support Events defined in MTK models. So we need to split the implementation:
    The new parameters need to be introduced to the `MTKModel` (extending CustomPiBranch), the callbacks need to be defined
    for the *compiled VertexModel*.

### Extension of the CustomPiBranch MTKModel
Let's add the new parameters to the `CustomPiBranch` model by *extending* the model.
Extend means that we essentially copy-paste the whole model definitions and are able to add
new parameters, equations, variables and so on.

We add an additional "observed" state `I_mag`, which allways contains the current magnitude at the src or dst terminal (whatever is higher).
=#
@mtkmodel ProtectedPiBranch begin
    @extend CustomPiBranch()
    @parameters begin
        I_max=Inf, [description="Maximum current magnitude"]
        t_cutoff=Inf, [description="Time when the line should be switched off"]
        t_delay=0.1, [description="Delay time after which the line should be switched off"]
    end
    @variables begin
        I_mag(t), [description="Current magnitude at src or dst terminal"]
    end
    @equations begin
        I_mag ~ max(sqrt(src.i_r^2 + src.i_i^2), sqrt(dst.i_r^2 + dst.i_i^2))
    end
end
nothing #hide

#=
Once the model is defined, we can go through the building hierarchy outlined in [Modeling Concepts](@ref).
First, we need to form something satisfying the [MTKLine Interface](@ref).

## Creating the Dual-Branch MTKLine

Here we implement our dual-branch architecture by creating two separate `ProtectedPiBranch` instances and combining them into a single `MTKLine`. This creates a line model with two parallel branches:

```
 ┌───────────────────────────────────────────┐
 │MTKLine   ┌─────────────────────┐          │
 │         ┌┤ ProtectedPiBranch A ├┐         │
 │┌───────┐│└─────────────────────┘│┌───────┐│
 ││LineEnd├o                       o┤LineEnd││
 │└───────┘│┌─────────────────────┐│└───────┘│
 │  :src   └┤ ProtectedPiBranch B ├┘  :dst   │
 │          └─────────────────────┘          │
 └───────────────────────────────────────────┘
```

The end terminals of both branches are connecte to the same physical line end. However, the branches
operate independently:
- Each branch monitors its own current magnitude (`pibranchA₊I_mag`, `pibranchB₊I_mag`)
- Each has independent protection parameters (`I_max`, `t_delay`, `t_cutoff`)
- Each can be individually switched off (`pibranchA₊active`, `pibranchB₊active`)
- Electrical parameters are adjusted so that parallel combination matches the original single-branch behavior

=#
branchA = ProtectedPiBranch(; name=:pibranchA)
branchB = ProtectedPiBranch(; name=:pibranchB)
mtkline = MTKLine(branchA, branchB)
nothing #hide #md
#=
Then, we take the mtkline and put it into a compiled [`EdgeModel`](@extref NetworkDynamics.EdgeModel-Tuple{}) by
calling the [`Line`](@ref) constructor
```

       ╔═══════════════════════════════════════════════╗
       ║ EdgeModel (compiled)                          ║
       ║ ┌───────────────────────────────────────────┐ ║
   src ║ │MTKLine   ┌─────────────────────┐          │ ║ dst
vertex ║ │         ┌┤ ProtectedPiBranch A ├┐         │ ║ vertex
   u ───→│┌───────┐│└─────────────────────┘│┌───────┐│←─── u
       ║ ││LineEnd├o                       o┤LineEnd││ ║
   i ←───│└───────┘│┌─────────────────────┐│└───────┘│───→ i
       ║ │  :src   └┤ ProtectedPiBranch B ├┘  :dst   │ ║
       ║ │          └─────────────────────┘          │ ║
       ║ └───────────────────────────────────────────┘ ║
       ╚═══════════════════════════════════════════════╝
```
=#
protected_template = Line(mtkline; name=:protected_piline)

#=
### Definition of the Callbacks

We implement the callbacks as outlined in the NetworkDynamic docs on Callbacks.

For robust overcurrent protection, we need **two complementary callbacks**:
- A **continuous callback** that detects smooth threshold crossings using root-finding
- A **discrete callback** that catches instantaneous jumps above the threshold

This dual approach is necessary because discrete events (like short circuits) can cause
the current to jump above the threshold without crossing it smoothly, which continuous
callbacks might miss.

#### Overcurrent Detection Callbacks

For [`ComponentCondition`](@extref NetworkDynamics.ComponentCondition), we need to specify
which symbols to monitor. By checking the available "observed" symbols (i.e. derive symbols that are not directly part of the states or parameters):
=#
obssym(protected_template)
#=
We see current magnitudes `:src₊i_mag` and `:dst₊i_mag` are available.
We'll monitor both ends and compare the maximum to the `I_max` parameter.

**Condition Definitions:**

The **continuous condition** uses root-finding, returning the difference between
limit and current magnitude (zero when the limit is reached):
=#
function continuous_overcurrent_condition(branchname)
    I_mag = Symbol(branchname, "₊", :I_mag) # pilineX₊I_mag
    I_max = Symbol(branchname, "₊", :I_max) # pilineX₊I_max
    t_cutoff = Symbol(branchname, "₊", :t_cutoff) # pilineX₊t_cutoff
    ComponentCondition([I_mag], [I_max, t_cutoff]) do u, p, t
        ## return max if allready sheduled
        p[t_cutoff] != Inf && return Inf
        p[I_max] - u[I_mag]
    end
end
nothing #hide #md

#=
The **discrete condition** uses a boolean check that triggers whenever
the current exceeds the threshold:
=#
function discrete_overcurrent_condition(branchname)
    I_mag = Symbol(branchname, "₊", :I_mag) # pilineX₊I_mag
    I_max = Symbol(branchname, "₊", :I_max) # pilineX₊I_max
    t_cutoff = Symbol(branchname, "₊", :t_cutoff) # pilineX₊t_cutoff
    ComponentCondition([I_mag], [I_max, t_cutoff]) do u, p, t
        ## return false if allready sheduled
        p[t_cutoff] != Inf && return false
        u[I_mag] ≥ p[I_max]
    end
end
nothing #hide #md

#=
**Shared Affect Function:**

Both callbacks use the same affect function. When triggered, it schedules
the line cutoff by setting `t_cutoff` and tells the integrator to step to that time:
=#
function overcurrent_affect(branchname)
    t_cutoff = Symbol(branchname, "₊", :t_cutoff) # pilineX₊t_cutoff
    t_delay = Symbol(branchname, "₊", :t_delay)   # pilineX₊t_delay
    ComponentAffect([], [t_cutoff, t_delay]) do u, p, ctx
        p[t_cutoff] != Inf && return # return early if allready schedule for cutoff
        tcutoff = ctx.t + p[t_delay]
        println("$branchname of line $(ctx.src)→$(ctx.dst) overcurrent at t=$(ctx.t), scheduling cutoff at t=$tcutoff")
        p[t_cutoff] = tcutoff
        ## tell the integrator to explicitly step to the cutoff time
        add_tstop!(ctx.integrator, tcutoff)
    end
end
nothing #hide #md

#=
#### Line Cutoff Callback

The cutoff callback switches off the line when the scheduled cutoff time is reached.
Since we expect the solver to explicitly hit the cutoff time (via `add_tstop!`),
we only need a discrete callback:
=#
function cutoff_condition(branchname)
    t_cutoff = Symbol(branchname, "₊", :t_cutoff) # pilineX₊t_cutoff
    ComponentCondition([], [t_cutoff]) do u, p, t
        t == p[t_cutoff]
    end
end
function cutoff_affect(branchname)
    active = Symbol(branchname, "₊", :active) # pilineX₊active
    ComponentAffect([], [active]) do u, p, ctx
        println("$branchname of line $(ctx.src)→$(ctx.dst) cutoff at t=$(ctx.t)")
        p[active] = 0 # switch off the line
    end
end
function cutoff_callback(branchname)
    DiscreteComponentCallback(cutoff_condition(branchname), cutoff_affect(branchname))
end
nothing #hide #md

#=
#### Adding Callbacks to Template
We build both callbacks by combining their respective conditions and affects.
Finally, we add all three callbacks to the protected template:
=#
function branch_callbacks(branchname)
    oc_affect = overcurrent_affect(branchname)
    oc1 = ContinuousComponentCallback(
        continuous_overcurrent_condition(branchname),
        oc_affect
    )
    oc2 = DiscreteComponentCallback(
        discrete_overcurrent_condition(branchname),
        oc_affect
    )
    cut = DiscreteComponentCallback(
        cutoff_condition(branchname),
        cutoff_affect(branchname)
    )
    (oc1, oc2, cut)
end
set_callback!(protected_template, branch_callbacks(:pibranchA))
add_callback!(protected_template, branch_callbacks(:pibranchB))
protected_template #hide #md

#=
## Simulate the IEEE39 Grid with the ProtectedLine

In the last part of this tutorial, we want to see our protected line in action.
The third part of the [IEEE39 Grid Tutorial](@ref ieee39-part3) simulates a short circuit on a line.
To do so, it uses two callbacks: one to enable the short circuit and one to disable the line. We can do this
much more elegantly now by just using the `ProtectedPiBranch` model.

Lets load the first part of that tutorial to get the IEEE39 Grid model.
Also, we initialize the model (the quintessence of [part II](@ref ieee39-part2)).
=#
EXAMPLEDIR = joinpath(pkgdir(PowerDynamics), "docs", "examples")
# include(joinpath(EXAMPLEDIR, "ieee39_part1.jl"))
# formula = @initformula :ZIPLoad₊Vset = sqrt(:busbar₊u_r^2 + :busbar₊u_i^2)
# set_initformula!(nw[VIndex(31)], formula)
# set_initformula!(nw[VIndex(39)], formula)
# s0 = initialize_from_pf!(nw; verbose=false)
nothing #hide #md

#=
### Derive Network with new line models
We need to build our own network model by replacing the line models with our `ProtectedPiBranch`.
For that, we create a helper function that takes an edge model from the old network and creates
a protected line model with equivalent electrical parameters.

Our protected line model uses two parallel branches (A and B), so we need to adjust the parameters.
For two parallel branches to behave like the original single branch:
- Impedances (R, X): 2× original (parallel combination gives original)
- Shunt admittances (G, B): 0.5× original (parallel combination gives original)
- Transformation ratios (r): same as original

=#
function protected_line_from_line(e::EdgeModel)
    new = copy(protected_template)
    ## copy src and destination information
    src_dst = get_graphelement(e)
    set_graphelement!(new, src_dst)
    for branch in [:pibranchA, :pibranchB]
        ## Impedances: double them (2× original)
        set_default!(new, Symbol(branch, "₊", :R), 2 * get_default(e, :piline₊R))
        set_default!(new, Symbol(branch, "₊", :X), 2 * get_default(e, :piline₊X))

        ## Shunt admittances: halve them (0.5× original)
        set_default!(new, Symbol(branch, "₊", :G_src), 0.5 * get_default(e, :piline₊G_src))
        set_default!(new, Symbol(branch, "₊", :B_src), 0.5 * get_default(e, :piline₊B_src))
        set_default!(new, Symbol(branch, "₊", :G_dst), 0.5 * get_default(e, :piline₊G_dst))
        set_default!(new, Symbol(branch, "₊", :B_dst), 0.5 * get_default(e, :piline₊B_dst))

        ## Transformation ratios: keep same
        set_default!(new, Symbol(branch, "₊", :r_src), get_default(e, :piline₊r_src))
        set_default!(new, Symbol(branch, "₊", :r_dst), get_default(e, :piline₊r_dst))
    end
    new
end
old_edgemodels = [nw[EIndex(i)] for i in 1:ne(nw)];
vertexmodels = [nw[VIndex(i)] for i in 1:nv(nw)];
new_edgemodels = protected_line_from_line.(old_edgemodels);
nothing #hide #md

#=
We can then build a new Network with those modified edgemodels:
=#
nw_protected = Network(vertexmodels, new_edgemodels)
#=
... and initialize it!
=#
s0_protected = initialize_from_pf!(nw_protected; verbose=false)

#=
As a short sanity check, lets compare the initialized values of both networks:
we do so by extracting the [`interface_values`](@extref NetworkDynamics.interface_values)
for both solutions (a dictionary of all currents and voltages (inputs and outputs of the models)))
Then we compair their values.
=#
@assert collect(values(interface_values(s0))) ≈ collect(values(interface_values(s0_protected))) #hide #md
collect(values(interface_values(s0))) ≈ collect(values(interface_values(s0_protected)))
#=
They are identical! If we would have made an error in our line model, the steady state would be most certainly different.

### Simulate with the Protected Line Models
Now that we have our protected line models ready, we need to configure them for the simulation.
First, we set the current threshold `I_max` for overcurrent protection.

We set the threshold to 130% of the power flow solution:
=#
AFFECTED_LINE = 24

for i in 1:46
    i_at_steadys = s0_protected[EIndex(i, :pibranchA₊I_mag)]
    s0_protected[EIndex(i, :pibranchA₊I_max)] = 1.3*i_at_steadys
    i_at_steadys = s0_protected[EIndex(i, :pibranchB₊I_mag)]
    s0_protected[EIndex(i, :pibranchB₊I_max)] = 1.3*i_at_steadys
end

#=
Next, we need to introduce a perturbation to test our protection system. We'll
introduce a shortcircuit with $Y_\mathrm{fault}=1\,\mathrm{pu}$ on branch A of
line 24.
Notably, we only need to start the short circuit, as the protection is
now "baked into" the line model.
=#
_enable_short = ComponentAffect([], [:pibranchA₊shortcircuit]) do u, p, ctx
    @info "Short circuit activated on branch A of line $(ctx.src)→$(ctx.dst) at t = $(ctx.t)s"
    p[:pibranchA₊shortcircuit] = 1
end
shortcircuit_cb = PresetTimeComponentCallback(0.1, _enable_short)
known_cbs = filter(cb -> !(cb isa PresetTimeComponentCallback), get_callbacks(nw_protected[EIndex(AFFECTED_LINE)]))
set_callback!(nw_protected, EIndex(AFFECTED_LINE), (shortcircuit_cb, known_cbs...))
nw_protected[EIndex(AFFECTED_LINE)] #hide #md

#=
With all those callbacks set, we can go ahead simulating the system.
=#
prob = ODEProblem(
    nw_protected,
    uflat(s0_protected),
    (0.0, 15),
    copy(pflat(s0_protected));
    callback=get_callbacks(nw_protected),
    dtmax=0.01,
)
sol = solve(prob, Rodas5P());

#=
When we run this simulation, the console output will show that the short circuit on Branch A activates at t=0.1s
and leads to a line shutdown at t=0.2s (after the 0.1s delay), which clears the fault.
However, due to the introduced dynamics in the system, branch B of the affected line also experiences
a current magnitude exceeding the threshold, which leads to a shutdown of the line at around 1.5 seconds.

Let's look at the current magnitude evolution during the simulation:
=#

fig = let fig = Figure()
    ax = Axis(fig[1, 1];
        title="Current Magnitude Across All Lines",
        xlabel="Time [s]",
        ylabel="Current Magnitude (reltive to steady state)")

    ## Full simulation time range
    ts = range(sol.t[begin], sol.t[end], length=3000)

    ## Plot current magnitude for but the failing line
    for i in 1:46
        i == AFFECTED_LINE && continue
        current = 2*sol(ts, idxs=EIndex(i, :pibranchA₊I_mag)).u
        current = current ./ current[begin]
        lines!(ax, ts, current)
    end

    A_current = sol(ts, idxs=EIndex(AFFECTED_LINE, :pibranchA₊I_mag)).u
    B_current = sol(ts, idxs=EIndex(AFFECTED_LINE, :pibranchB₊I_mag)).u
    A_current = A_current ./ A_current[begin]
    B_current = B_current ./ B_current[begin]

    lines!(ax, ts, A_current; linewidth=2, color=:blue, label="Branch A")
    lines!(ax, ts, B_current; linewidth=2, color=:red, label="Branch B")

    hlines!(ax, [1.3]; color=:black, linestyle=:dot)
    xlims!(ax, ts[begin], ts[end])
    axislegend(ax; pos=:tr)
    fig
end

#=
We can also zoom into the time range around the short circuit to see how the current of branch B
crosses the threshold and the branch is disabled shortly after.
=#
xlims!(0, 2)
ylims!(0.8, 1.5)
fig
