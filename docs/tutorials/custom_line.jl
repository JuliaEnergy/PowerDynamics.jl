#=
# Tutorial on custom Line Models

In this tutoral we'll implement a custom line model: a line with over power protection.
To make it more interesting, this kind of line will monitor the current magnitude continuously.
The safty protection is triggered whenever the current magnitude exceeds a threshold,
however, to make the model more realistic, the protection is not triggered immediately but instead
after a delay.

This script can be downloaded as a normal Julia script [here](@__NAME__.jl). #md

To use it for some initial parturbations too, we'll also implement a basic short circuit model.

=#
using PowerDynamics
using ModelingToolkit
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

The PiLine we want to describe looks lie this.
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
                     ↓ i₁       │       i₂ ↓
                     ┴      (breaker)      ┴
Y_src = G_src+jB_src ┬          │          ┬  Y_dst = G_dst+jB_dst
                     │          │          │
                     ⏚          ⏚          ⏚
```

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

To calculate the middle voltage $V_\mathrm{m}$, we need to consider two cases:
- without shortcut, we have a simple voltage divider while
- during shortcut, the voltage is set to zero.
```math
\begin{aligned}
\bar{V}_\mathrm{m} &= \frac{V_1 \, Z_b + V_2 \, Z_a}{Z_a + Z_b} = V_1 \, (1-p_\mathrm{fault}) + V_2 \, p_\mathrm{fault}\\
V_\mathrm{m} &= \bar{V}_\mathrm{m} \cdot (1-\mathrm{shortcircuit}) + \mathrm{shortcircuit} \cdot (0 + 0 \cdot j)
\end{aligned}
```

Once we have the middle voltage definde, we can calculate the currents $i_\mathrm{a}$ and $i_\mathrm{b}$:
```math
\begin{aligned}
i_\mathrm{a} &= \frac{V_1 - V_\mathrm{m}}{Z_a}\\
i_\mathrm{b} ^= \frac{V_\mathrm{m} - V_2}{Z_b}
```
Finally, we can calculate the terminal currents using circhhoff law and the transformation ratios:
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

However, there is a "hack" to provent this issue. Lets say we want to model the compelx equation
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
## Implement the CustomPiLine MTKModel

With the equations and the knowlege on how to use complex terms within MTK
Models the definition is relatively straigh forward:
=#
@mtkmodel CustomPiLine begin
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
        ## nominal volage in the middle of the line
        V_mnormal = V₁*(1-pos) + V₂*pos #(V₁*Z_b + V₂*Z_a)/(Z_a+Z_b)
        ## if shortcircuit is active, set the middle voltage to zero
        V_m = V_mnormal * (1-shortcircuit) + shortcircuit * (0+0*im) #ifelse(shortcircuit>0, 0.0, V_mnormal)
        ## current through the two Z parts
        i_a = (V₁ - V_m) / Z_a
        i_b = (V_m - V₂) / Z_b
        ## terminal currents
        isrc = (-i_a - i₁)*r_src
        idst = (i_b - i₂)*r_dst
    end
    @equations begin
        src.i_r ~ active * simplify(real(isrc))
        src.i_i ~ active * simplify(imag(isrc))
        dst.i_r ~ active * simplify(real(idst))
        dst.i_i ~ active * simplify(imag(idst))
    end
end
nothing #hide #md
#=
Additionaly to the equations defined above, we multiply the currents by `active`.
This is equivalent of opening two ideal breakers on both ends of the line when `active=false`.

Lastly lets ensure that our model satisfies the [Branch Interface](@ref):
=#
@named piline = CustomPiLine()
isbranchmodel(piline)

#=
## Extending the model for dynamic over-current Protection

In order to implement the overcurrent proection, we need to make a plan in terms of callbacks.
Callbacks are a neat [feature of DifferentialEquations.jl](https://docs.sciml.ai/DiffEqDocs/stable/features/callback_functions/),
which allow you to stop the solver under certain conditions and trigger a user-defined affect function to change the state of the system.
Their general capability is [extended in NetworkDynamics](@extref NetworkDynamics.Callbacks).

We want to implement the following behavior:
1. Continuously monitor the current magnitude and compare to the maximal current threashhold.
2. If the maximum current is reached at time $t$, mark the line as to be switched off at time $t_\mathrm{cutoff} = t + \Delta t$.
3. Continuously monitor time of the simulation and switch off the line at $t_\mathrm{cutoff}$.

The way to implement this is by introducing 3 new parameters:
- `I_max`, the maximum current magnitude,
- `t_cutoff=Inf`, the time when the line should be switched off, *which defaults to infitiny* and
- `t_delay`, the delay time after which the line should be switched off.

Additional, we need to define 2 callbacks:
- A [ContinuousComponentCallback](@extref
  NetworkDynamics.ContinuousComponentCallback) to monitor the current
  magnitude, if the current magnitude exceeds the threshold, we set the
  `t_cutoff` to the current time plus the delay time. It also tells the solver to explicitly step to this cutoff time.
- A [DiscreteComponentCallback](@extref NetworkDynamics.DiscreteComponentCallback) which compares the current time to the cutoff time
  and switches off the line by setting the `active` parameter to `false`.

!!! note
    NetworkDynamics currently does not support Events Defined in MTK models. So we need to split the implementation:
    The new parameters need to be introduced to the `MTKModel` (extending CustomPiLine), the callbacks need to be defined
    for the *compield VertexModel*.

### Extension of the CustomPiLine MTKModel
Lets add the new parameters to the `CustomPiLine` model by *extending* the model.
Extend means, that we essentially copy-past the whole model definitions and are able to add
new parameters, equations, variables and so on.
=#
@mtkmodel ProtectedPiLine begin
    @extend CustomPiLine()
    @parameters begin
        I_max=Inf, [description="Maximum current magnitude"]
        t_cutoff=Inf, [description="Time when the line should be switched off"]
        t_delay=0.1, [description="Delay time after which the line should be switched off"]
    end
end
nothing #hide

#=
Once the model is defined, we can go through the building hierarchy outlined in [Modeling Concepts](@ref).
First, we need to form something satisfying the [MTKLine Interface](@ref)

```
 ┌───────────────────────────────────────────────┐
 │ MTKLine                                       │
 │┌─────────┐                         ┌─────────┐│
 ││ LineEnd │   ┌─────────────────┐   │ LineEnd ││
 ││  :src   ├─o─┤ ProtectedPiLine ├─o─┤  :dst   ││
 ││         │   └─────────────────┘   │         ││
 │└─────────┘                         └─────────┘│
 └───────────────────────────────────────────────┘
```
=#
branch = ProtectedPiLine(; name=:piline)
mtkline = MTKLine(branch)
nothing #hide #md
#=
Then, we take the mtkline and put it into a compiled [`EdgeModel`](@extref NetworkDynamics.EdgeModel-Tuple{}) by
calling the [`Line`](@ref) constructor
```

       ╔═════════════════════════════════════════╗
       ║ EdgeModel (compiled)                    ║
   src ║ ┌─────────────────────────────────────┐ ║ dst
vertex ║ │MTKLine                              │ ║ vertex
   u ───→│┌───────┐ ┌───────────────┐ ┌───────┐│←─── u
       ║ ││LineEnd├o┤ProtectedPiLine├o┤LineEnd││ ║
   i ←───│└───────┘ └───────────────┘ └───────┘│───→ i
       ║ └─────────────────────────────────────┘ ║
       ╚═════════════════════════════════════════╝
```
=#
protected_template = Line(mtkline; name=:protected_piline)

#=
### Definition of the Callbacks

We implement the callbacks as outlined in the NetworkDynamic docs on Callbacks.

Each callback is split in two parts: the *condition* and the *affect* function.
The condition is called at each timestep. If the condition is met, the affect is calld.
We start with implementing the overcurrent condition.
For a [`ComponentCondition`](@extref NetworkDynamics.ComponentCondition), we need to specify
which symbols we want to monitor. We need two: the parameter `I_max` and the current magnitudes at both ends.
By checking
=#
obssym(protected_template)
#=
We see, that we have lots of "observed" symbols available, for example the
current maginutes `:src₊i_mag` and `:dst₊i_mag`.
In our condition, we just need to compare the bigger of the two too the `I_max` parameter.

Note that the callback is a rootfind process, so we need to return the
difference between limit and current magnitude (which its zero once the limit is
reached).
=#
continuous_overcurrent_condition = ComponentCondition([:src₊i_mag, :dst₊i_mag], [:piline₊I_max]) do u, p, t
    p[:piline₊I_max] - max(u[:src₊i_mag], u[:dst₊i_mag])
end
nothing #hide #md
#=
Next up, we need to implement the affect function. Here we check, if the cutoff is allready scheduled.
If not, it sets the `t_cutoff` parameter to the current time plus the delay time, and tells the integrator
to explicitly step to this cutoff time (see docs on the [Integrator Interface](https://docs.sciml.ai/DiffEqDocs/stable/basics/integrator/)).
=#
overcurrent_affect = ComponentAffect([], [:piline₊t_cutoff, :piline₊t_delay]) do u, p, ctx
    if p[:piline₊t_cutoff] == Inf # otherwise, it is allready sheduled for cutoff
        tcutoff = ctx.t + p[:piline₊t_delay]
        println("Line $(ctx.src)→$(ctx.dst) overcurrent at t=$(ctx.t), sheduling cutoff at t=$tcutoff")
        p[:piline₊t_cutoff] = tcutoff
        ## tell the integrator to explicitly step to the cutoff time
        add_tstop!(ctx.integrator, tcutoff)
    end
end
nothing #hide #md
#=
We then build the overall callback by combining the condition and the affect function.
=#
continuous_overcurrent_callback = ContinuousComponentCallback(continuous_overcurrent_condition, overcurrent_affect)

#=
NEEDS TEXT
=#
discrete_overcurrent_condition = ComponentCondition([:src₊i_mag, :dst₊i_mag], [:piline₊I_max]) do u, p, t
    max(u[:src₊i_mag], u[:dst₊i_mag]) ≥ p[:piline₊I_max]
end
discrete_overcurrent_callback = DiscreteComponentCallback(discrete_overcurrent_condition, overcurrent_affect)

#=
Nextup, we need to define the cutoff callback. Since we expect the solver to explicitly hit the
cutoff time, we don't need implement a continuouse callback but instead use a discrete callback.
=#
cutoff_condition = ComponentCondition([], [:piline₊t_cutoff]) do u, p, t
    t == p[:piline₊t_cutoff]
end
cutoff_affect = ComponentAffect([], [:piline₊active]) do u, p, ctx
    println("Line $(ctx.src)→$(ctx.dst) cutoff at t=$(ctx.t)")
    p[:piline₊active] = 0 # switch off the line
end
cutoff_callback = DiscreteComponentCallback(cutoff_condition, cutoff_affect)
#=
Finally, we can add the callbacks to the protected template.
=#
set_callback!(protected_template,
    (continuous_overcurrent_callback, discrete_overcurrent_callback, cutoff_callback))
protected_template #hide #md

#=
## Simulate the IEEE39 Grid with the ProtectedLine

In the last part of this tutorial, we want to see our protected line in action.
The third part of the [IEEE39 Grid Tutorial](@ref ieee39-part3) simulates a short circuit on a line.
To do so, it uses two callbacks: one to enable the short circuit and one to disable the line. We can do this
much more elegantly now by just using the `ProtectedPiLine` model.

Lets load the first part of that tutorial to get the IEEE39 Grid model.
Also, we initialize the model (the quintessence of [part II](@ref ieee39-part2)).
=#
# EXAMPLEDIR = joinpath(pkgdir(PowerDynamics), "docs", "examples")
# include(joinpath(EXAMPLEDIR, "ieee39_part1.jl"))
# formula = @initformula :ZIPLoad₊Vset = sqrt(:busbar₊u_r^2 + :busbar₊u_i^2)
# set_initformula!(nw[VIndex(31)], formula)
# set_initformula!(nw[VIndex(39)], formula)
# s0 = initialize_from_pf!(nw; verbose=false)
nothing #hide #md

#=
### Derive Network with new line models
We need to build our own model by replacing the line model with our `ProtectedPiLine`.
For that, we need to create a small helper function, which taks the edge model
from the olde network, and creates a protected line model with similar
parameters.
This is relativly straigh forward, as we use the exact same variable names. So essentially
we just need to copy all the `:default` metadata from the edge models to copy their parameters.
=#
function protected_line_from_line(e::EdgeModel)
    new = copy(protected_template)
    ## copy src and destination information
    src_dst = get_graphelement(e)
    set_graphelement!(new, src_dst)
    ## copy parameter values
    for p in setdiff(psym(new), [:piline₊I_max, :piline₊t_cutoff, :piline₊t_delay])
        set_default!(new, p, get_default(e, p))
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

### Simualate with the line models
There are two things left to do: firstoff, we need to set the `I_max` parameter for the lines we want to protec.

We set the threashold to 120% of the powerflow solution:
=#
AFFECTED_LINE = 11
i_at_steadys = max(s0[EIndex(AFFECTED_LINE, :src₊i_mag)], s0[EIndex(AFFECTED_LINE, :dst₊i_mag)])
s0_protected[EIndex(AFFECTED_LINE, :piline₊I_max)] = 1.1 * i_at_steadys

# for i in 1:46
#     local i_at_steadys = max(s0[EIndex(i, :src₊i_mag)], s0[EIndex(i, :dst₊i_mag)])
#     s0_protected[EIndex(i, :piline₊I_max)] = 2.5 * i_at_steadys
# end

#=
Secondly, we neet to introduce some perturbation, for that we'll recreate the short circuit scenario from the IEEE39 example.
Notably, this time we only need to start the short circuit, as the protection is now "baked into" the line model.
=#
_enable_short = ComponentAffect([], [:piline₊shortcircuit]) do u, p, ctx
    @info "Short circuit activated on line $(ctx.src)→$(ctx.dst) at t = $(ctx.t)s"
    # p[:piline₊active] = 0
    p[:piline₊shortcircuit] = 1
end
shortcircuit_cb = PresetTimeComponentCallback(0.1, _enable_short)
known_cbs = filter(!Base.Fix2(isa, PresetTimeComponentCallback), get_callbacks(nw_protected[EIndex(AFFECTED_LINE)]))
set_callback!(nw_protected, EIndex(AFFECTED_LINE), (shortcircuit_cb, known_cbs...))
nw_protected[EIndex(AFFECTED_LINE)] #hide #md

prob = ODEProblem(
    nw_protected,
    uflat(s0_protected),
    (0.0, 15),
    copy(pflat(s0_protected));
    callback=get_callbacks(nw_protected)
)
break
@time sol = solve(prob, Rodas5P());
~45 seconds?

let fig = Figure(; size=(800, 600))
    ax = Axis(fig[1, 1];
        title="Current Magnitudes Across All Lines",
        xlabel="Time [s]",
        ylabel="Voltage Magnitude [pu]")

    ## Full simulation time range
    ts = range(sol.t[begin], sol.t[end], length=1000)

    ## Plot voltage magnitude for all buses
    for i in 1:46
        src_current = sol(ts, idxs=EIndex(i, :src₊i_mag)).u
        dst_current = sol(ts, idxs=EIndex(i, :dst₊i_mag)).u
        current = max.(src_current, dst_current)
        current = current ./ current[begin]

        lines!(ax, ts, current; linewidth=2)
    end
    # hlines!(ax, [2.5])
    # ylims!(ax, 0.85, 1.15)
    fig
end

let fig = Figure(; size=(800, 600))
    ax = Axis(fig[1, 1];
        title="Current Magnitudes Across All Lines",
        xlabel="Time [s]",
        ylabel="Voltage Magnitude [pu]")

    ## Full simulation time range
    ts = range(sol.t[begin], sol.t[end], length=1000)

    ## Plot voltage magnitude for all buses
    i=11
    src_current = sol(ts, idxs=EIndex(i, :src₊i_mag)).u
    dst_current = sol(ts, idxs=EIndex(i, :dst₊i_mag)).u
    current = max.(src_current, dst_current)
    # current = current ./ current[begin]

    lines!(ax, ts, current; linewidth=2)
    hlines!(ax, [3.56])
    # ylims!(ax, 0.85, 1.15)
    fig
end
