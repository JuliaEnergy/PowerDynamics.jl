# Advanced Component Initialization

[Powergrid Initialization](@ref) covered the standard case: the power flow prescribes the
voltage and current at each terminal, and every component has exactly as many free symbols as
that leaves equations. This page is about the cases where it does not work out that way —
where a component has too many unknowns, or needs information the interface values do not
carry.

Everything here is optional. If your models initialize, you never need any of it.

## Why a component can be underconstrained

Consider a bus that carries *two* injectors — a machine and a voltage-dependent load:

```@example adv
using PowerDynamics
using PowerDynamics.Library

function build_network()
    slackbus = compile_bus(SlackAlgebraic(; name=:slack); vidx=1)

    # both injectors go into one `MTKBus`, so they share a single busbar
    machine = Swing(; name=:machine)
    load = ZIPLoad(; name=:load, Pset=-0.3, Qset=-0.1, KpZ=0.5, KqZ=0.5, KpI=0, KqI=0)
    # the power flow sees only the *net* injection of both: 0.8 (machine) - 0.3 (load)
    genbus = compile_bus(MTKBus(machine, load); name=:genbus, vidx=2, pf=pfPQ(P=0.5, Q=0.1))

    line = compile_line(PiLine(; name=:piline, R=0.02, X=0.1, B_src=0.01, B_dst=0.01);
                        src=1, dst=2)

    Network([slackbus, genbus], [line])
end
nw = build_network()
nothing # hide
```

!!! note "Why a function?"
    Initialization is mutating: the values it settles on are written back into the component
    metadata. That is what makes `NWState(nw)` return the initialized state afterwards, but it
    also means a network that has been through initialization once is no longer in its
    original condition. When trying out several variants, as we do below, it is easiest to
    build a fresh network for each.

The bus has two states and two outputs, so four equations — but five unknowns. `free_u` lists
the states that carry no default and therefore have to be solved for:

```@example adv
free_u(nw[VIndex(2)])
```

and `free_p` does the same for the parameters:

```@example adv
free_p(nw[VIndex(2)])
```

Four of those five belong to the machine and are exactly the ones the previous page solved
for. The fifth, `load₊Vset`, is the voltage at which the load draws its nominal power. Nothing
in the power flow solution says what it should be — the interface only knows the *total*
current at the bus, not how it splits between machine and load. Initialization says so and
gives up:

```@example adv
try#hide
initialize_from_pf!(nw)
catch e#hide
showerror(stdout, e)#hide
end#hide
nothing # hide
```

The line to look for is the warning: *"Initialization problem is underconstrained (5 vars for
4 equations)"*, followed by the list of the five. The error after it is only the fallout — a
nonlinear solve with more unknowns than equations has no unique answer.

This is a property of the model, not of the scenario, and the fixes below differ mainly in how
much of the fix lives in the model itself.

## Level 1: pin the value by hand

The most direct fix: decide what `Vset` should be and give it a `default`, so it stops being
an unknown. The sensible choice is the voltage the bus actually settles at, and we can look
that up — solving the power flow is a separate step, and its result is an `NWState` we can
index symbolically:

```@example adv
nw = build_network()
pfs = solve_powerflow(nw)
vmag = pfs[VIndex(2, :busbar₊u_mag)]   # voltage magnitude at bus 2
```

[`set_default!`](@extref NetworkDynamics.set_default!) writes that number into the component's
metadata, which moves `load₊Vset` from the "free" to the "fixed" side:

```@example adv
set_default!(nw[VIndex(2)], :load₊Vset, vmag)
free_p(nw[VIndex(2)])
```

Four unknowns for four equations again, so the initialization goes through. We hand the power
flow solution we already have back in as `pfs`, instead of letting `initialize_from_pf!`
compute it a second time:

```@example adv
s0 = initialize_from_pf!(nw; pfs)
(; Pm = s0[VIndex(2, :machine₊Pm)], P_load = s0[VIndex(2, :load₊P)])
```

The machine came out at 0.8 pu and the load at −0.3 pu, which together make the 0.5 pu net
injection the power flow was given.

This is a perfectly good answer for a one-off study. Its cost is that the value is now
hard-wired: change a load, a line or a setpoint anywhere in the network, and the bus voltage
moves while `Vset` stays where you put it.

## Level 2: let the component compute it

The relation "the load's reference voltage is the voltage at its bus" holds at *every*
operating point. Attaching it to the component instead of writing out its value keeps it true
when the operating point changes. NetworkDynamics has two ways to express such a rule:

An [`InitFormula`](@extref NetworkDynamics.InitFormula) **computes** the value and thereby
removes an unknown:

```@example adv
nw = build_network()
# read as: at initialization, set Vset to the magnitude of the busbar voltage
formula = @initformula :load₊Vset = sqrt(:busbar₊u_r^2 + :busbar₊u_i^2)
set_initformula!(nw[VIndex(2)], formula)

s0 = initialize_from_pf!(nw)
s0[VIndex(2, :load₊Vset)]
```

Same value as before, but nobody had to look it up. The formula runs before the solve, using
the terminal voltage the power flow just prescribed.

An [`InitConstraint`](@extref NetworkDynamics.InitConstraint) instead **adds an equation**,
leaving the symbol free but tying it down:

```@example adv
nw = build_network()
# a residual: whatever the solver picks must make this expression zero
constraint = @initconstraint :load₊Vset - sqrt(:busbar₊u_r^2 + :busbar₊u_i^2)
set_initconstraint!(nw[VIndex(2)], constraint)

s0 = initialize_from_pf!(nw)
s0[VIndex(2, :load₊Vset)]
```

Same answer, different route: the formula shrinks the problem to four unknowns before the
solve; the constraint grows it to five unknowns and five equations. Prefer a formula whenever
you can solve for the symbol explicitly — it is cheaper and cannot fail to converge. Reach for
a constraint when the relation cannot be rearranged into an assignment.

## Level 3: rules that live in the model

Levels 1 and 2 attach something to a *compiled* component from the outside. When you write
your own models, the same rules can be declared next to the variable they belong to, using the
`initf` variable option:

```julia
@parameters begin
    v_ref, [guess=1, initf = v_meas]   # zero control error in steady state
end
```

and [`set_initf`](@extref NetworkDynamics.set_initf) for symbols owned by a subcomponent. The
payoff is that such rules **chain**: each block states how to recover its own state from its
own output, and once the blocks are wired together the chain runs backwards from the terminal
through the machine, the exciter and the governor, until every unknown in a deeply nested bus
model is determined without any nonlinear solve at all.

The guess-side twins `guessf` / [`set_guessf`](@extref NetworkDynamics.set_guessf) do the same
thing more softly: they only seed the solver, so approximate expressions are fine.

The tutorial [Backward Initialization of Nested Models](@ref nested-init) builds such a model
from scratch; the mechanism itself is documented in the
[NetworkDynamics initialization guide](@extref backward-flow-init).

## Level 4: rules that need the power flow solution

All of the above can only use values of the component itself. Occasionally you need something
that is only known in the *power flow* model — a parameter of the static stand-in, or a
quantity that exists nowhere in the dynamic model. PowerDynamics extends both mechanisms for
that case:

| method | purpose | use it when |
|:--|:--|:--|
| [`PFInitFormula`](@ref) | sets values directly, without adding equations | you can compute the value from the power flow results |
| [`PFInitConstraint`](@ref) | adds constraint equations that must be satisfied | you need to enforce a relationship you cannot solve for |

These are the power-flow-aware versions of NetworkDynamics'
[InitFormulas and InitConstraints](@extref
Advanced-Component-Initialization:-Formulas-and-Constraints), and the difference between them
is the same one as in level 2: a formula **reduces the number of free variables**, a constraint
**increases the number of equations**.

Inside [`@pfinitformula`](@ref) and [`@pfinitconstraint`](@ref), `@pf(:x)` refers to the
symbol `x` of the *power flow* model of the same component — its states, parameters and
observables are all available, not just the interface values. A plain `:x` still refers to the
component's own symbol, so both worlds can be mixed in one expression. Either macro also takes
a `begin ... end` block to define several rules at once:

```julia
# a single formula: take the load's setpoint from what the pf model was told to consume
formula = @pfinitformula :dynamicload₊Pset = @pf(:pq₊Pset)

# ... or several at once
formulas = @pfinitformula begin
    :dynamicload₊Pset = @pf(:pq₊Pset)
    :dynamicload₊Qset = @pf(:pq₊Qset)
end
set_pfinitformula!(my_load, formulas)
```

```julia
# the same relation, written as a residual instead
constraint = @pfinitconstraint :dynamicload₊P - @pf(:pq₊Pset)

# blocks work here too, and `@pf` is optional per symbol
constraints = @pfinitconstraint begin
    :pibranch₊X - @pf(:pibranch₊X)   # copy a parameter over from the pf model
    :P_gen - @pf(:P_load)            # power balance against the pf solution
    :AVR₊Vset - :busbar₊u_mag        # purely component-internal, no `@pf` needed
end
set_pfinitconstraint!(my_generator, constraints)
```

Both are handled automatically by [`initialize_from_pf!`](@ref): the power flow is solved
first, then every `@pf(…)` is replaced by its value — "specializing" the rules — which turns
them into ordinary `InitFormula`s and `InitConstraint`s that are passed down to NetworkDynamics'
[component initialization](@extref Single-Component-Initialization). The process is transparent:
define the rules, then use `initialize_from_pf[!]` as usual.

If a dynamic model and its power flow stand-in share the same parameters and should keep the
same values, [`copy_pf_parameters`](@ref) builds the corresponding formula for all of them at
once.

## The full pipeline

For reference, this is what [`initialize_from_pf!`](@ref) does, written out:

```julia
pfnw = powerflow_model(nw)                   # static network
pfs0 = NWState(pfnw)                         # starting guess
pfs  = solve_powerflow(nw; pfnw, pfs0)       # power flow solution

interf = interface_values(pfs)               # u/i at every terminal
s0 = initialize_componentwise!(              # solve each component
    nw;
    default_overrides = interf,
    additional_initconstraint = specialize_pfinitconstraints(nw, pfs),
    additional_initformula = specialize_pfinitformulas(nw, pfs),
)
check_base_consistency(s0)                   # per-unit sanity check
```

Useful keyword arguments:

- **`pfnw`, `pfs0`, `pfs`** — replace any stage with your own. Passing `pfs` skips the power
  flow entirely, which is how you reuse one solution for several dynamic networks.
- **`default_overrides`** — merged *on top of* the interface values; use it to pin additional
  values for one particular initialization without touching the components.
- **`subverbose`** — print the per-component log, either for everything (`true`) or for
  selected components (`[VIndex(2), EIndex(1)]`).
- **`tol`, `nwtol`** — residual tolerances for the individual components and for the assembled
  network. Loosening them is occasionally justified for stiff models, where a physically
  negligible mismatch shows up as a large residual.
- **`check`** — severity of the per-unit base check (`:error`, `:warn`, `:none`).

Note that `additional_initconstraint` and `additional_initformula` are *not* keyword arguments
of `initialize_from_pf!` — they carry the specialized power flow rules. If you need them, call
[`initialize_componentwise!`](@extref NetworkDynamics.initialize_componentwise!) directly.

## Summary

| you want to … | use | lives in |
|:--|:--|:--|
| pin a value for one scenario | `default_overrides` | the init call |
| pin a value on a component | [`set_default!`](@extref NetworkDynamics.set_default!) | the component |
| compute a value from others | `initf` / [`InitFormula`](@extref NetworkDynamics.InitFormula) | the model |
| impose a relation you cannot solve for | [`InitConstraint`](@extref NetworkDynamics.InitConstraint) | the model |
| only help the solver start | `guessf` / [`GuessFormula`](@extref NetworkDynamics.GuessFormula) | the model |
| any of the above, from power flow values | [`PFInitFormula`](@ref) / [`PFInitConstraint`](@ref) | the model |
