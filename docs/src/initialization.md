# Powergrid Initialization

A dynamic simulation has to start somewhere, and that somewhere is almost always a **steady
state**: all derivatives zero, the grid at rest, so that whatever happens later is caused by
the disturbance you apply and not by a badly chosen initial condition.

The difficulty is that you know the desired operating point in *power system* terms — this
bus holds 1.02 pu, that generator feeds in 0.8 pu of active power — while the simulation
needs it in *model state* terms: rotor angles, integrator states, controller setpoints.
Initialization is the translation between the two, and it happens in two steps:

1. **Power flow** — solve a purely algebraic version of the grid to obtain the voltage and
   the current at every bus.
2. **Component initialization** — for each component *separately*, find the internal states
   and the free parameters such that the component sits in steady state **and** produces
   exactly the voltage and current the power flow prescribed at its terminals.

This is the [Two-Step Initialization Pattern](@extref Component-wise-Network-Initialization)
of NetworkDynamics.jl, specialized for power grids. The rest of this page walks through the
steps on a small example; for the general machinery see the
[NetworkDynamics Initialization Guide](@extref initialization-guide).

!!! tip "Hands-on version"
    [Part II of the IEEE39 Bus Example](@ref ieee39-part2) runs the same process on a full
    network and pokes at every intermediate object.

## A small example

Two buses: a stiff grid connection (the slack) and a generator modeled by the swing equation,
joined by a pi-line.

```@example init
using PowerDynamics
using PowerDynamics.Library

# bus 1: an algebraic slack, i.e. a stiff connection to the rest of the grid
slackbus = compile_bus(SlackAlgebraic(; name=:slack); vidx=1)

# bus 2: a swing-equation machine. `pf=` attaches its power flow stand-in, see step 1
machine = Swing(; name=:machine)
genbus = compile_bus(machine; name=:genbus, vidx=2, pf=pfPV(P=0.8, V=1.02))

# the line between them: R and X in pu, B the shunt susceptance of each end
line = compile_line(PiLine(; name=:piline, R=0.02, X=0.1, B_src=0.01, B_dst=0.01);
                    src=1, dst=2)

nw = Network([slackbus, genbus], [line])
```

Note what we did *not* specify: the mechanical power of the machine, its internal voltage
magnitude, and its rotor angle. Those are outcomes of the initialization, not inputs to it.

## Step 1: The power flow model

The power flow model is a second, purely algebraic network with the same topology as the
dynamic one. Each dynamic component needs a static stand-in, and there are three ways to get
one:

1. **Attach one**, either with the `pf` keyword of [`compile_bus`](@ref) (as above) or
   afterwards with [`set_pfmodel!`](@ref). PowerDynamics ships the classical bus types:

   | constructor | fixes | free |
   |:--|:--|:--|
   | [`pfSlack`](@ref) | voltage magnitude and angle | P, Q |
   | [`pfPV`](@ref) | active power and voltage magnitude | Q, angle |
   | [`pfPQ`](@ref) | active and reactive power | voltage magnitude and angle |
   | [`pfShunt`](@ref) | a static admittance | voltage |

2. **Nothing** — if the component has no dynamics of its own it is already a valid power flow
   model and is used as-is. This covers all static loads, pi-lines, transformers and the
   algebraic slack.

3. A **hand-written static model** for anything the four constructors do not cover, attached
   with [`set_pfmodel!`](@ref) like any other.

[`powerflow_model`](@ref) assembles the static network from those pieces. Per component it
checks the `:pfmodel` metadata — which is what the `pf` keyword and [`set_pfmodel!`](@ref)
write — and falls back to the component itself when there is none. The result is an ordinary
[`Network`](@extref NetworkDynamics.Network) in the NetworkDynamics sense, with the same graph
as the dynamic one:

```@example init
pfnw = powerflow_model(nw)
```

!!! note "Voltage sources and current sources"
    A power flow model must have the same inputs and outputs as the component it replaces. If
    your bus is a *current source* rather than the usual voltage source, its power flow model
    has to be one too — the `pf*` constructors take a `current_source=true` keyword for that.
    See [On Voltage and Current Sources](@ref vc-and-cs).

## Step 2: Solving the power flow

[`solve_powerflow`](@ref) does two things. It builds the starting guess `pfs0 = NWState(pfnw)`
from the component defaults — filling in a flat start (`u_r=1`, `u_i=0`) for busbars that have
none, and falling back to `guess` metadata where there is no default. Then it calls
NetworkDynamics' [`find_fixpoint`](@extref NetworkDynamics.find_fixpoint), which solves the
algebraic system with [NonlinearSolve.jl](https://github.com/SciML/NonlinearSolve.jl):

```@example init
pfs = solve_powerflow(nw)
```

The result is an [`NWState`](@extref NetworkDynamics.NWState) of the *power flow* network —
the same kind of object a simulation works with, so all of the usual symbolic indexing
applies. Note how little is left to solve: the slack pins its own voltage, so the entire
problem is the two voltage components of bus 2.

[`show_powerflow`](@ref) collects the four classical quantities per bus into a table:

```@example init
show_powerflow(pfs)
```

The generator bus holds the 1.02 pu we prescribed and exports 0.8 pu; the slack takes up the
remainder including the losses. Individual values are reachable by name, observables such as
the voltage magnitude included:

```@example init
pfs[VIndex(2, :busbar₊u_mag)]
```

## Step 3: Interface values

The power flow state and the dynamic network share the same *network interface*: buses take
currents and produce voltages, lines take voltages and produce currents. So the solution can
be handed over one component at a time — this is what [`interface_values`](@extref
NetworkDynamics.interface_values) extracts:

```@example init
interface_values(pfs)
```

Four numbers per bus (`busbar₊u_r`, `busbar₊u_i`, `busbar₊i_r`, `busbar₊i_i`) and eight per
line (the same, for the `src` and `dst` ends). **Nothing else crosses over.** Whatever a
component needs beyond these four numbers, it has to determine on its own — which is exactly
what the next step does.

## [Step 4: Initializing the components](@id component-init-core)

Take a single bus. Its model is an input-output system

```math
\begin{aligned}
M\,\frac{\mathrm{d}}{\mathrm{d}t}x &= f(x, i, p)\\
u &= g(x, p)
\end{aligned}
```

Now fix the terminal quantities to the power flow values and demand that nothing moves:

```math
\begin{aligned}
\color{red}{0} &= f(x, \color{red}{i}, p)\\
\color{red}{u} &= g(x, p)
\end{aligned}
```

Red is known. That leaves ``\mathrm{dim}(x) + 2`` equations, and a component is initializable
when it has exactly that many unknowns. Which symbols count as unknown is decided per symbol
by its [metadata](@extref NetworkDynamics Metadata):

- a symbol with a **`default`** is **fixed** — it keeps that value,
- a symbol with only a **`guess`** is **free** — it is solved for, starting from the guess.

This problem is NetworkDynamics' [Single Component Initialization](@extref
Single-Component-Initialization), where it is described in full.

Both states *and* parameters take part. That is the key idea of the whole process: a
controller setpoint with no default is not a missing input, it is an unknown to be
**back-computed** from the operating point. Library models are built this way — in the
printout below, `=` marks a default and `≈` a guess:

```@example init
genbus
```

Two states and two parameters carry only a guess: `machine₊ω`, `machine₊θ`, `machine₊Pm`,
`machine₊V`. Four unknowns for ``2 + 2`` equations — exactly determined.

The terminal voltages do carry defaults (the flat start `1` and `0`) and the terminal currents
carry none — either way the interface values from step 3 take over. They enter as overrides
and show up as `Additional defaults` in the log below.

[`initialize_from_pf!`](@ref) does steps 1–3 internally and then solves this problem for
every component. With `subverbose=true` it reports what it did:

```@example init
s0 = initialize_from_pf!(nw; subverbose=true)
nothing # hide
```

Read the log per component:

- **`VIndex(1)`** (the slack) — *"No free variables!"*. Everything was already fixed by
  defaults, so there was nothing to solve; the residual only confirms that the power flow
  values are consistent with the model.
- **`VIndex(2)`** (the machine) — *"fully constrained"*: four unknowns, four equations, solved
  to machine precision.
- **`EIndex(1)`** (the line) — no free variables either, but note the `InitFormulas set`
  entry: this is where the line ends pick up their voltage bases from the buses they are
  attached to.

And these are the values that came out — the same symbolic indexing as before, now on the
initialized state of the *dynamic* network:

```@example init
(; Pm = s0[VIndex(2, :machine₊Pm)],
   V  = s0[VIndex(2, :machine₊V)],
   θ  = s0[VIndex(2, :machine₊θ)])
```

The mechanical power had to match the 0.8 pu the power flow assigned to the bus, the internal
voltage magnitude had to match the 1.02 pu terminal voltage, and the rotor angle came out of
the bus voltage angle. None of the three was ever specified by hand.

!!! note "Component by component, not all at once"
    Every component is initialized on its own, which keeps the nonlinear problems small and
    makes failures local: an error message names the one component that could not be solved.
    Afterwards the residual of the *complete* network is checked, which is what the closing
    `Initialized network with residual …` line reports.

## Step 5: Per-unit base consistency

Finally [`check_base_consistency`](@ref) verifies that the per-unit bases actually fit
together: one `Sbase`/`ωbase`/`ωframe` for the whole network, and a `Vbase` on each line end
matching the bus it is attached to. Mismatched bases do not fail loudly on their own — they
just silently produce wrong numbers — which is why this runs on every initialization. Use
`check=:warn` or `check=:none` if you have a good reason to allow a mismatch.

## The steps, and what wraps them

Written out, with the braces showing which wrapper covers which steps:

```julia
# extract the powerflow model                # ⎫                 ⎫
pfnw = powerflow_model(nw)                   # │                 │
# starting guess for the powerflow           # ⎬ solve_powerflow │
pfs0 = NWState(pfnw)                         # │                 │
# solve the algebraic network                # │                 │
pfs = find_fixpoint(pfnw, pfs0)              # ⎭                 ⎬ initialize_from_pf[!]
# extract the interface (u/i) values         #                   │
interf = interface_values(pfs)               #                   │
# initialize every component around them     #                   │
s0 = initialize_componentwise[!](            #                   │
    nw; default_overrides = interf           #                   │
)                                            #                   │
# verify the per-unit bases                  #                   │
check_base_consistency(s0)                   #                   ⎭
```

Every one of these is a public function, so the step-wise interface is there whenever you want
full control: swap out one stage, inspect an intermediate result, or reuse a single power flow
solution across several dynamic networks.

The sketch is slightly simplified — `initialize_from_pf[!]` also specializes
power-flow-dependent initialization rules, if a component carries any. See
[The full pipeline](@ref) for that version.

## Using the result

`s0` is an `NWState` of the dynamic network, which is all you need to start simulating —
[`uflat`](@extref NetworkDynamics.uflat) and [`pflat`](@extref NetworkDynamics.pflat) hand the
state and the parameters to the solver as plain vectors:

```@example init
using OrdinaryDiffEqRosenbrock
prob = ODEProblem(nw, uflat(s0), (0, 1.0), pflat(s0))
sol = solve(prob, Rodas5P())
sol.retcode
```

Nothing happens over that second, which is the point: a steady state stays where it is until
something disturbs it.

The mutating version (`!` at the end) writes the results back into the component metadata, so
`NWState(nw)` returns the initialized state again and later calls start from it. The non-mutating
[`initialize_from_pf`](@ref) leaves the network untouched and only returns the state — useful
when you initialize the same network for several scenarios.

## When initialization fails

The two steps fail in recognizably different ways.

**The power flow does not converge.** Then no dynamic component is even attempted. Usually
the setpoints are infeasible (too much load for the network, a voltage that cannot be held)
or the flat start is too far from the solution. You can pass a better starting point via
`pfs0`, or solve the power flow separately with [`solve_powerflow`](@ref) and hand the result
in as `pfs=`.

**A component cannot be initialized.** The log names it and the message says what went wrong:

- *"Initialization problem is underconstrained (N vars for M equations)"* — the component has
  more unknowns than the interface values determine. Typical for a bus carrying several
  injectors: nothing tells the solver how the power splits between them. This is a property of
  the model, and the fix belongs in the model — see
  [Advanced Component Initialization](@ref).
- *the residual is too large* — a solution was found but it is not a steady state. Often the
  operating point is not one the model can actually hold (a limiter is active, a setpoint is
  outside its range), or the solver converged to a different root than intended. Inspect the
  component with [`dump_initial_state`](@extref NetworkDynamics.dump_initial_state) and give
  the free symbols better guesses.
- *the bases are inconsistent* — see step 5 above.

!!! tip "Handing a failed initialization to an LLM"
    LLMs are good at finding initialization problems in user-defined models, provided they get
    everything the solver saw. Let's provoke such a problem on purpose by pinning
    `machine₊V = 1.0`, which takes the symbol out of the set of free variables:
    ```@example init
    set_default!(nw[VIndex(2)], :machine₊V, 1.0)
    try#hide
    initialize_from_pf!(nw)
    @assert false "Expected to throw"#hide
    catch e#hide
    e isa AssertionError && rethrow()#hide
    showerror(stdout, e)#hide
    end#hide
    ```
    Initialization fails, and the error already reports the residual of every model equation.
    To inspect the failed state instead of just the error, run it again with the tolerance
    checks disabled — the log then also shows that the problem became overconstrained, three
    unknowns for four equations:
    ```@example init
    initialize_from_pf!(nw; tol=Inf, nwtol=Inf, subverbose=[VIndex(2)])
    nothing # hide
    ```
    Now that initialization ran through, its result is stored in the metadata of
    `nw[VIndex(2)]`. Let's dump the residual together with the full equations:
    ```@example init
    comp = nw[VIndex(2)]
    init_residual(comp; verbose=true)
    println("Full Equations")
    println.(comp.metadata[:full_equations])
    println("Full Output Equations")
    println.(comp.metadata[:full_outputeqs])
    nothing #hide
    ```
    Output like this becomes unreadable for a human very quickly, but it is exactly what an LLM
    needs to spot a modeling error. The last piece is the full initial state:
    ```@example init
    dump_initial_state(comp; obs=false)
    nothing # hide
    ```
    Hand all of it over together — equations, residuals and state — and describe the operating
    point you expected. That is usually enough to identify the equation that cannot be
    satisfied.

For everything that needs more than a better guess, continue with
[Advanced Component Initialization](@ref).
