# Per-Unit Systems

Everything a PowerDynamics model computes is in per unit: voltages, currents, powers,
impedances and time constants. The bases exist for the two moments where that is not enough —
entering data that comes in SI units, and reading results back in SI units — and for attaching
devices that bring their own per-unit system.

There are three primitive bases. Everything else is derived from them.

| base    | unit  | scope       |                                                                         |
|:--------|:------|:------------|:------------------------------------------------------------------------|
| `Sbase` | MVA   | global      | one power base for the whole network                                    |
| `ωbase` | rad/s | global      | the frequency at which reactances were evaluated; converts SI time ↔ pu |
| `Vbase` | kV    | **per bus** | the voltage level the bus' `busbar` belongs to                          |

Derived and never set directly: `Ibase = Sbase/(√3·Vbase)`, `Zbase = Vbase²/Sbase`,
`Ybase = Sbase/Vbase²`, and `fbase = ωbase/2π` [Hz].

`Sbase` is a **three-phase** power and `Vbase` a **line-to-line** voltage — the convention
every data set and every power flow case you are likely to import uses. The derived bases
follow from it: `Ibase` is the line current, so `Sbase = √3·Vbase·Ibase`, and `Zbase` is the
per-phase (star-equivalent) impedance that transmission line data is given in, so
`Zbase = (Vbase/√3)/Ibase`. The `√3` appears only in `Ibase`; `Zbase` and `Ybase` happen to
take the same value whether you read the pair as three-phase/line-to-line or as
per-phase/line-to-neutral.

## Where the bases live

The bases live in the model, not in a global configuration. Every bus and every line carries a
[`SystemBase`](@ref) component named `systembase`, holding `Sbase`, `ωbase` and the frame speed
`ωframe`. It is added by `MTKBus`/`MTKLine`. Note that it is *duplicated* per bus and
per line rather than shared — which is why the values have to be checked for agreement, see
[Consistency](@ref) below.

`Vbase` is not on it. It belongs to the bus' [`BusBar`](@ref), and each of a line's two
[`LineEnd`](@ref)s carries its own.

```asciiart
        bus 1                       line 1→2                      bus 2
┌──────────────────────┐  ┌───────────────────────────┐  ┌──────────────────────┐
│ SystemBase           │  │ SystemBase                │  │ SystemBase           │
│  Sbase ωbase ωframe  │  │  Sbase ωbase ωframe       │  │  Sbase ωbase ωframe  │
│  (reads global)      │  │  (reads global)           │  │  (reads global)      │
├──────────────────────┤  ├─────────────┬─────────────┤  ├──────────────────────┤
│ BusBar               │  │ LineEnd src │ LineEnd dst │  │ BusBar               │
│  Vbase ──────────────────→ Vbase      ╵  Vbase ←───────── Vbase               │
│  (manual)            │  │       (inherited)         │  │  (manual)            │
├──────────────────────┤  └───────────────────────────┘  ├──────────────────────┤
│ Machine              │                                 │ Load                 │
│  Sn  Vn (manual)     │                                 │  (no local base)     │
└──────────────────────┘                                 └──────────────────────┘
```

Those are the names a model refers to: `:systembase₊Sbase`, `:busbar₊Vbase`, `:src₊Vbase`. The
arrows are the inheritance covered [below](@ref vbase-inheritance); `Sn` and `Vn` on the machine
are a device-local system, not a base of the bus. The "Internals" section of
[Modeling Concepts](@ref) describes each of the three components in full.

## Setting the bases

`Sbase` and `ωbase` are set through module-level globals, *before* you build anything:

```@example pu
using PowerDynamics
using PowerDynamics.Library

set_Sbase!(300)   # MVA
set_fbase!(60)    # Hz — same as set_ωbase!(2π*60)
nothing # hide
```

!!! warning "The setters affect future construction only"
    Every component copies the globals into its own defaults **when it is constructed**.
    Setting a base afterwards does not reach a model that already exists, and there is no
    network-level setter. Set the bases at the top of a script, before the first model.

    For the same reason the globals leak between models built in the same session. Every
    setter restores its default when called without an argument, so a script that changes a
    base should undo that at the end — a bare [`set_fbase!`](@ref)`()` is "back to 50 Hz".

`Vbase` is per bus, so it is a keyword rather than a global:

```julia
bus = MTKBus(machine; Vbase=16.5)         # on the symbolic model
vertex = compile_bus(bus; Vbase=16.5)     # or when compiling
set_default!(vertex, :busbar₊Vbase, 16.5) # or afterwards
```

[`set_Vbase!`](@ref) exists too, but it is only the *fallback* for buses that were given
nothing — useful when the whole network sits at one voltage level, pointless otherwise. Its
default is `1.0`, which makes every SI observable an identity of the pu value: a
self-announcing "no voltage base was set" rather than a plausible fiction.

## [Inheritance: who takes `Vbase` from whom](@id vbase-inheritance)

Only buses own a voltage base. Their neighbours inherit it:

- **line ends** take it from the bus they are attached to — per end, so a transformer is
  simply a line whose two ends resolve to different values, with the turns ratio falling out
  of the two bases;
- a **satellite bus** of a current injector takes it from its hub (see
  [Current Injector Bus](@ref current-injector-bus)).

Both use NetworkDynamics' [parameter sharing](@extref parameter-sharing) and both are *weak*:
a value you set yourself always wins over the inherited one.

## Device-local per-unit systems: `Sn` and `Vn`

A device is usually rated differently from the network it sits in, and its data comes on its
own rating: a machine's reactances and time constants are given on the machine's MVA rating,
not on the system base. So a device may open a **local per-unit system** and convert at exactly
one boundary — its terminal. The convention for spelling it is a nominal power and voltage:

```julia
@parameters begin
    # shadows of the bases that surround the device
    Sbase, [bound_to = :systembase₊Sbase]
    Vbase, [bound_to = :busbar₊Vbase]
    # the device's own ratings, weakly initialized from those shadows
    Sn, [initf_weak = Sbase]   # machine power rating [MVA]
    Vn, [initf_weak = Vbase]   # machine voltage rating [kV]
end
```

First the device declares local copies of the bus' `Vbase` and the system's `Sbase` with
`bound_to`, which essentially inherit their values from the surrounding bus. The local ratings
are then declared `initf_weak` *on those shadows*, which reads as "if nobody gives me a rating,
take the surrounding base".

Because that initialization is weak, a rating you do set wins over it. And because it is the
surrounding base it falls back to, leaving `Sn`/`Vn` out is a useful default rather than a
special case: the local and the system per-unit systems become identical, every ratio is
exactly `1`, and the device ends up on the system base.

Everything behind the terminal boundary is in the local system — not just the internal states,
but also the measurements the device publishes (e.g. `v_mag`, `P`, `Q`) and the control inputs
it accepts (e.g. `τ_m`, `vf`). The AVR and the governor wired to those ports therefore need no
base of their own; they work in the machine's per unit, which is why the control blocks in the
library carry no `Sbase` or `Vbase` at all.

!!! note "`n` is for *nominal*"
    `Sn`/`Vn` are the device's nominal ratings — what is written on its nameplate, and what a
    pu value can be *on*. They are not setpoints: a machine's nominal voltage `Vn` and its
    voltage setpoint `Vset` are different quantities and may legitimately differ.

    The two names are a convention, not a requirement. Models ported from elsewhere keep their
    source's spelling (the PSS/E machines call the same quantity `M_b`), and nothing stops a
    device from opening more than one local system — a turbine rated separately from the
    alternator is the classical example.

## Frequency: base, frame, setpoint, state

Four quantities in a power system model look like a frequency. Three of them are in pu, and in
an undisturbed system at default settings all three read `1`: the frame does not rotate, nothing
commands anything but nominal, and the machines settle there. That coincidence is what makes
them so easy to conflate — and it ends as soon as any one of them moves. The fourth is the unit
the other three are measured in:

| symbol   | unit  | kind     |                                                   |
|:---------|:------|:---------|:--------------------------------------------------|
| `ωbase`  | rad/s | unit     | converts SI time ↔ pu; global, on `systembase`    |
| `ωframe` | pu    | gauge    | speed of the global dq frame, pinned to `1`       |
| `ωset`   | pu    | setpoint | commanded frequency of a droop, a damping term, … |
| `ω`      | pu    | state    | rotor or PLL speed                                |

`ωframe` is what makes the reference frame explicit in equations that would otherwise carry a
bare literal `1`, as in `Dt(δ) ~ ωbase*(ω - ωframe)`. It is a settable parameter, but in
PowerDynamics 5.0 it should stay at `1.0`; it exists to make the frame speed explicit, so that
if we ever introduce simulation in a COI or reference-machine frame, the models are already
written for it.

Note that the frame speed is spelled `ωframe` rather than the sometimes-used `ωref`, to avoid
confusion with the `ref` of a reference value, i.e. a setpoint.

## Reading results in SI units

Every busbar and line end carries the SI observables `u_kV`, `P_MW`, `Q_MVAr` and `i_kA` — the
pu quantity times the matching base, and therefore exactly as meaningful as the bases you set:

```@example pu
slack = compile_bus(SlackAlgebraic(; name=:slack); vidx=1, Vbase=230)
gen   = compile_bus(Swing(; name=:machine); vidx=2, pf=pfPV(P=0.8, V=1.02), Vbase=18.0)
line  = compile_line(PiLine(; name=:piline, R=0.02, X=0.1); src=1, dst=2)

nw = Network([slack, gen], [line])
s0 = initialize_from_pf!(nw)

(; u_kV = s0[VIndex(2, :busbar₊u_kV)],
   P_MW = s0[VIndex(2, :busbar₊P_MW)])
```

0.8 pu on a 300 MVA base is 240 MW, and 1.02 pu on an 18 kV bus is 18.36 kV. Matching the
convention above, `u_kV` is a line-to-line voltage, `P_MW`/`Q_MVAr` are three-phase, and
`i_kA` is the line current. The two buses sit
at different levels, so the `PiLine` between them is a transformer in the sense of the previous
section — its two ends resolve to different `Vbase` values.

## Consistency

Bases that disagree do not fail loudly on their own — they quietly produce wrong numbers. So
[`check_base_consistency`](@ref) runs on every [`initialize_from_pf!`](@ref) and verifies that

1. all components agree on `Sbase`,
2. all components agree on `ωbase`,
3. all components agree on `ωframe`,
4. each line end matches the `Vbase` of the bus it is attached to,
5. each satellite bus matches the `Vbase` of its hub.

Components that carry none of these symbols — a hand-rolled pure-pu model without a
`systembase` — are skipped rather than failed; mixing based and baseless components is
legitimate. All findings are collected into a single message. Pass `check=:warn` or
`check=:none` to `initialize_from_pf!` if you have a good reason to allow a mismatch.

What a mismatch actually costs:

| situation                     | consequence                                                                 |
|:------------------------------|:----------------------------------------------------------------------------|
| no `Vbase` set anywhere       | pu physics exact, SI observables meaningless — the tool degrades to pure pu |
| `Vbase` wrong on one line end | wrong implied transformer ratio, wrong SI values at that end                |
| `Sbase` disagrees             | power silently mis-converted at every device boundary                       |
| `ωbase` disagrees             | time constants silently rescaled per component                              |
| `ωframe` disagrees            | frames rotate apart, no steady state exists                                 |

The voltage side degrades softly, the power and frequency side does not.

Finally, this page changed two globals, so it puts them back — the rule from the top applied
to the page itself:

```@example pu
set_Sbase!()
set_fbase!()
nothing # hide
```
