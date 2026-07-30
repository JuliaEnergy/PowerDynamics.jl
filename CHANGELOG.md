# PowerDynamics.jl Changelog
## Version 5.0.0 Changelog

PowerDynamics v5 is the downstream half of
[NetworkDynamics v1.0](https://github.com/JuliaDynamics/NetworkDynamics.jl/blob/8628b8a/NEWS.md):
ModelingToolkit v11, the reworked initialization
pipeline, and the SciML v3 stack. Most of the breakage comes from there rather than from
PowerDynamics itself.

### Breaking

- **All NetworkDynamics v1.0 breaking changes apply unchanged.** PowerDynamics
  `@reexport`s NetworkDynamics, so its API *is* part of the PowerDynamics API. Affected:
  `set_mtk_defaults!` → non-mutating `set_mtk_defaults`, symbolic expressions as a `guess`
  or bound to an unknown now error (use `initf`/`guessf`), `find_fixpoint` takes an
  `NWState`, MTK models are no longer simplified by `mtkcompile`, `ODEProblem` passes
  `initializealg=BrownFullBasicInit()`, and `VectorContinuousComponentCallback` lost
  `affect_neg!`. See the
  [NetworkDynamics v1.0 release notes](https://github.com/JuliaDynamics/NetworkDynamics.jl/blob/8628b8a/NEWS.md)
  for the details and the rewrites.
- **`ModelingToolkit` → `ModelingToolkitBase` + `SciCompDSL`**
  ([#261](https://github.com/JuliaEnergy/PowerDynamics.jl/pull/261)). MTK v11 split into
  `ModelingToolkitBase` (MIT) and `ModelingToolkit` (AGPL); PowerDynamics depends only on
  the MIT half, and `@mtkmodel` moved out into `SciCompDSL`. Along with it: `Symbolics` ≥7, `SciMLBase` ≥3, `NetworkDynamics` ≥1, and the
  `OrdinaryDiffEq*` subpackages ≥2.
- **Per-unit bases are now part of the model, not per-device constructor arguments.** On
  v4 every base was an ordinary parameter you had to pass to each device
  (`SauerPaiMachine(; S_b=100, V_b=18, ω_b=2π*60, Sn=100, Vn=18)`), and `BusBar`/`LineEnd`
  carried no bases at all. In 5.0 the bases live on the bus/line and the devices inherit
  them. Requires NetworkDynamics ≥ 1.1.0. The pieces:
  - **Module-level globals with accessors.** A global power base (100 MVA), frequency base
    (2π·50 rad/s) and *fallback* voltage base (1.0 kV), reached through the new exported
    `get_Sbase`/`set_Sbase!`, `get_ωbase`/`set_ωbase!`, `get_fbase`/`set_fbase!` (Hz) and
    `get_Vbase`/`set_Vbase!`. They are read at **model construction**, so set them *before*
    building anything; changing them afterwards does not touch models that already exist.
    Each setter defaults its argument to the original value, so a script that changes a base
    can restore it with a bare `set_Sbase!()`.
  - **New `SystemBase` component** carrying the *global* quantities: `Sbase` [MVA], `ωbase`
    [rad/s], `ωframe` [pu], plus the observable `fbase` [Hz]. Exactly one instance named
    `systembase` sits at the top level of every bus and every line. `MTKBus`/`MTKLine` add
    it; `compile_bus`/`compile_line` inject one if a hand-rolled model does not bring its
    own.
    (The name is recycled: an unused Dynawo-style `SystemBase` model with `SnRef`/`fNom`/
    `ωNom`/`ωRef0Pu`/`ω0Pu` has been **removed** from `PowerDynamics.Library`.)
  - **`BusBar` gained `Vbase`** (the one genuinely per-bus base) and **`LineEnd` gained a
    per-end `Vbase`**, together with derived `Ibase`/`Zbase`/`Ybase` observables and SI
    observables `u_kV`, `P_MW`, `Q_MVAr`, `i_kA`. Both also carry bound shadows of
    `Sbase`/`ωbase` (observables `busbar₊Sbase`, `src₊Sbase`, …).
  - A compiled bus exposes `busbar₊Vbase` plus `systembase₊Sbase/ωbase/ωframe` as its only
    base parameters; a line exposes `src₊Vbase`, `dst₊Vbase` and the same three.
  - **`ωframe`** is the speed of the global dq reference frame in pu. It is **pinned to `1`**
    in 5.0 and exists so that frame-dependent terms name the frame explicitly instead of
    hiding an invisible `1`. Do not change it.
- **Devices no longer take `S_b`/`V_b`/`ω_b`.** `SauerPaiMachine`, `ClassicalMachine`, the
  PSS/E machines (`PSSE_GENCLS`, `PSSE_BaseMachine` and its extenders `PSSE_GENSALIENT`/
  `PSSE_GENROUND`) and `PSSE_Load` lost those constructor arguments. They now declare
  `Sbase`/`Vbase`/`ωbase` parameters that are structurally bound (`bound_to`) —
  `Sbase`/`ωbase`/`ωframe` to the container's `systembase`, `Vbase` to the `busbar` — and
  eliminated at `compile_bus`. **At a call site you simply drop them**; set the system
  power/frequency base globally (`set_Sbase!`, `set_fbase!`/`set_ωbase!`) and the per-bus
  voltage base with the new `MTKBus(…; Vbase=…)` / `compile_bus(…; Vbase=…)` keyword.
  Related: the PSS/E `CoB = M_b/Sbase` and `fn = ωbase/2π` are now parameter bindings
  (observables) rather than a hidden local and a free `fn=50` parameter respectively.
- **`Sn`/`Vn` are optional and weakly default to the bus base** (`SauerPaiMachine`,
  `ClassicalMachine`; the PSS/E `M_b` likewise defaults to `Sbase`). On v4 these were
  required parameters. Now the nameplate rating carries a *weak* default (`initf_weak`) to
  the enclosing bus's `Sbase`/`Vbase`: leave it unset and the machine runs on the system
  base, or set it explicitly for a genuine `Sn ≠ Sbase`. This restores the convenience that
  the symbolic-default removal
  ([#261](https://github.com/JuliaEnergy/PowerDynamics.jl/pull/261)) had to drop — now
  without a permanent parameter binding, since NetworkDynamics ≥ 1.1.0 supports weak
  defaults via `initf_weak`.
- **`Vbase` is inherited across components.** You normally set the voltage base once per bus
  and everything attached to it follows:
  - a **line end** takes `Vbase` from the bus it is connected to
    (`default_from = (:src|:dst, :busbar₊Vbase)`). A transformer is just a line whose two
    ends resolve to different values, with the ratio falling out of the two bases.
  - a **satellite/injector bus** (`compile_bus(…; current_source=true)`, connected through
    a `LoopbackConnection`) takes `Vbase` from its hub bus.

  In both cases an explicitly set value wins over the inherited one. Mechanically, the
  `VBASE[]` fallback is attached to `BusBar` as a *weak init formula* rather than a plain
  default, so that it stays distinguishable from a value you set; `compile_bus` then
  materializes it into a real default on an ordinary bus, or drops it in favour of the hub
  inheritance on a satellite.
- **`LineEnd` requires a `side` keyword.** `LineEnd(; name, side)` with `side ∈ (:src,
  :dst)` — it selects which incident bus the end inherits `Vbase` from. `MTKLine` passes it
  automatically; **hand-rolled line models must now write `LineEnd(; side=:src)` /
  `LineEnd(; side=:dst)`** instead of a bare `LineEnd()`.
- **Setpoint naming rule: `…set`, except in ported models.** A quantity a user or an outer
  controller *commands* is spelled with a `set` suffix — `Pset`, `Qset`, `Vset`, `δset`,
  and now `ωset`. The suffixes `n`/`nom` are reserved exclusively for genuine nameplate
  ratings (`Sn`, `Vn`), which may legitimately differ from the corresponding setpoint — a
  machine's nominal voltage `Vn` is not its voltage setpoint `Vset`. The exception is a
  model ported from elsewhere, which keeps its source's names so it can be checked against
  the original: hence `vref` on the AVRs, `ω_ref`/`p_ref` on the governors and `n_ref` on
  `PSSE_HYGOV` are unchanged. The rule renames a symbol only where its old name was actually
  ambiguous; an unambiguous `…_ref` on a ported model is left alone rather than churned.
- **`VoltageDependentLoad`: `Vn` → `Vset`.** It is the ZIP reference voltage — a setpoint,
  not a device rating.
- **Frame frequency `ω0` on dynamic shunts and inverters is now the bus `ωbase`.**
  `DynamicCShunt`, `DynamicParallelRCShunt` (shunts) and the `ComposableInverter` models
  (`VoltageSource`, `CurrentSource`, `SimpleGFL`, `SimpleGFLDC`, `DroopOuter`/`DroopInverter`)
  no longer take an `ω0 = 2π*50` keyword; they carry `ωbase` bound to `:systembase₊ωbase`,
  so the frame frequency comes from the bus (default still `2π*50`). Set it globally with
  `set_fbase!`/`set_ωbase!`. `DynamicSeriesRLBranch` lost its `ω0` keyword the same way.
- **Frequency vocabulary split: `ωbase` vs `ωframe` vs `ωset`.** One symbol used to do
  three unrelated jobs. They are now separated everywhere in the library:

  | symbol | kind | unit | meaning |
  |---|---|---|---|
  | `ωbase` | **unit** | rad/s | converts SI time ↔ pu; the frequency at which `X = ωbase·L`, `B = ωbase·C` were evaluated |
  | `ωframe` | **gauge** | pu | speed of the global dq frame; pinned to `1` |
  | `ωset` | **setpoint** | pu | frequency a droop/damping law is commanded to hold |

  Every place where the reference *frame* was previously an invisible literal `1` now names
  it: `Dt(δ) ~ ωbase*(ω - ωframe)` in `SauerPaiMachine`, `ClassicalMachine`, `PSSE_GENCLS`,
  `PSSE_BaseMachine`, `VariableFrequencySlack`, `Swing`, `IdealDroopInverter` and
  `DroopOuter`; and the dq cross-coupling terms of the EMT components
  (`DynamicCShunt`, `DynamicParallelRCShunt`, `DynamicSeriesRLBranch`, `LFilter`/`LCFilter`/
  `LCLFilter`) gained the explicit `ωframe` factor they were missing. Since `ωframe == 1`,
  every one of these is algebraically the same equation as before — the OpenIPSL comparison
  tests still pass at their 1e-5 tolerance. (`DynamicSeriesRLBranch` was additionally
  rewritten from `Dt(i) ~ ωbase/L*Δu - …` into the classical `L/ωbase*Dt(i) ~ Δu - R*i +
  ωframe*L*i` form, so that the only `ωbase` left is the coefficient of the derivative;
  same equation, different floating-point association.) The genuine behaviour changes are
  the three models listed below. Renames:
  - **`ω_ref` → `ωset` in `Swing`**, where the one symbol was doing both jobs at once
    (`Dt(θ) ~ ω - ω_ref` *and* `D*(ω - ω_ref)`): the frame half became `ωframe`, the
    setpoint half `ωset`. The governors (`TurbineGovTypeI`, `TGOV1`) keep `ω_ref` — there it
    only ever meant the setpoint, so there was nothing to disambiguate and no reason to
    break call sites. The IEEE39 example's `gov.csv` column stays `ω_ref` accordingly.
  - **`Swing`: the `ω_ref_input`/`ωset_input` structural parameter is removed**; `ωset` is
    always a plain parameter. It entered the model only through the damping term
    `D*(ω - ωset)`, so driving it was indistinguishable from driving `Pm` by `D*Δωset` — a
    redundant spelling of the already existing `Pm_input` port. Use `Pm_input=true` to
    actuate the machine from a controller.
  - **`PSSE_GENCLS`: output port `ωout` → `SPEED_out`**, matching `PSSE_BaseMachine`. It
    carries a speed *deviation* (OpenIPSL convention) while `ωout` on `SauerPaiMachine`/
    `ClassicalMachine` carries an absolute pu speed, so the shared name was actively
    misleading. Nothing consumed it.
  - **The two PLLs' states are renamed and unified.** See docstring and model sources for details.
    Main reason: highlight that tracking frequency is in rad/s and unify output states θ and ω.
- **`Swing` and `IdealDroopInverter` had no `ωbase` at all** — their angle equations
  (`Dt(θ) ~ ω - ω_ref`) implicitly assumed normalized time, so with `ω` in pu the rotor
  angle advanced ~314× too slowly. Both now use `Dt(θ) ~ ωbase*(ω - ωframe)`, and their
  **parameters therefore change meaning and default value**:
  - `Swing`: `M` is now the inertia `M = 2H` in **seconds** (default `0.005` → `6.0`, i.e.
    `H = 3 s`) and `D` is a damping power coefficient in pu power per pu frequency
    (default `0.0001` → `2.0`).
  - `IdealDroopInverter`: `Kp` is a droop in pu frequency per pu power (default `1` →
    `0.05`, a conventional 5 % droop) and the measurement filters `τ_p`/`τ_q` default to
    `0.1 s` instead of `1 s`.

  Existing `Swing`/`IdealDroopInverter` parameter sets do **not** carry over; retune them
  against the units above. A physically-sized machine on a stiff line is lightly damped, so
  demo/test call sites that want a quickly settling transient now pass a `D` well above the
  model default.
- **`DroopOuter`/`DroopInverter` droop is in pu.** `ω` used to be an absolute angular
  frequency in rad/s (with `ω0 = 2π*50` doubling as the frame speed), which made `Kp` a
  rad/s-per-pu-power coefficient. Now `ω` is pu, the setpoint is `ωset = 1` and
  `Dt(δ) ~ ωbase*(ω - ωframe)`, so **`Kp` is a dimensionless droop** (default `0.4` →
  `0.05`, i.e. 5 %). Convert an existing value with `Kp_new = Kp_old/ωbase` — with that
  substitution the rewrite is an exact identity (`Dt(δ)` was `-Kp_old·ΔP` and is now
  `-ωbase·Kp_new·ΔP`), and since `ω` is an algebraic observable rather than a state, the
  linearization of `DroopInverter` is unchanged.
- **`IdealDroopInverter` and `DroopOuter`/`DroopInverter` now share the same droop defaults**,
  `Kp = Kq = 0.05` (the two models previously disagreed on `Kq`, so swapping one for the other
  silently changed the voltage stiffness). `Kq = 0.05` means rated reactive output costs a 5 %
  voltage deviation; it acts as a virtual reactance in series with the coupling reactance `X`
  (`Q - Qset = (Vset - V_grid)/(X + Kq)`), and 0.05 pu is the same order as a typical converter
  transformer. Since the frame rework already forces `Kp` to be retuned, retune `Kq` with it
  rather than relying on the default.

### Fixed

- **`TurbineGovTypeI` could not be constructed at all.** Its droop equation used the
  `ω_meas` *connector* instead of `ω_meas.u` (`MethodError: no method matching -(::Num,
  ::System)`), and used the raw `p_ref` where the `p_ref_input`-aware `_p_ref` was intended,
  so a `RealInput` power reference was silently ignored. Found while auditing the frequency
  symbols; the model is now constructible, but note it remains **untested** — there is no
  reference trajectory for it in the test suite.

### Initialization

- **Backward initialization through nested components.** NetworkDynamics' new
  `initf`/`guessf` metadata and `set_initf`/`set_guessf` let each model carry its own init
  recipe; composed together they chain into a DAG that can fix *every* free variable from
  the powerflow alone, with no nonlinear solve at all. The new tutorial "Backward
  Initialization of Nested Models" builds a generator bus (machine + AVR + governor + inner
  blocks) that initializes to "No free variables!".
- **`@pfinitconstraint` / `@pfinitformula` hygiene fixed**: both macros mishandled
  escaping, so a constraint could not close over local runtime variables
  (`@pfinitconstraint :x * scale - @pf(:z)` inside a loop). They now capture locals
  correctly, and `show` prints the macro form the user wrote instead of leaking
  `Expr(:escape, …)` into the output.
- **`initialize_from_pf` no longer ignores `pfs0`** — the supplied start state was
  overwritten with the powerflow *network*, so a hand-tuned powerflow guess had no effect.
- `tol` and `nwtol` are forwarded through `initialize_from_pf` to the residual check.

### Library

- **New: optional stator dynamics in `SauerPaiMachine`** via
  `SauerPaiMachine(; stator_dynamics=true)`, which replaces the algebraic stator
  formulation with `1/ωbase · Dt(ψ_d,ψ_q)`.
- **New: `mtkcompile` keyword on `compile_bus` and `compile_line`**, forwarded to
  `VertexModel`/`EdgeModel`. Pass `true` for MTK's (AGPL) simplification pipeline or
  `:compare` to print both side by side.
- **Bug fix: `ClassicalMachine` per-unit conversion.** The Park transform applied the
  current base ratio upside down (`Ibase(Sn,Vn)/Ibase(S_b,V_b)`), giving wrong currents
  whenever the machine rating differed from the system base.
- **Bug fix: `ComposableInverter` in the disconnected state.** With `connected=0` the LCL
  filter's grid-side current was still driven by `V_C` and built up to ≈70 pu internally.
  Both `V_C` and `terminal.u` are now gated, so `i_g` decays to zero as an open circuit
  should.
- **Bug fix: limited-integrator callbacks** now resolve the integrator state through the
  observed-alias chain instead of assuming the name `<ns>₊x` survives compilation. MTK may
  demote the state to an alias of the block output, which silently broke the saturation
  callbacks (visible in `PSSE_HYGOV` and with `mtkcompile=true`).
- **Saturation functions are AD-safe**: `QUAD_SE`, `EXP_SE` and `PSSE_ESST4B`'s
  `FEX_function` are now type-generic and use `NaNMath`, so `sqrt`/`log` on a slightly
  out-of-range argument no longer throws a `DomainError` mid-initialization.
- `SlackDifferential` is built with `set_initf` rather than symbolic parameter defaults.

### Documentation

- **New "Getting Started with Julia" section**
  ([#276](https://github.com/JuliaEnergy/PowerDynamics.jl/pull/276)): the setup page was
  rewritten and joined by "Environment Management" and an opinionated "How to Structure
  Research Projects" guide, aimed at users who arrive at PowerDynamics before they arrive
  at Julia.
- **New tutorial** "Backward Initialization of Nested Models" (see above).

## Version 4.4.1 Changelog
- fixes a precompile issue on Windows machines
- small changes to documentation + test system

## Version 4.4.0 Changelog
- **New component: ComposableInverter**: Added composable grid-forming (GFM) and grid-following (GFL) inverter models with filter blocks (`LFilter`, `LCFilter`, `LCLFilter`), controller blocks (`CC1`, `VC`, `CC2`), PLL models, and composite models (`DroopInverter`, `SimpleGFL`, `SimpleGFLDC`)
- **New component: Breaker**: Added `Breaker` component for modeling switchable zero-impedance connections between buses. Includes example demonstrating breaker opening and synchronized reclosing.
- **New component: DynamicSeriesRLBranch**: Dynamic RL transmission line in rotating dq-frame with optional transformer ratios
- **New shunt models**: Added `StaticShunt`, `DynamicCShunt`, and `DynamicParallelRCShunt` for bus shunt modeling
- **New load model: ConstantCurrentLoad**: Added constant current magnitude load model with configurable phase offset relative to voltage
- **New compile_bus option**: Added `current_source` keyword argument to `compile_bus` to support injector nodes via loopback connections (switches from current-in/voltage-out to voltage-in/current-out)
- **New utility functions**: Added `unwrap_deg` and `unwrap_rad` for continuous phase trajectories
- **New saturation methods**: `LimIntegrator` and `SimpleLagLim` now support four saturation methods: `:callback` (original), `:complementary`, `:rhs_hard`, and `:rhs_soft`, configurable via `SaturationConfig`
- **New documentation**: Added linear analysis tutorial (eigenvalue analysis, Bode plots, 4-bus validation against SimplusGT), breaker tutorial, and conceptual page on voltage/current source modeling
- **Improved building blocks documentation**: Added comprehensive docstrings with ASCII art diagrams for all control system building blocks (`SimpleLag`, `SimpleLead`, `LeadLag`, `Derivative`, `SimpleGain`, `SimpleLagLim`, `LimIntegrator`, `DeadZone`)
- **New building block features**: Added `allowzeroT` structural parameter to `SimpleLag` and `LeadLag` blocks to optionally bypass lag/lead when time constants are zero
- **Exported building blocks**: All building blocks and saturation functions are now properly exported and documented in the Library API documentation
- **Enhanced limited integrator callbacks**: Refactored to use NetworkDynamics v0.10.12's `ComponentPostprocessing` metadata for cleaner callback attachment
- **Improved power flow solving**: Added sparsity support for more efficient power flow computation
- **Bug fix: EXP_SE formula**: Corrected the exponential saturation function to properly interpolate through both data points
- **Bug fix: EMT toymodel sign convention**: Fixed capacitor and RL line equations to use correct injector sign convention
- **Added saturation function tests**: New tests verify that `QUAD_SE` and `EXP_SE` correctly pass through specified points
- **Renamed saturation functions** (Library): `PSSE_QUAD_SE` → `QUAD_SE`, `PSSE_EXP_SE` → `EXP_SE`. The old names are no longer available. Users should update their code to use the new names.
- **Removed global postprocessing system**: The global `POSTPROCESSING_FUNCTIONS` array and `register_postprocessing_function!()` function have been removed. Component postprocessing is now handled via `ComponentPostprocessing` metadata in ModelingToolkit models. This change only affects users who were manually registering postprocessing functions.

## Version 4.3.0 Changelog
- [#232](https://github.com/JuliaDynamics/PowerDynamics.jl/pull/232): Added lots of new models based on the great OpenIPSL Library.
  Those models are validated against OpenIPSL by reproducing their component test-harness and comparing trajectories of internal
  states.

## Version 4.2.0 Changelog
- [#230](https://github.com/JuliaDynamics/PowerDynamics.jl/pull/230):
  - deprecate `Bus(...)` → `compile_bus(...)` and `Line(...)` → `compile_line(...)`
  - remove `PowerDynamicsTesting` as separate package and just load it as module (less env hassle)
  - add new `asciiart` code style for documentation

## Version 4.1.0 Changelog
- [#221](https://github.com/JuliaDynamics/PowerDynamics.jl/pull/221) update for ModelingToolkit.jl v10 compatibility:
  - Update minimum ModelingToolkit.jl requirement from to v10
  - Update minimum NetworkDynamics.jl requirement from v0.10.3 to v0.10.4
  - Remove internal `pin_parameters` function - was never part of public API
  - Rename all `ODESystem` → `System` throughout codebase (follows MTK v10 API)
  - Replace `structural_simplify` with `mtkcompile` for model compilation
  - Replace custom `_to_zero()` hack with `implicit_output()` from NetworkDynamics
  - Update custom metadata system to use `CustomMetadata` wrapper for safer metadata handling

## Version 4.0.0 - Major Breaking Release
In Q2 2024, we began a complete rewrite of PowerDynamics.jl, bringing it much closer in alignment with the modern SciML stack. This rewrite heavily leverages ModelingToolkit.jl for equation-based models and includes a vastly modernized version of our backend, NetworkDynamics.jl.

The 4.0.0 update incorporates all of these changes. We consider the modeling concepts and simulation tools stable enough for release. The library, however, is marked as experimental for now and may change in upcoming minor versions — but you can copy the model definitions into your own code if you rely on an “old” model.

The new version has improved significantly in terms of modeling, initialization, and solution analysis. However, some models and tools previously available are not yet available. If you want to continue using the (unmaintained) old version, stick with PowerDynamics@v3.

--------------------------
## Version 2.4.0

* ![enhancement](https://img.shields.io/badge/PD-enhancement-%23a2eeef.svg) [added Architecture.md](https://github.com/JuliaEnergy/PowerDynamics.jl/issues/52)
* ![enhancement](https://img.shields.io/badge/PD-enhancement-%23a2eeef.svg) [using github actions for CI]()

## Version 2.3.3

* ![enhancement](https://img.shields.io/badge/PD-enhancement-%23a2eeef.svg) [refactored the fault design](https://github.com/JuliaEnergy/PowerDynamics.jl/issues/87)
* ![enhancement](https://img.shields.io/badge/PD-enhancement-%23a2eeef.svg) [added OrderDicts for nodes and lines definition](https://github.com/JuliaEnergy/PowerDynamics.jl/issues/86)

## Version 2.3.0

* ![enhancement](https://img.shields.io/badge/PD-enhancement-%23a2eeef.svg) [added a voltage dependent load model](https://github.com/JuliaEnergy/PowerDynamics.jl/pull/109)
* ![enhancement](https://img.shields.io/badge/PD-enhancement-%23a2eeef.svg) [added the feature to find operationpoints using SteadyStateProblem](https://github.com/JuliaEnergy/PowerDynamics.jl/pull/97)
* ![enhancement](https://img.shields.io/badge/PD-enhancement-%23a2eeef.svg) [added new implementation of NodeShortCircuit](https://github.com/JuliaEnergy/PowerDynBase.jl/pull/93)
* ![enhancement](https://img.shields.io/badge/PD-enhancement-%23a2eeef.svg) [added dynamic RL Line](https://github.com/JuliaEnergy/PowerDynamics.jl/pull/96)
* ![enhancement](https://img.shields.io/badge/PD-enhancement-%23a2eeef.svg) [added composition of node types](https://github.com/JuliaEnergy/PowerDynamics.jl/pull/75)

## Version 0.8

* ![bugfix](https://img.shields.io/badge/PD-bugfix-%23d73a4a.svg) [fixed state conversion](https://github.com/JuliaEnergy/PowerDynBase.jl/pull/62)
* ![bugfix](https://img.shields.io/badge/PD-bugfix-%23d73a4a.svg) [Voltage measurement of the exciter in grid reference frame](https://github.com/JuliaEnergy/PowerDynBase.jl/pull/63)
* ![bugfix](https://img.shields.io/badge/PD-bugfix-%23d73a4a.svg) [removed pi-model function from Transformer.jl and added export statement](https://github.com/JuliaEnergy/PowerDynamics.jl/pull/26)

## Version 0.7

* ![enhancement](https://img.shields.io/badge/PD-enhancement-%23a2eeef.svg) [added CHANGELOG.md and check whether CHANGELOG.md has been modified to ci/travis and proper splitting of ci on travis](https://github.com/JuliaEnergy/PowerDynBase.jl/pull/36)
* ![bugfix](https://img.shields.io/badge/PD-bugfix-%23d73a4a.svg) [symbolsof is now defined on the class level via the @DynamicNode Macro](https://github.com/JuliaEnergy/PowerDynBase.jl/pull/35)
* ![bugfix](https://img.shields.io/badge/PD-bugfix-%23d73a4a.svg) & ![enhancement](https://img.shields.io/badge/PD-enhancement-%23a2eeef.svg) [add Julia 1.1. to travis/ci and fixed wrong coverage reporting (thus)](https://github.com/JuliaEnergy/PowerDynBase.jl/pull/38)
* ![enhancement](https://img.shields.io/badge/PD-enhancement-%23a2eeef.svg) [cleaning .travis.yml](https://github.com/JuliaEnergy/PowerDynBase.jl/pull/39)
* ![enhancement](https://img.shields.io/badge/PD-enhancement-%23a2eeef.svg) [fixing node docs](https://github.com/JuliaEnergy/PowerDynBase.jl/pull/46)
* ![miscellaneous](https://img.shields.io/badge/PD-miscellaneous-lightgrey.svg) [add Sabine Auer to AUTHORS](https://github.com/JuliaEnergy/PowerDynBase.jl/pull/45)
* ![bugfix](https://img.shields.io/badge/PD-bugfix-%23d73a4a.svg) [temporarily pin MbedTLS version to 0.6.6 to fix travis ci issue](https://github.com/JuliaEnergy/PowerDynBase.jl/pull/49)
* ![bugfix](https://img.shields.io/badge/PD-bugfix-%23d73a4a.svg) [adding and package update before the tests in travis to ensure the latest packges are used](https://github.com/JuliaEnergy/PowerDynBase.jl/pull/50)
* ![enhancement](https://img.shields.io/badge/PD-enhancement-%23a2eeef.svg) [add current source inverter to NodeDynamics](https://github.com/JuliaEnergy/PowerDynBase.jl/pull/52)
* ![enhancement](https://img.shields.io/badge/PD-enhancement-%23a2eeef.svg) [Add Exponential Recovery Load Model](https://github.com/JuliaEnergy/PowerDynBase.jl/pull/54)
* ![enhancement](https://img.shields.io/badge/PD-enhancement-%23a2eeef.svg) [Fourth Order Equation with Exciter, AVR and Governor](https://github.com/JuliaEnergy/PowerDynBase.jl/pull/53)
* ![enhancement](https://img.shields.io/badge/PD-enhancement-%23a2eeef.svg) [`make` has now an `open-docs` target that builds the docs and then opens it](https://github.com/JuliaEnergy/PowerDynBase.jl/pull/55)
* ![bugfix](https://img.shields.io/badge/PD-bugfix-%23d73a4a.svg) [removing accidentally added build file `.DS_Store`](https://github.com/JuliaEnergy/PowerDynBase.jl/pull/49)
* ![bugfix](https://img.shields.io/badge/PD-bugfix-%23d73a4a.svg) [adding include statement for PiModel.jl](https://github.com/JuliaEnergy/PowerDynamics.jl/pull/29)
* ![bugfix](https://img.shields.io/badge/PD-bugfix-%23d73a4a.svg) [enabling simulations without slack bus](https://github.com/JuliaEnergy/PowerDynamics.jl/pull/34)
* ![bugfix](https://img.shields.io/badge/PD-bugfix-%23d73a4a.svg) [docs fixed, inertia constant in SM model corrected](https://github.com/JuliaEnergy/PowerDynamics.jl/pull/37)
