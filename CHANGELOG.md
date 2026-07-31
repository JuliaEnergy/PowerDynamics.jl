# PowerDynamics.jl Changelog
## Version 5.0.0 Changelog

Two themes: the ecosystem move that comes down from
[NetworkDynamics v1.0/v1.1](https://github.com/JuliaDynamics/NetworkDynamics.jl/blob/main/NEWS.md)
(ModelingToolkit v11, the reworked initialization pipeline, SciML v3), and a per-unit
overhaul that turns base values from per-device constructor arguments into a property of
the bus and line.

The library saw a broad naming and unit cleanup in the process. Renamed parameters fail
loudly at construction, so the practical upgrade path is to build your models and follow
the errors — **check the docstrings for anything that no longer takes the argument you
pass it**. The parts that change *numbers* rather than names are called out explicitly
below.

### Breaking: ecosystem

- **All NetworkDynamics v1.0 breaking changes apply unchanged.** PowerDynamics
  `@reexport`s NetworkDynamics, so its API *is* part of the PowerDynamics API — most
  notably `set_mtk_defaults!` → non-mutating `set_mtk_defaults`, symbolic expressions as a
  `guess` now error in favour of `initf`/`guessf`, `find_fixpoint` takes an `NWState`, and
  MTK models are no longer simplified by `mtkcompile`. See the
  [NetworkDynamics release notes](https://github.com/JuliaDynamics/NetworkDynamics.jl/blob/main/NEWS.md)
  for the details and the rewrites. NetworkDynamics ≥ 1.1.0 is required.
- **`ModelingToolkit` → `ModelingToolkitBase` + `SciCompDSL`**
  ([#261](https://github.com/JuliaEnergy/PowerDynamics.jl/pull/261)). MTK v11 split into
  `ModelingToolkitBase` (MIT) and `ModelingToolkit` (AGPL); PowerDynamics depends only on
  the MIT half, and `@mtkmodel` moved out into `SciCompDSL`. Along with it: `Symbolics` ≥7,
  `SciMLBase` ≥3 and the `OrdinaryDiffEq*` subpackages ≥2.
- The `Bus(...)`/`Line(...)` deprecations are removed; use `compile_bus`/`compile_line`.

### Breaking: per-unit system

Base values are now part of the model structure rather than something you pass to each
device. There is a new documentation chapter, **"Per-Unit Systems"**, that explains the
whole scheme; the short version:

- **Global bases with accessors**: a system power base (100 MVA), frequency base (2π·50
  rad/s) and *fallback* voltage base, set through the new exported `set_Sbase!`,
  `set_ωbase!`/`set_fbase!` and `set_Vbase!` (plus `get_*` counterparts). They are read at
  **model construction**, so set them before building anything.
- **New `SystemBase` component** carries the global quantities (`Sbase`, `ωbase`, `ωframe`).
  Exactly one instance named `systembase` sits at the top level of every bus and line;
  `MTKBus`/`MTKLine`/`compile_bus`/`compile_line` place it for you.
- **`BusBar` and `LineEnd` carry `Vbase`**, the one genuinely local base, along with derived
  `Ibase`/`Zbase`/`Ybase` and SI observables (`u_kV`, `P_MW`, `Q_MVAr`, `i_kA`). `Sbase` is a
  three-phase power and `Vbase` line-to-line, hence `Ibase = Sbase/(√3·Vbase)` — the exported
  helper `Ibase(S, V)` changed accordingly.
- **Devices inherit their bases and no longer take `S_b`/`V_b`/`ω_b`.** Drop those arguments
  at the call site: set the power/frequency base globally and the per-bus voltage base with
  the new `MTKBus(…; Vbase=…)` / `compile_bus(…; Vbase=…)` keyword. `Vbase` propagates
  outward from there — a line end inherits from the bus it connects to (a transformer is
  just a line whose two ends resolve differently, ratio included), and a satellite/injector
  bus inherits from its hub. An explicitly set value always wins.
- **`Sn`/`Vn` are now optional**, weakly defaulting to the bus base, so a machine running on
  the system base needs no rating at all.
- **`LineEnd` requires a `side` keyword** (`:src`/`:dst`) — it selects which bus the end
  inherits from. `MTKLine` passes it automatically; hand-rolled line models must be updated.
- **New `check_base_consistency(s::NWState)`**, run automatically by `initialize_from_pf[!]`
  (`check=:error` by default, `:warn`/`:none` to downgrade), verifying that the global bases
  agree everywhere and that every line end and satellite bus matches the bus it hangs off.

### Breaking: naming and units in the library

- **Setpoints are spelled `…set`** (`Pset`, `Qset`, `Vset`, `ωset`, …); `n`/`nom` is reserved
  for genuine nameplate ratings (`Sn`, `Vn`). Models ported from elsewhere keep their
  source's names so they can still be checked against the original — the AVRs' `vref` and
  the governors' `ω_ref`/`p_ref` are unchanged.
- **EMT components name parameters for what they hold**: reactance `X` and susceptance `B`
  instead of `L` and `C`, which stored reactance/susceptance values all along. The topology
  is in the model name, so there is no ambiguity.
- **Frequency vocabulary is split three ways.** One symbol used to do three unrelated jobs:

  | symbol | kind | unit | meaning |
  |---|---|---|---|
  | `ωbase` | **unit** | rad/s | converts SI time ↔ pu; the frequency at which `X = ωbase·L`, `B = ωbase·C` were evaluated |
  | `ωframe` | **gauge** | pu | speed of the global dq frame; pinned to `1` |
  | `ωset` | **setpoint** | pu | frequency a droop/damping law is commanded to hold |

  Consequences throughout the library: an `ω0 = 2π*50` keyword on the dynamic shunts,
  branches and `ComposableInverter` models is gone (they take `ωbase` from the bus instead),
  every place where the reference frame was an invisible literal `1` now writes `ωframe`,
  and the dq cross-coupling terms of the EMT components gained the `ωframe` factor they were
  missing. Since `ωframe == 1` these are all algebraically identical to before — the OpenIPSL
  comparison tests pass unchanged.
- **Three models genuinely change behaviour** and their parameter sets do *not* carry over:
  - **`Swing`** and **`IdealDroopInverter`** had no `ωbase` at all, so with `ω` in pu the
    angle advanced ~314× too slowly. Both now use `Dt(θ) ~ ωbase*(ω - ωframe)`. `Swing`'s
    `M` is the inertia `M = 2H` in seconds (default `6.0`, i.e. `H = 3 s`) and `D` a damping
    power coefficient in pu/pu (default `2.0`); `IdealDroopInverter`'s `Kp` is a pu droop
    (default `0.05`). Retune existing parameter sets against these units. A physically sized
    machine is lightly damped, so call sites that want a quickly settling transient now pass
    a `D` well above the default.
  - **`DroopOuter`/`DroopInverter`**: `ω` is now pu rather than rad/s, making `Kp` a
    dimensionless droop (default `0.05`). Convert an existing value with
    `Kp_new = Kp_old/ωbase` — with that substitution the rewrite is an exact identity, and
    the linearization is unchanged.

  `IdealDroopInverter` and `DroopInverter` now also share the same droop defaults
  `Kp = Kq = 0.05`; they previously disagreed on `Kq`, so swapping one model for the other
  silently changed the voltage stiffness.
- Smaller renames in the same spirit — `VoltageDependentLoad`'s `Vn` → `Vset`,
  `PSSE_GENCLS`'s output port `ωout` → `SPEED_out`, and unified state/output names on the
  two PLLs. Consult the docstrings.

### Initialization

- **Backward initialization through nested components.** NetworkDynamics' new
  `initf`/`guessf` metadata lets each model carry its own init recipe; composed together
  they chain into a DAG that can fix *every* free variable from the powerflow alone, with no
  nonlinear solve at all. The new tutorial "Backward Initialization of Nested Models" builds
  a full generator bus that initializes to "No free variables!", and a new **"Advanced
  Initialization"** documentation chapter covers the mechanics.
- **`@pfinitconstraint`/`@pfinitformula` hygiene fixed**: both macros mishandled escaping, so
  a constraint could not close over local runtime variables. They now capture locals
  correctly and `show` prints the macro form you wrote.
- **`initialize_from_pf` no longer ignores `pfs0`** — the supplied start state was
  overwritten with the powerflow network, so a hand-tuned guess had no effect. `tol`/`nwtol`
  are now forwarded to the residual check.

### Library

- **Optional stator dynamics in `SauerPaiMachine`** via `SauerPaiMachine(; stator_dynamics=true)`.
- **`compile_bus`/`compile_line` accept a bare injector/branch**, wrapping it in
  `MTKBus`/`MTKLine` automatically instead of erroring.
- **New `mtkcompile` keyword on `compile_bus`/`compile_line`**, forwarded to
  `VertexModel`/`EdgeModel`.
- Bug fixes: `ClassicalMachine`'s Park transform applied the current base ratio upside down
  (wrong currents whenever `Sn ≠ Sbase`); a disconnected `ComposableInverter` kept driving
  its LCL filter and wound up to ≈70 pu internally; limited-integrator callbacks now resolve
  the integrator state through the observed-alias chain instead of assuming a name survives
  compilation (visible in `PSSE_HYGOV`); `TurbineGovTypeI` could not be constructed at all
  (it is constructible now, but remains untested — there is no reference trajectory for it).
- Saturation functions (`QUAD_SE`, `EXP_SE`, `PSSE_ESST4B`'s `FEX_function`) are AD-safe and
  use `NaNMath`, so a slightly out-of-range argument no longer throws a `DomainError`
  mid-initialization.

### Documentation

- **New chapters**: "Per-Unit Systems" and "Advanced Initialization" (see above).
- **New "Getting Started with Julia" section**
  ([#276](https://github.com/JuliaEnergy/PowerDynamics.jl/pull/276)): the setup page was
  rewritten and joined by "Environment Management" and an opinionated "How to Structure
  Research Projects" guide, aimed at users who arrive at PowerDynamics before they arrive at
  Julia.
- **New tutorial** "Backward Initialization of Nested Models"; the EMT toymodel example was
  overhauled, and the initialization and modeling-concepts chapters were reworked for v5.

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
