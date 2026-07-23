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
- **`Sn` and `Vn` no longer default to `S_b`/`V_b`**
  ([#261](https://github.com/JuliaEnergy/PowerDynamics.jl/pull/261)) in `SauerPaiMachine`
  and `ClassicalMachine`. They were *symbolic* defaults, which under MTK v11 would turn
  into permanent parameter bindings rather than a one-time value. Pass the machine rating
  explicitly, e.g. `SauerPaiMachine(; S_b=100, V_b=1, Sn=100, Vn=1, …)`.

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
  formulation with `1/ω_b · Dt(ψ_d,ψ_q)`.
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
