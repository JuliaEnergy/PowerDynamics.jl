# [Part III: Structuring Research Projects in Julia](@id research-projects)

This is Part III of our Julia introduction series. It assumes you already have a working setup ([Part I](@ref julia-setup)) and understand environments and `Pkg` ([Part II](@ref env-management)).

Once a research project grows beyond a single file, you'll quickly find yourself wanting to:

- reuse helper functions across multiple scripts,
- share well-tested utilities with collaborators,
- come back to a project after six months and have it still work,
- keep your high-level analysis flow readable rather than buried in setup code.

This guide describes the pattern we recommend for handling all of those at once: the **companion package**.

!!! tip "Caveat: Opinionated Guide"
    This is a very opinionated guide on one specific development workflow which worked out well for me for quite a few years now. Your mileage may vary, and there are many other perfectly good ways to structure research code in Julia! 
    The key point here is, that Julia has its own set of unique tradeoffs. Blindly following a Python or Matlab workflow in Julia is a recipe for frustration.


## The script-only problem

Most research code starts as a single script. Maybe it looks something like this:

```julia
using CSV
using DataFrames
using NetworkDynamics

# load the data
data = CSV.read("my_data.csv", DataFrame)

# build the model
nw = build_network_from_data(data)

# simulate model and generate trajectory
sol = simulate_specific_scenario(nw)

# analyze the trajectory and plot results
results = perform_analysis(sol)
fig = plot_results(results)
```

So where should `build_network_from_data`, `simulate_specific_scenario`, `perform_analysis`, and `plot_results` actually *live*?

You have a few options. You could keep stuffing them into the script (it grows unwieldy). You could put them in a separate `.jl` file and `include` (or `Revise.includet`) it (a step in the right direction). Or you could promote them into a proper Julia **package** — a *companion package* — that lives next to your scripts in the same repository.

## The companion package pattern

A companion package is a small, project-specific Julia package that holds the reusable code for a particular piece of research. The repository looks like this:

```
MyPaperCompanion/
├╴ src/
│  ├╴ MyPaperCompanion.jl        # top-level module file
│  ├╴ network_construction.jl    # subfiles with the actual implementation
│  ├╴ analysis.jl
│  └╴ plotting.jl
├╴ examples/
│  ├╴ scenario_a.jl              # high-level scripts that USE the package
│  └╴ scenario_b.jl
├╴ test/
│  └╴ runtests.jl                # optional file for unit tests
├╴ Project.toml                  # both the env AND the package metadata
└╴ Manifest.toml                 # exact resolved versions
```

The key idea: **scripts in `examples/` import the companion package and call its functions**. The companion package owns the implementation; the scripts own the narrative.

```julia
# examples/scenario_a.jl
using MyPaperCompanion
using CairoMakie

data = MyPaperCompanion.load_default_data()
nw   = MyPaperCompanion.build_network(data)
sol  = MyPaperCompanion.simulate(nw; scenario=:base)
fig  = MyPaperCompanion.plot_results(sol)
```

## Packages vs environments

In [Part II](@ref env-management) you learned that an environment is "any directory with a `Project.toml`". A **package** is a slightly stronger thing: it's a directory with a `Project.toml` *plus* a `src/PackageName.jl` file *plus* a `name`, `uuid`, and `version` field in the `Project.toml`.

In other words: every package is an environment, but not every environment is a package.

```toml
# Project.toml of a package
name = "MyPaperCompanion"
uuid = "12345678-1234-1234-1234-123456789abc"
version = "0.1.0"

[deps]
NetworkDynamics = "..."
DataFrames = "..."

[compat]
NetworkDynamics = "1.2"
julia = "1.11"
```

When you `activate` a package's directory, two things happen at once:

1. The environment (with its dependencies) becomes active.
2. The package itself becomes available — it's automatically `dev`'d (used directly from the source on disk), since it *is* the active project.

So when you're working *inside* `MyPaperCompanion`, you can just write `using MyPaperCompanion` from any script in `examples/` and it works — no separate install step, and any change you make in `src/` is picked up by Revise on the next call.

You can still `dev` *other* packages alongside (e.g. when you're simultaneously fixing a bug in `NetworkDynamics`), and `add` registered packages as usual.

## Setting up your companion package

The fastest way to scaffold a package is to let `Pkg` do it.

```julia-repl
pkg> generate MyPaperCompanion
```

Run this in the parent directory where you want the package to live. It creates the basic `MyPaperCompanion/` folder with a `Project.toml`, `src/MyPaperCompanion.jl`, and the auto-generated UUID. You can then `cd` into the new folder, run `git init`, and add things to your taste.

The auto-generated `src/MyPaperCompanion.jl` looks roughly like:

```julia
module MyPaperCompanion

greet() = print("Hello World!")

end # module MyPaperCompanion
```

This file is the **entry point** of your package. When someone writes `using MyPaperCompanion`, this is the file Julia loads first. Anything else you want in your package needs to be reachable from here.

Now, activate the package's environment, add the dependencies you need (just like in any other environment), and start writing code. Activate it explicitly:

```julia-repl
pkg> activate path/to/MyPaperCompanion

(MyPaperCompanion) pkg> add NetworkDynamics DataFrames
```

Or — and this is the better workflow — open the `MyPaperCompanion` folder in VSCode, set the Julia env to it via the bottom-left status bar, and start the REPL there. From now on, any `examples/*.jl` script you `SHIFT+ENTER` from will run in the package's environment.

## How to organize code inside the package

Once you have more than one file's worth of code, you'll split things across files. There's a convention to this, and it's important to get right.

### One flat module, multiple files

The recommended structure for a research companion package is **a single, flat module** spread across multiple files using `include`.

```julia
# src/MyPaperCompanion.jl
module MyPaperCompanion

using NetworkDynamics   # all dependencies go here, at the top of the main module file
using DataFrames        # loaded dependencies are visible to all included files
using CairoMakie

export build_network, simulate, perform_analysis, plot_results

include("network_construction.jl")
include("analysis.jl")
include("plotting.jl")

end # module
```

That's typically all your top-level file does: `using` the dependencies, declare which names should be visible to users (`export`), and `include` the implementation files in some sensible order.

### What `include` actually does

`include("file.jl")` is essentially **copy-paste of the file's contents at that point**. It's not "import" in the Python sense, and it doesn't create a separate namespace.

This has a few consequences worth internalizing:

- Everything inside `src/` lives in the same module (`MyPaperCompanion`) and can see everything else there. There's no need for relative imports between your own files. I.e. if you define function `foo()` in file A, you can use `foo()` directly in a different function in file B.
- `using` statements at the top of `MyPaperCompanion.jl` apply to *all* included files. You don't (and shouldn't) re-`using NetworkDynamics` at the top of `analysis.jl`.
- **Including the same file twice is almost always a bug.** You'll get every function defined twice, redefinition warnings, and confused behavior. If you find yourself reaching for "include this file from multiple places", stop — you want a function call, not a re-include.
- The order of definitions only matters for types. I.e. you need to define `struct Foo` befor you can define a function `myfun(a::Foo)`. You don't need to worry about the order of function definitions, or the order of `include`s that contain them.

### `export` controls what users see

Names listed in `export` are made available unqualified when someone writes `using MyPaperCompanion`. Names *not* exported are still accessible — just qualified with the module name:

```julia
using MyPaperCompanion

build_network(data)              # exported, no qualification needed
MyPaperCompanion.internal_thing()  # not exported, but still reachable
```

A reasonable rule of thumb: export the few high-level functions a script user actually calls (`build_network`, `simulate`, `plot_results`); don't export internal helpers (`_normalize_admittance_matrix`, `_pivot_dataframe`). Internal names by convention often start with an underscore, but Julia doesn't enforce this.

### Don't reach for nested submodules

Julia's `module ... end` syntax lets you define modules within modules. **You almost never want to do this in a research companion package.**

Nested submodules add a real cost:

- They create separate namespaces, so you have to `using ..ParentModule.X` everywhere.
- They interact awkwardly with `export` — exports from a submodule aren't automatically re-exported from the parent.
- They make the import story for downstream users more confusing.

The legitimate use cases for submodules (extensions, conditional code, large libraries with multiple distinct subsystems) are not what research companion packages need. Stick to one flat module per package; if your codebase grows large enough to *want* submodules, that's usually a signal it should be split into multiple packages instead.

## The day-to-day workflow

Once your companion package exists, your typical session looks like this:

1. Open the package folder in VSCode. Pick its environment via the bottom-left status bar.
2. Start the REPL. (Revise is auto-loaded by the VSCode extension if it's in your global env — see [Part II](@ref env-management).)
3. Open an `examples/*.jl` script. `SHIFT+ENTER` through it block by block.
4. When you find a bug or want a new helper, edit a file in `src/`. Save.
5. Re-run the relevant lines in the script. Revise has already picked up your edit. No restart.

This loop is the productive heart of Julia work. The persistent REPL keeps everything compiled and warm; Revise keeps your library code in sync with what's on disk; the companion package gives all your reusable logic a clean home.

## From companion package to shared library

One of the nicest properties of the companion-package pattern is that it sets you up for a smooth transition when your "private helpers" turn out to be useful to other people too.

The typical arc looks like this:

**Stage 1 — your own research companion.** You're working on a paper. You've put your helpers, your custom component models, and your plotting routines into `MyPaperCompanion`. It's a normal package; it just happens to be *paper-specific* and lives in your paper's repository.

**Stage 2 — coworkers get interested.** Maybe in the course of your paper work you developed a clean implementation of grid-forming droop inverters with LCL filters, and a colleague now wants to use the same models for an unrelated study. The droop inverter implementation inside `MyPaperCompanion` is genuinely reusable — it's not paper-specific at all.

Because that code already lives in a proper package, extracting it is mostly mechanical:
1. Create a new package, e.g. `GridFormingModels`, in its own Git repository.
2. Move the relevant files from `MyPaperCompanion/src/` into `GridFormingModels/src/` and adjust imports/exports.
3. In each research project that needs them, run `] dev path/to/GridFormingModels` (or `] dev https://github.com/your-org/GridFormingModels.jl` once it's pushed).
4. Replace the now-extracted definitions in each companion package with `using GridFormingModels`.

The crucial point: if you'd put everything in loose `include`d files in your scripts directory, this extraction would be a manual untangling exercise. Because you started with package scaffolding from day one, it's just *moving files between two packages*.

At this stage you should probably start adding unit tests. Start with a single small
**`test/runtests.jl`** test file. This is one of the cheapest ways to catch regressions when you `update` your environment. Even just a few `@test`s on the highest-level functions pays for itself. See the [Julia `Test` stdlib documentation](https://docs.julialang.org/en/v1/stdlib/Test/) for what you can write there.

**Stage 3 — going public.** If `GridFormingModels` proves genuinely useful beyond your group, you can register it in the [General registry](https://github.com/JuliaRegistries/General) and downstream users can `add` it like any other package. By that point, you've already worked out the API shape through real internal use, you have tests, and you have multiple downstream "consumers" (your various research projects) keeping you honest about backward compatibility.

You don't have to plan for any of this in advance. The pattern just doesn't get in the way when it happens — which, in our experience, is the most you can ask for from a code-organization scheme.

## Gotchas and tips

A few things that trip people up:

**The companion package is its own environment.** Don't accidentally activate the global env in a session that was meant to be working on the companion. (See the [warning about runtime environment switching](@ref start-in-project) in Part II.)

**Don't worry about dependency hygiene in a research companion.** For *real* libraries you intend others to use, you want the package's environment to be lean — every dependency you add becomes a constraint imposed on every downstream user. A research companion is different. It's not meant to be consumed by anyone; it's a private home for your project's code. So go ahead and add `CairoMakie`, `DataFrames`, `CSV`, whatever your scripts need, directly to the companion's `Project.toml`. Putting them in a nested environment under `examples/` would introduce back unnecessary complexity. The hygiene principle reasserts itself the moment you start spinning out genuinely shared packages — those should keep their dependencies tight.

**Pay attention to your REPL prompt color.** If the `julia>` prompt turns yellow, Revise is yelling at you. Maybe you introduce a syntax error while modifying source files or something like that. You can call `Revise.retry()` to retry loading or at least resurfacing the error. If you ignore it, you might not be running the code you think you're running!

**Commit the `Manifest.toml` for research projects.** Yes, it's auto-generated, but committing it is what guarantees a future-you (or a collaborator) can reproduce your numbers. This is good safeguard against accidentally breaking your environment by package updates. For a *library* package intended for downstream use, you'd `.gitignore` it instead — a different situation.

**`dev` other packages locally when you need to.** If you find a bug in `NetworkDynamics` while working on your project, run `dev NetworkDynamics` in your companion env. This clones it to `~/.julia/dev/NetworkDynamics` and points your env at the source. Edit, test from your scripts, push fixes upstream, then `free NetworkDynamics` to switch back to the registered version.

**Resist the temptation of nested modules.** When you feel the pull, the right answer is almost always "split into a second package" or "just keep everything flat".

## Where to go from here

You now have the full picture: a working Julia setup ([Part I](@ref julia-setup)), an understanding of environments and `Pkg` ([Part II](@ref env-management)), and a structure for organizing real research code (this part). From here:

- The [Modern Julia Workflows blog](https://modernjuliaworkflows.org/) covers many topics adjacent to what's here — testing, CI, profiling, debugging.
- The [Pkg.jl documentation](https://pkgdocs.julialang.org/) is excellent when you need depth on a specific Pkg feature. Specific, notably and quite recent features of Pkg are the [`[sources]`-sections](https://pkgdocs.julialang.org/v1/toml-files/#The-%5Bsources%5D-section) and the [`[workspace]`-section](https://pkgdocs.julialang.org/v1/toml-files/#Workspaces) and how you might [use those features across your packages](https://discourse.julialang.org/t/how-to-use-workspaces-to-test/133670/16).
- Check out how to run your tests automatically in CI with GitHub actions.
- Check out [`Documenter.jl`](https://juliadocs.github.io/Documenter.jl/stable/) for how to write documentation for your package and publish it online.
