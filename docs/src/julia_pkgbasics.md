# [Part II: Environments and Package Management](@id env-management)

This is Part II of our Julia introduction series. If you haven't yet, start with [Part I: Julia Setup](@ref julia-setup) to get Julia and VSCode running.

Many newcomers to Julia stumble over package and environment management. The good news: Julia's package manager (`Pkg`) is genuinely excellent. The bad news: package management is a hard problem, so the tools that solve it can't be entirely simple.

This guide is a high-level overview of the *concepts*, aimed at people from non-computer-science backgrounds.
It is **not** a reference manual — for that, see the [official Pkg documentation](https://pkgdocs.julialang.org/).

## Why care about environments?

Reusing code is a core pillar of software development. You don't write your own matrix multiplication — you use `LinearAlgebra`. You don't write your own ODE solver — you use `OrdinaryDiffEq`. You don't write your own plotting code — you use `Makie`.

The catch: all of those packages are in flux. `OrdinaryDiffEq` improves with every release. So today you might have to use it slightly differently than five years ago. But you still want to **reproduce** results from five years ago while simultaneously working on a brand-new project — without one breaking the other.

This is what environments are for. You can have two folders, one for your hot new paper and one for your old research. Each folder pins down exactly which packages, at exactly which versions, that project uses. Switching between projects becomes a matter of changing folders — not reinstalling everything.

## `Project.toml` and `Manifest.toml`

Concretely, a Julia **environment** is a directory containing two files:

- **`Project.toml`** — short and human-readable. Lists *your* direct dependencies, plus optional `[compat]` bounds (more on those below). This is the file you'll actually edit (or have `Pkg` edit for you).
- **`Manifest.toml`** — a complete snapshot of *every* package in the environment, including transitive dependencies, with exact versions and source information. Auto-generated. **Never edit it by hand.**

A minimal `Project.toml` might look like this:

```toml
[deps]
PowerDynamics = "10c5d281-7529-5104-9c6c-2e54f9d2c00e"

[compat]
PowerDynamics = "1.0"
julia = "1.10"
```

The `[deps]` section says *what* you depend on (with each package's UUID as a unique identifier). The `[compat]` section says *which versions* are acceptable. The `Manifest.toml` then says *exactly which version of every package* (yours and theirs) was selected.

## Why versions matter: semantic versioning

Why should you care about versions? Imagine you're using `NetworkDynamics`. The maintainers occasionally release new versions, and you want them — for example, maybe there was a bug in solution interpolation that just got fixed. You should be able to run `] update` without fear of your code suddenly breaking.

There are roughly three kinds of changes a release can contain:

- **Bug fixes** — small, targeted, shouldn't change behavior except where it was wrong.
- **New features** — additions that don't affect anything you were already doing.
- **Breaking changes** — e.g. a function's argument order was swapped, so old call sites no longer work.

You want bug fixes and features automatically. You want to **opt in** explicitly to breaking changes. **Semantic versioning** ([semver](https://semver.org/)) makes this distinction visible in the version number itself:

```asciiart
    SomePackage@1.2.4
                ╷ ╷ ╷
                │ │ ╰───╴ Patch — small non-breaking fix
                │ ╰─────╴ Minor — non-breaking new feature
                ╰───────╴ Major — breaking change
```

In short: a bump in the *third* number is a fix, the *second* is a feature, the *first* is a break.

!!! note "Julia's special rule for `0.x.y`"
    Julia treats pre-`1.0` versions specially. Below `1.0`, the *second* number `x` is treated as breaking and the *third* number `y` covers everything non-breaking (both fixes and features). This matches the [Julia compat convention](https://pkgdocs.julialang.org/v1/compatibility/) and lets young, rapidly-evolving packages signal breakage without bumping to `1.0` prematurely.

In a `Project.toml`, you express acceptable versions through **compat bounds**:

```toml
[compat]
ModelingToolkit = "10.1.0"
```

This says: *I require **at least** ModelingToolkit `10.1.0`, but I'm also fine with any later non-breaking version.* So `Pkg` is allowed to pick `10.1.1`, `10.2.0`, even `10.999.999` — but **never** `11.0.0`, because the major bump from 10→11 might break your code, and that needs your explicit say-so.

If you don't write any compat bounds, `Pkg` defaults to permissive behavior. This is fine while exploring, but for any project you want to come back to later, you should add bounds. The shorthand command `] compat` helps you set them.

## What does the package manager actually do?

When you write a `Project.toml`, you're describing *constraints* — top-level dependencies and (optionally) version bounds. But each of those packages has its own `Project.toml` with its own dependencies and bounds, and so on recursively.

`Pkg`'s job is to walk this **dependency graph** and find a single set of versions where:

- Each package appears exactly once.
- Each package's version satisfies the compat bounds of every package that depends on it.

```asciiart
(Subset of a) Package Dependency Graph of an Environment:

 Top-level Dependencies  ┃ Transitive Dependencies
 from Project.toml       ┃ resolved into Manifest.toml
          ╭──────────────────────▶───────────────────────╮
 ╔════════╧════════════╗ ┃                               │
 ║ NetworkDynamics@1.2 ╟─────▶────╮                      │
 ╚═════════════════════╝ ┃        │                      │
 ╔════════════╗          ┃ ╭──────┴────────╮     ╭───────┴────────────────╮
 ║ DiffEq@1.3 ╟────▶───────┤ SciMLBase@2.1 ├─▶───┤ OrderedCollections@1.3 │
 ╚════════╤═══╝          ┃ ╰───────────────╯     ╰───────┬────────────────╯
          ╰──────────────────────▶───────────────────────╯
                         ┃
```

Note that each package can appear only once in the resolved graph. If `OrderedCollections` is a dependency of three different packages (as above), the version chosen must satisfy *all three* sets of compat bounds simultaneously. When this is impossible, `Pkg` errors out and tells you what conflicts.

Once a consistent set of versions is found, `Pkg` writes the result to `Manifest.toml`. That's why the manifest is the complete picture for reproducing an environment: it pins every version of every (transitive) dependency. The `Project.toml` only states your top-level intent.

## Working with `Pkg`

Now let's go through the basic commands. Create a folder you want to work in — say `~/tmp/testenv` — and open it in VSCode. Inside, create `test.jl` and start a Julia REPL by hitting `SHIFT + ENTER` somewhere in the file.

### Creating an environment and adding packages

In the REPL, hit `]` to enter Pkg mode. The prompt changes:

```julia-repl
(@v1.12) pkg>
```

The `(@v1.12)` tells you you're currently in the **global** environment for Julia 1.12. Let's activate the current folder instead:

```julia-repl
(@v1.12) pkg> activate .

(testenv) pkg>
```

The prompt changing to `(testenv)` confirms we're in the new environment. Add a small package to test:

```julia-repl
(testenv) pkg> add OrderedCollections
    Updating registry at `~/.julia/registries/General.toml`
   Resolving package versions...
    Updating `path/to/testenv/Project.toml`
  [bac558e1] + OrderedCollections v1.8.1
    Updating `path/to/testenv/Manifest.toml`
  [bac558e1] + OrderedCollections v1.8.1

(testenv) pkg> status
Status `path/to/testenv/Project.toml`
  [bac558e1] OrderedCollections v1.8.1
```

`add` automatically pulls the latest compatible version. The `Project.toml` and `Manifest.toml` are now both updated.

You can use the package in `test.jl`:

```julia
using OrderedCollections
dict = OrderedDict(1 => 2)
dict[1]
```

### Useful variants of `add`

The basic form `add SomePackage` covers most cases, but `add` has a few important variants:

**Add without disturbing the rest of the environment.**
```julia-repl
pkg> add --preserve=all SomePackage
```
This tells `Pkg` not to upgrade or change anything else when resolving. Useful when an environment is in a known-good state and you just want to bolt on one more thing.

**Add a specific version.**
```julia-repl
pkg> add SomePackage@1.2.3
```
Useful for reproducing a known-good setup, or for debugging version conflicts.

**Add directly from a Git URL.**
```julia-repl
pkg> add https://github.com/someuser/SomePackage.jl
```
Useful for unregistered packages.

**Add a specific branch, tag, or commit.**
```julia-repl
pkg> add SomePackage#feature-branch
pkg> add SomePackage#v1.2.3-rc1
pkg> add SomePackage#a1b2c3d
```
Useful for trying out a not-yet-released feature, or pinning to a known commit.

### Other commands you'll use constantly

You can list every Pkg command with `]?` and read the docs for a specific one with e.g. `]?add`. The most important ones beyond `add`:

**`status` (or `st`)** — show the current environment. So central it has a shortcut. `st --outdated` is also very useful: it shows you which packages have newer versions available and, importantly, *what's holding them back*.

**`instantiate`** — given a `Project.toml` and `Manifest.toml`, download and install everything specified. This is the command you run when you clone someone else's repo (or your own from a year ago): it reconstructs the exact environment from the manifest.

**`update` (or `up`)** — upgrade packages within the bounds in `Project.toml`. Without compat bounds set, `update` can move you to wildly newer versions. **Tip**: copy your `Manifest.toml` to `Manifest.toml.bak` before running `update` on an important project. If something breaks, just restore the backup.

**`rm` (or `remove`)** — remove a top-level dependency.

**`activate`** — switch the active environment. With a path: `activate /some/folder` or `activate .` for the current directory. **Without** a path: returns to the global environment. (See the [section on switching environments at runtime](@ref start-in-project) for an important gotcha.)

**`compat`** — set version bounds interactively for the active environment. You can also edit `[compat]` in `Project.toml` directly.

**`dev`** — clone a package's source code and use *that* as the package, instead of the registered, immutable version. This is essential for editing package code, and we'll come back to it in [Part III](@ref research-projects).

## Revise.jl — or how Julia becomes usable

Without [Revise.jl](https://github.com/timholy/Revise.jl), Julia is borderline unusable for any serious work. Revise's job is to **reload code in the background** as you edit files.

The usage pattern: you write a function `do_complex_simulation()` in a file, call it from the REPL, find a bug, edit the function in your editor, save — and your *next* call to `do_complex_simulation()` from the REPL automatically picks up the change. No restart, no recompile of the whole world.

!!! tip "Install Revise globally"
    Revise is one of the few packages that should live in your **global** environment. It's a tool — like your editor — not a project dependency. The VSCode extension auto-loads Revise if it's installed and will offer to install it for you on first launch.

    To check or install manually:
    ```julia-repl
    pkg> activate                # no path: returns to global env

    (@v1.12) pkg> status
    ```
    If `Revise` shows up, you're done. Otherwise, `add Revise`. While you're here, look at what else is in your global environment — there's rarely a good reason to have many packages there. Restart Julia after adding Revise.

Revise can't track every package on your system — that would be far too much. Instead, it tracks any package you've `dev`'d, plus files you've explicitly told it to watch via `Revise.includet("myfile.jl")` (note the **t** — "include and track").

We'll dig into the `dev` workflow in [Part III](@ref research-projects). For now, the basic story is: when you're hacking on `NetworkDynamics` itself (or any other package), running `dev NetworkDynamics` clones its source to `~/.julia/dev/NetworkDynamics`, makes your environment use *that copy*, and Revise watches the files. Edits → reload → rerun, all in the same session.

Sometimes, Revise can't pick up a change. The `julia >` prompt of the REPL will turn yellow. Sometimes it helps to call
```julia-repl
julia> Revise.retry()
```
explicitly to force revise to retry and -- at the very least -- redisplay the error.

!!! tip "Restart the REPL after big environment changes"
    Switching a package from `add`'d to `dev`'d (or vice versa), or running a major `update`, are good moments to close the REPL and start fresh. Julia's loaded-package state is per-session and can become inconsistent with the on-disk environment after large changes.

## [Pitfall: don't switch active environments at runtime](@id start-in-project)

This is a very common foot-gun for newcomers. The recommendation is simple:

> **Always start Julia in the project environment you intend to use, and avoid `activate`-ing a different environment in the middle of a session.**

In VSCode, click the `Julia env: ...` indicator at the bottom-left and pick your project folder *before* starting the REPL. From the command line, use `julia --project=/path/to/folder` or `julia --project=@.` (the latter activates the closest enclosing `Project.toml`).

Why does this matter? Julia **stacks** environments — your active project env sits on top of the global env. This is what makes globally-installed Revise visible from any project. But it also means the **resolved dependency graph differs** between "I started Julia in my project" and "I started Julia globally and then activated my project".

Consider a concrete scenario. Your global environment has `Revise` installed, and `Revise` happens to depend on `OrderedCollections`. Last time the global env was resolved, that pinned `OrderedCollections@1.2`:


```asciiart
 ╔════════╗                                      ╭────────────────────────╮
 ║ Revise ╟────────────────────▶─────────────────┤ OrderedCollections@1.2 │
 ╚════════╝                                      ╰────────────────────────╯
```

Now you start Julia globally. Either through `startup.jl` or VSCode's auto-Revise, `using Revise` runs — which loads `OrderedCollections@1.2` into the session. However, you actual Project was previously resolved (and precompiled!) for `OrderedCollections@1.3`.

A Julia session can't unload a package. Per semver, 1.2 and 1.3 are compatible, so things will *appear to work*. But your entire precompile cache for everything downstream of `OrderedCollections` — including big packages like `DifferentialEquations` — was built against `@1.3` and is now invalidated. So using those packages in this environment can lead to huge amounts of redudant precompilation which is slow and surprising:

```asciiart
 ╔════════╗                                      ╭────────────────────────╮
 ║ Revise ╟────────────────────▶─────────────────┤ OrderedCollections@1.2 │
 ╚════════╝                                      ╰─┬──┬──┬────────────────╯
          ╭──────────────────────▶─────────────────╯  │  │
 ╔════════╧════════════╗                     ╭────▶───╯  │
 ║ NetworkDynamics@1.2 ╟─────▶────╮          │ ╭─────────╯
 ╚═════════════════════╝          │          │ │
 ╔════════════╗            ╭──────┴────────╮ │ │ ╭ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─╮
 ║ DiffEq@1.3 ╟────▶───────┤ SciMLBase@2.1 ├─╯ │   OrderedCollections@1.3
 ╚════════╤═══╝            ╰───────────────╯   │ ╰ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─╯
          ╰──────────────────────▶─────────────╯

```


If instead you start Julia *directly in the project environment*, the manifest pins `@1.3` from the start, Revise's reference to `OrderedCollections` resolves into that same `@1.3`, and only Revise's small precompile cache (in the global env) gets invalidated:
That's much better!

The takeaway: pick the right environment **before** starting Julia, not after.

## A note on auto-precompilation

By default, `Pkg` automatically precompiles packages after every operation that changes the environment (`add`, `update`, `rm`, ...). This is usually what you want — but if you're doing a sequence of related operations in one go, the repeated precompiles add up.

You can disable auto-precompilation for the duration of a Julia session by setting

```julia-repl
julia> ENV["JULIA_PKG_PRECOMPILE_AUTO"] = 0
```

before doing your environment changes, then triggering a single final precompile via `] precompile` once you're done. See the [Pkg precompilation reference](https://pkgdocs.julialang.org/v1/api/#Pkg.precompile) for details.

For the official, exhaustive Pkg documentation, see the [Pkg.jl docs](https://pkgdocs.julialang.org/). The community-run [Modern Julia Workflows blog](https://modernjuliaworkflows.org/) is also excellent for practical advice.

## Next Steps
With environments and `Pkg` under your belt, you're ready for [Part III: Structuring Research Projects](@ref research-projects), which covers how to organize a growing codebase into a *companion package* — the pattern we recommend for any non-trivial PowerDynamics work.
