
[![Build Status](https://travis-ci.org/JuliaEnergy/PowerDynamics.jl.svg?branch=master)](https://travis-ci.org/JuliaEnergy/PowerDynamics.jl)
[![Chat on Slack.](https://img.shields.io/badge/chat%20on-slack-yellow.svg)](https://julialang.slack.com/messages/CDAGL4T09/)
[![Get your Slack invitation.](https://img.shields.io/badge/get%20invitation-slack-yellow.svg)](https://slackinvite.julialang.org/)
[![Code on Github.](https://img.shields.io/badge/code%20on-github-blue.svg)](https://github.com/JuliaEnergy/PowerDynamics.jl)

# PowerDynamics.jl - Dynamic Power System Analysis in Julia

This package provides all the tools you need to create a dynamic power grid model
and analyze it.

The source code is licensed under [GPLv3](https://www.gnu.org/licenses/gpl-3.0.en.html) and [published on github](https://github.com/JuliaEnergy/PowerDynamics.jl).

These Docs have been built with the following version of the sub-packages:
```@setup versions
using Pkg
env = Pkg.Types.EnvCache()
subpackages = copy(env.project["deps"])
pop!(subpackages, "Reexport") # ignore Reexport here
if haskey(subpackages, "Documenter")
  pop!(subpackages, "Documenter") # ignore Documenter here
end
uuids = map(Base.UUID, values(subpackages))
function printversions()
  for uuid in uuids
    pkg = Pkg.Types.manifest_info(env, uuid)
    println(rpad(" * $(pkg["name"]) ", 30, "."), " $(pkg["version"])")
  end
end
```
```@example versions
printversions() # hide
```

## Installation

The installation can be done via the new package manager. Either use
```
]add PowerDynamics
```
or copy
```Julia
using Pkg; Pkg.add("PowerDynamics")
```

Please note that `PowerDynamics.jl` is a fast developing library whose API is not settled yet.
In order to ensure that your old code will still work in the future while using the latest version of
`PowerDynamics.jl` for your new code, **we strongly recommend the usage of environments**. [Please check out
this video from the introduction of Pkg3, where environments are introduced, too.](https://www.youtube.com/watch?v=HgFmiT5p0zU)

### Compatibility

`PowerDynamics.jl` is written for Julia 1.0 and above.
We will quickly switch to new Julia version as they come out, but support older versions and enable long transition periods for users.
Julia versions 0.x are not supported.

## Usage

Generally, we distinguish three types of user for `PowerDynamics.jl`:
- [Grid Modeler](@ref)
- [Grid Component Developer](@ref)
- [`PowerDynamics.jl` Developer](@ref)

### Grid Modeler

**Your Goal** is to use `PowerDynamics.jl` to model your grid of preference. You don't
want to implement new types of nodes.

We recommend you to choose your favorite example from [`PowerDynamicsExamples`](https://github.com/JuliaEnergy/PowerDynamicsExamples),
read [Node Types](@ref) and try to understand it. That should give you the kickstart you need. If you
have any questions, [contact us](@ref Contact).

### Grid Component Developer

**Your Goal** is to use `PowerDynamics.jl` to develop types of nodes, e.g. new control schemes for inverters or
new descriptions of synchronous machines.

After going through the introduction for a [Grid Modeler](@ref), we recommend that you read
through [Dynamics Types](@ref) and [Custom Node Types](@ref) and try to implement
a new node type for an example grid. With that, you should have all the tools you need.
If you have any questions, [contact us](@ref Contact).

### `PowerDynamics.jl` Developer

**Your Goal** is to extend `PowerDynamics.jl` with new fundamental functionalities.

After going throught the introduction for a [Grid Modeler](@ref) and a [Grid Component Developer](@ref),
read through the code where hopefully all of this documentation will helpful for you.
Afterwards, it's probably best to open an issue explainng the idea you want to implement
and we can discuss how you can transform your idea into a pull request.
