```@meta
CurrentModule = OpPoDyn
```

# OpPoDyn

Documentation for [OpPoDyn](https://github.com/JuliaEnergy/OpPoDyn.jl).

## Project Structure
The project is structured as follows

```
OpPoDyn/
├── assets: contains asses for reference tests
├── docs: contains this documentation
├── OpPoDynTesting: helper package for testing, defines test utilities like reference tests
├── Project.toml
├── src: source code
│   ├── Library: submodule for library, all the models live here
│   └── ...
└── test: test code
```
At this stage, this project is meant to be used with the `main` branch from NetworkDynamics.
Unfortunately, it also depends on the unregistered subpackage `OpPoDynTesting` which makes instantiating the environment a bit tricky (because you can neither add `NetworkDynamics#main` nor `OpPoDyntesting#../OpPoDynTesting` without it complaining about the other dependency.
Thanks to the `[sources]` block in `Project.toml` in Julia v1.11, this shouldn'te be a problem anymore.

If you want to use the realse version of Julia v1.10 I suggest to create a new development environment:

```julia
julia> pwd() # make sure you're in the right folder
".../.julia/dev/OpPoDyn"

(v1.10) pkg> activate devenv

(devenv) pkg> dev NetworkDynamics

(devenv) pkg> dev ./OpPoDyntesting

(devenv) pkg> dev .
```

```@index
```

## Funding
Development of this project was in part funded by the *German Federal Ministry for Economic Affairs and Climate Action* as part of the *OpPoDyn*-Project (Project ID 01258425/1, 2024-2027).

```@raw html
<img src="assets/bmwk_logo_en.svg" width="300"/>
```
