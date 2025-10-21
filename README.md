![PowerDynamics Banner](./docs/src/assets/banner-dark.png#gh-dark-mode-only)
![PowerDynamics Banner](./docs/src/assets/banner.png#gh-light-mode-only)

[![codecov](https://codecov.io/gh/JuliaEnergy/PowerDynamics.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaEnergy/PowerDynamics.jl)
[![Stable Docs](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliaenergy.github.io/PowerDynamics.jl/stable/)
[![Dev Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://juliaenergy.github.io/PowerDynamics.jl/dev/)

PowerDynamics.jl: An Open-Source Framework Written in Julia for Dynamic Power Grid Modeling and Analysis.
The main idea of this package is to turn dynamic power grid models into a right-hand-side function of a DAE system, which can then be solved using [OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl).

The main features of PowerDynamics are:
- **[ModelingToolkit.jl](https://github.com/SciML/ModelingToolkit.jl) based Frontend**: Component models such as machines, controllers, power lines, ... are all defined as symbolic, equation-based MTK models.
- **[NetworkDynamics.jl](https://github.com/JuliaDynamics/NetworkDynamics.jl) based Backend**: PD.jl generates Julia code for each of the MTK models as NetworkDynamics.jl compatible vertex- and edge models. NetworkDynamics is used to interconnect those bus and edge models.
- **Component-Based Initialisation**: Annotate your dynamic models with simplified powerflow models. Solve the powerflow and initialize free dynamic states and parameters of the node and edge models to fit the powerflow result.
- **Symbolic Indexing and Observables**: Use your typical SciML syntax like `plot(sol, idxs=VIndex(2,:busbar₊u_mag))` to access arbitrary states, parameters, observables of your MTK models (in this case, the variable named `busbar₊u_mag` of vertex model 2).

... and there is much more! Please check out the [Documentation](https://juliaenergy.github.io/PowerDynamics.jl/stable/).

## Citation

If you use PowerDynamics.jl in your research publications, [please cite our paper](https://www.sciencedirect.com/science/article/pii/S2352711021001345).

```latex
@article{PowerDynamics2022,
  title={PowerDynamics.jl--An experimentally validated open-source package for the dynamical analysis of power grids},
  author={Plietzsch, Anton and Kogler, Raphael and Auer, Sabine and Merino, Julia and Gil-de-Muro, Asier and Li{\ss}e, Jan and Vogel, Christina and Hellmann, Frank},
  journal = {SoftwareX},
  volume = {17},
  pages = {100861},
  year = {2022},
  publisher={Elsevier}
}
```

## Funding
Development of this project was in part funded by the *German Federal Ministry for Economic Affairs and Climate Action* as part of the *OpPoDyn*-Project ([Project ID 01258425/1](https://www.enargus.de/pub/bscw.cgi/?op=enargus.eps2&q=%2201258425/1%22), 2024-2027).

<img src="docs/src/assets/bmwk_logo_en.svg" width="300"/>
