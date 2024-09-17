# Modeling Concepts

## Terminal

The `Terminal`─Connector is an important building block for every model.
It represents a connection point with constant voltage in dq─cordinates `u_r` and `u_i` and enforces the kirchoff constraints `sum(i_r)=0` and `sum(i_i)=0`.

## Modeling of Buses
### Model class `Injector`

An injector is a class of components with a single `Terminal` (called `:terminal`).
Examples for injectors might be Generators, Shunts, Loads.
```
     ┌───────────┐
(t)  │           │
 o───┤  Injector │
     │           │
     └───────────┘
```
The current for injectors is allways in injector convention, i.e. positive currents flow *out* of the injector *towards* the terminal.

### Model: `Busbar`
A busbar is a concrete model used in bus modeling.
It represents the physical connection within a bus, the thing where all injectors and lines attach.

```
           ┌──────────┐
i_lines ──>│          │  (t)
           │  Busbar  ├───o
  u_bus <──│          │
           └──────────┘
```
It has special input/output connectors which handle the network interconnection.

### Model class `BusModel`
A `BusModel` is a clase of models, which contains a single `busbar` (named `:busbar`).
Optionally, it may contain various injectors and or branches (see below).
It is the basis for creating Node models for NetworkDynamics.
```
           ┌───────────────────────────────────┐
           │BusModel             ┌───────────┐ │
           │   ┌──────────┐   o──┤ Generator │ │
i_lines ──────>│          │   │  └───────────┘ │
           │   │  Busbar  ├───o                │
  u_bus <──────│          │   │  ┌───────────┐ │
           │   └──────────┘   o──│ Load      │ │
           │                     └───────────┘ │
           └───────────────────────────────────┘
```
Sometimes it is not possible to connect all injectors directly but instead one needs or wants `Branches` between the busbar and injector terminal.

You can either create a model which satisfy the `BusModel` interface manually. For simple models (plain connections of a few injectors) it is possible to use the convenience method `BusModel(injectors...)` to create the composite model based on provide injector models.

## Line Modeling
### Model class: `Branch`
A branch is the two─port equivalent to an injector.
I needs to have two terminals, one is called `:src`, the other `:dst`.

Examples for branches are: PI─Model branches, dynamic RL branches or transformers.

```
      ┌───────────┐
(src) │           │ (dst)
  o───┤  Branch   ├───o
      │           │
      └───────────┘
```
*Both* ends follow the injector interface, i.e. current leaving the device towards the terminals is always positive.

### Model: `LineEnd`
A `LineEnd` model is very similar to the `BusBar` model.
It represents one end of a transmission line.
```
          ┌───────────┐
 u_bus ──>│           │  (t)
          │ LineEnd   ├───o
i_line <──│           │
          └───────────┘
```
The main differenc beeing the different input/output conventions for the network interface.

### ModelClass: `LineModel`
Similar to the `BusModel`, a `LineModel` is a model class representing a model which can be used to instantiate a transmission line model for NetworkDynamics.

It musst contain two `LineEnd` instances, one called `:src`, one called `:dst`.

```
         ┌────────────────────────────────────────────────┐
         │ BranchModel      ┌──────────┐                  │
         │  ┌─────────┐  o──┤ Branch A │──o  ┌─────────┐  │
 u_bus ────>│ LineEnd │  │  └──────────┘  │  │ LineEnd │<──── u_bus
         │  │  :src   ├──o                o──┤  :dst   │  │
i_line <────│         │  │  ┌──────────┐  │  │         │────> i_line
         │  └─────────┘  o──┤ Branch B │──o  └─────────┘  │
         │                  └──────────┘                  │
         └────────────────────────────────────────────────┘
```

Simple line models, which consist only of valid `Branch` models can be instantiated using the `LineModel(branches...)` constructor.

More complex models can be created manually.
For example if you want to chain multiple branches between the `LineEnds`, for example something like

```
LineEnd(:src) ──o── Transformer ──o── Pi─Line ──o── LineEnd(:dst)
```

## From MTK Models to NetworkDynamics
Valid `LineModels` and `BusModels` can be transformed into so called `Line` and `Bus` objects.

`Line` and `Bus` structs are no MTK models anymore, but rather containers.
Currently, they mainly contain a NetworkDynamics component function (`ODEVertex`, `StaticEdge`).

Eventually, those models will contain more metadata. For example

 - static representation for powerflow,
 - possibly local information about PU system (for transforming parameters between SI/PU),
 - meta information for initialization, for example initialization model or the information which parameters are considered "tunable" in order to initialize the dynamical model


The exact structure here is not clear yet!

The names are also up for disscussion. I am not to happy with the confusing names

```
Bus (Struct)  ⊃ BusModel (MTK)  ⊃ BusBar (MTK)
Line (Struct) ⊃ LineModel (MTK) ⊃ LineEnd(s) (MTK)
```

Also there are some important naming conventions which we might need to rethink:

- injectors *musst* contain terminal named `:terminal`
- in each `BusModel`, there *musst* be a `BusBar` called `:busbar`, which intern has a terminal called `:terminal`
- each branch *musst* have the terminals called `:src` and `:dst`, howver
- the `LineModel` posesses two `LineEnds` which are called `:src` and `:dst`, both of which contain a terminal called `:terminal`.
