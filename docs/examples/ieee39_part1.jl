#=
# [IEEE39 Bus Tutorial - Part I: Model Creation](@id ieee39-part1)

This tutorial can be downloaded as a normal Julia script [here](@__NAME__.jl). #md

This is the first part of a four-part tutorial series for the IEEE 39-bus test system:

- **Part I: Model Creation** (this tutorial) - Build the network structure with buses, lines, and components
- **Part II: Initialization** - Perform power flow calculations and dynamic initialization
- **Part III: Dynamic Simulation** - Run time-domain simulations and analyze system behavior
- **Part IV: Advanced Modeling & Parameter Optimization** - Create custom components and optimize system parameters

In this first part, we'll construct the complete IEEE 39-bus network model using PowerDynamics.jl,
including generators, loads, transmission lines, and control systems.

## System Structure
The system consists of 39 buses (with 10 generators and 19 loads) and 46 branches (12 of which are transformers).

The buses fall into the following categories:
- Junction: pure transient buses without dynamic components
- Load: buses with loads only
- Controlled Machine: buses with controlled machines (generators with AVR and GOV)
- Controlled Machine + Load: buses with controlled machines and loads
- Uncontrolled Machine + Load: buses with uncontrolled machines and loads

For the power flow solution, we have a slack bus, PV buses and PQ buses.

For the dynamic simulation, we will use the following models:
- ZIP Load for loads,
- 6th Order Sauer-Pai Machine and
- AVR Type I and TGOV1 for controlled machines.

## Setup and Data Loading

!!! warning "No Standardized Data Import"
    As of now, PowerDynamics.jl does not support any advanced import mechanisms for power grids.
    Therefore, this tutorial loads the data from some custom CSV files.

First, we'll load the required packages and read the system data from CSV files.
The IEEE 39-bus system data is organized into separate files for different components.
=#

using PowerDynamics
using PowerDynamics.Library
using ModelingToolkitBase
using NetworkDynamics
using DataFrames
using CSV

DATA_DIR = joinpath(pkgdir(PowerDynamics), "docs", "examples", "ieee39data")
nothing #hide #md

#=
The system data is stored in CSV files containing:

!!! details "bus.csv - Bus Configuration Data"
    | Parameter | Description |
    |:----------|:------------|
    | `bus` | Bus number (unique identifier) |
    | `bus_type` | Power flow bus type: "PQ" (load), "PV" (generator), "Slack" (reference) |
    | `category` | Component category: "junction", "load", "ctrld\_machine", "ctrld\_machine\_load", "unctrld\_machine\_load" |
    | `P` | Active power injection [pu] \(positive = generation, negative = load\) |
    | `Q` | Reactive power injection [pu] \(positive = generation, negative = load\) |
    | `V` | Voltage magnitude [pu] \(for PV and Slack buses\) |
    | `base_kv` | Base voltage level [kV] |
    | `has_load` | Boolean flag indicating presence of load component |
    | `has_gen` | Boolean flag indicating presence of generator component |
    | `has_avr` | Boolean flag indicating presence of automatic voltage regulator |
    | `has_gov` | Boolean flag indicating presence of turbine governor |

!!! details "branch.csv - Transmission Line and Transformer Data"
    | Parameter | Description |
    |:----------|:------------|
    | `src_bus` | Source bus number |
    | `dst_bus` | Destination bus number |
    | `transformer` | Transformer flag (0 = line, 1 = transformer) |
    | `r_src` | Source end transformation ratio [pu] |
    | `R` | Series resistance [pu] |
    | `X` | Series reactance [pu] |
    | `G_src` | Source end shunt conductance [pu] |
    | `G_dst` | Destination end shunt conductance [pu] |
    | `B_src` | Source end shunt susceptance [pu] |
    | `B_dst` | Destination end shunt susceptance [pu] |

!!! details "load.csv - ZIP Load Model Parameters"
    | Parameter | Description |
    |:----------|:------------|
    | `bus` | Bus number where load is connected |
    | `Pset` | Active power at operation point [pu] |
    | `Qset` | Reactive power at operation point [pu] |
    | `KpZ` | Active power constant impedance fraction |
    | `KqZ` | Reactive power constant impedance fraction |
    | `KpI` | Active power constant current fraction |
    | `KqI` | Reactive power constant current fraction |
    | `KpC` | Active power constant power fraction (1-KpZ-KpI) |
    | `KqC` | Reactive power constant power fraction (1-KqZ-KqI) |

    Note: ZIP loads combine constant impedance (Z), constant current (I), and constant power (P) components.

!!! details "machine.csv - Generator (Sauer-Pai Machine) Parameters"
    | Parameter | Description |
    |:----------|:------------|
    | `bus` | Bus number where generator is connected |
    | `Sn` | Machine power rating [MVA] |
    | `V_b` | Bus voltage base [kV] \(unused: duplicates `base_kv` in bus.csv\) |
    | `Vn` | Machine voltage rating [kV] |
    | `R_s` | Stator resistance [pu] |
    | `X_ls` | Stator leakage reactance [pu] |
    | `X_d` | d-axis synchronous reactance [pu] |
    | `X_q` | q-axis synchronous reactance [pu] |
    | `X′_d` | d-axis transient reactance [pu] |
    | `X′_q` | q-axis transient reactance [pu] |
    | `X″_d` | d-axis subtransient reactance [pu] |
    | `X″_q` | q-axis subtransient reactance [pu] |
    | `T′_d0` | d-axis transient time constant [s] |
    | `T′_q0` | q-axis transient time constant [s] |
    | `T″_d0` | d-axis subtransient time constant [s] |
    | `T″_q0` | q-axis subtransient time constant [s] |
    | `H` | Inertia constant [s] |
    | `D` | Direct shaft damping coefficient |

!!! details "avr.csv - Automatic Voltage Regulator (AVR Type I) Parameters"
    | Parameter | Description |
    |:----------|:------------|
    | `bus` | Bus number where AVR-controlled generator is located |
    | `Ka` | Amplifier gain |
    | `Ke` | Field circuit integral deviation |
    | `Kf` | Stabilizer gain |
    | `Ta` | Amplifier time constant [s] |
    | `Tf` | Stabilizer time constant [s] |
    | `Te` | Field circuit time constant [s] |
    | `Tr` | Measurement time constant [s] |
    | `vr_min` | Minimum regulator voltage [pu] |
    | `vr_max` | Maximum regulator voltage [pu] |
    | `E1` | First ceiling voltage [pu] |
    | `Se1` | First ceiling saturation factor |
    | `E2` | Second ceiling voltage [pu] |
    | `Se2` | Second ceiling saturation factor |

!!! details "gov.csv - Turbine Governor (TGOV1) Parameters"
    | Parameter | Description |
    |:----------|:------------|
    | `bus` | Bus number where governor-controlled generator is located |
    | `V_min` | Minimum valve position [pu] |
    | `V_max` | Maximum valve position [pu] |
    | `R` | Governor droop [Machine PU] |
    | `T1` | First transient time constant [s] |
    | `T2` | Second transient time constant [s] |
    | `T3` | Third transient time constant [s] |
    | `DT` | Turbine damping coefficient |
    | `ω_ref` | Reference frequency [pu] |
=#

branch_df = CSV.read(joinpath(DATA_DIR, "branch.csv"), DataFrame)
bus_df = CSV.read(joinpath(DATA_DIR, "bus.csv"), DataFrame)
load_df = CSV.read(joinpath(DATA_DIR, "load.csv"), DataFrame)
machine_df = CSV.read(joinpath(DATA_DIR, "machine.csv"), DataFrame)
avr_df = CSV.read(joinpath(DATA_DIR, "avr.csv"), DataFrame)
gov_df = CSV.read(joinpath(DATA_DIR, "gov.csv"), DataFrame)
nothing #hide #md

#=
System base values follow the IEEE 39-bus standard. `Sbase` and the frequency base are
*global* in PowerDynamics: they are set once, before any model is constructed, and every
component built afterwards bakes them in. Since they are read at construction time, changing
them later has no effect on models that already exist — and by the same token, a script that
changes them should restore the defaults when it is done, which we do at the very end of this
page.
=#

BASE_MVA = 100.0
BASE_FREQ = 60.0
set_Sbase!(BASE_MVA)
set_fbase!(BASE_FREQ)
nothing #hide #md


#=
## Subcomponent Definition
As stated above, our buses fall into 5 different categories.
We will define a "template" for each of those categories and then create the individual buses from those templates.
By doing so, we can reach substantial performance improvements, as we do not have to repeatedly **compile** the same models (the symbolic simplification is quite costly).
Instead, we copy the templates and adjust parameters.

However, before we can define the bus templates, we need to define the individual subcomponents.
Those subcomponents are MTK models and not yet compiled node models. See [Modeling Concepts](@ref) and the [custom bus tutorial](@ref custom-bus).

### Load Model
We use the ZIP load model to represent loads. This model satisfies the [Injector Interface](@ref).
```asciiart
(t) ┌──────────┐
 o──┤ ZIP Load │
    └──────────┘
```

Every parameter that we are going to read from the CSV files is declared **free** by passing
`nothing` instead of a value. A library model normally comes with sensible defaults; passing
`nothing` removes the default for that parameter, so the compiled model carries the parameter
but no value for it. That way the CSV data is the single source of truth: if a column is
missing, initialization fails loudly instead of silently falling back to a library default.
=#
load = ZIPLoad(;
    name=:ZIPLoad,
    Pset=nothing, Qset=nothing,
    KpZ=nothing, KqZ=nothing,
    KpI=nothing, KqI=nothing,
    KpC=nothing, KqC=nothing,
)
nothing #hide #md

#=
### Generator Models

For generators, we use the Sauer-Pai machine model, which is a 6th-order synchronous machine model.
We create two variants:

**Uncontrolled Machine**: No external control inputs for mechanical torque or field voltage.
This model satisfies the [Injector Interface](@ref) directly.
```asciiart
(t) ┌─────────┐
 o──┤ Machine │
    └─────────┘
```
=#

## all machine parameters come from machine.csv, so they are declared free
unset_machine_p = (;
    Sn=nothing, Vn=nothing, R_s=nothing, X_ls=nothing,
    X_d=nothing, X_q=nothing, X′_d=nothing, X′_q=nothing, X″_d=nothing, X″_q=nothing,
    T′_d0=nothing, T′_q0=nothing, T″_d0=nothing, T″_q0=nothing,
    H=nothing, D=nothing,
)

uncontrolled_machine = SauerPaiMachine(;
    τ_m_input=false,  ## No external mechanical torque input
    vf_input=false,   ## No external field voltage input
    name=:machine,
    unset_machine_p...,
)
nothing #hide #md

#=
**Controlled Machine**: Includes automatic voltage regulator (AVR) and turbine governor controls.

The controlled machine is modeled as a **composite injector**. It consists
of 3 subcomponents: the machine, the AVR and the governor.
The AVR receives the voltage magnitude measurement from the terminal of the machine and sets the field voltage.
The governor receives the frequency measurement and sets the mechanical torque.
Together, they satisfy the [Injector Interface](@ref).

```asciiart
      ┌───────────────────────────────┐
      │ CtrldMachine  u_mag_meas      │
      │              ╭─────→────╮     │
      │    ┌─────────┴─┐      ┌─┴───┐ │
  (t) │    │           ├───←──┤ AVR │ │
   o──┼────┤ Sauer-Pai │ vf   └─────┘ │
      │    │ Machine   │ τ_m  ┌─────┐ │
      │    │           ├───←──┤ Gov │ │
      │    └─────────┬─┘      └─┬───┘ │
      │              ╰─────→────╯     │
      │                 ω_meas        │
      └───────────────────────────────┘
```
=#

_machine = SauerPaiMachine(;
    name=:machine,
    unset_machine_p...,
)
_avr = AVRTypeI(;
    name=:avr,
    ceiling_function=:quadratic,
    ## from avr.csv
    Ka=nothing, Ke=nothing, Kf=nothing,
    Ta=nothing, Tf=nothing, Te=nothing, Tr=nothing,
    vr_min=nothing, vr_max=nothing,
    E1=nothing, Se1=nothing, E2=nothing, Se2=nothing,
)
_gov = TGOV1(;
    name=:gov,
    ## from gov.csv
    V_min=nothing, V_max=nothing, R=nothing,
    T1=nothing, T2=nothing, T3=nothing, DT=nothing, ω_ref=nothing,
)

controlled_machine = CompositeInjector(
    [_machine, _avr, _gov],
    name=:ctrld_gen
)
nothing #hide #md

#=
## Bus Template Creation

Now we have all the components (i.e., the MTK models) so we can combine them into full bus models and compile the methods.

### Junction Bus
Pure transmission buses with no generation or load
```asciiart
           ╔══════════════════════╗
           ║ Junction (compiled)  ║
 Network   ║  ┌─────────────────┐ ║
interface  ║  │MTKBus           │ ║
 current ────→│┌──────┐         │ ║
           ║  ││BusBar│(nothing)│ ║
 voltage ←────│└──────┘         │ ║
           ║  └─────────────────┘ ║
           ╚══════════════════════╝
```
=#

#=
Note that the per-unit bases need no attention here. `Sbase` and `ωbase` were set globally at
the top of this page and are picked up by every `compile_bus` call below; they end up on the
bus's `systembase`, and each device's own `Sbase`/`ωbase` are structurally bound to it rather
than being per-device parameters. The one genuinely per-bus base, `busbar₊Vbase`, is set from
`bus.csv` when we instantiate the individual buses further down.
=#

@named junction_bus_template = compile_bus(MTKBus())
junction_bus_template #hide #md

#=
### Load Bus
Buses with only load components

```asciiart
           ╔═════════════════════╗
           ║ Load (compiled)     ║
 Network   ║  ┌────────────────┐ ║
interface  ║  │MTKBus          │ ║
 current ────→│┌──────┐ ┌────┐ │ ║
           ║  ││BusBar├o┤Load│ │ ║
 voltage ←────│└──────┘ └────┘ │ ║
           ║  └────────────────┘ ║
           ╚═════════════════════╝
```
=#

@named load_bus_template = compile_bus(MTKBus(load))
load_bus_template #hide #md

#=
### Generator Bus (Controlled)
Buses with controlled generators (machine + AVR + governor)

```asciiart
            ╔════════════════════════════════════════════════╗
            ║ Ctrld Machine Bus (compiled)                   ║
            ║  ┌───────────────────────────────────────────┐ ║
            ║  │MTKBus                                     │ ║
            ║  │         ┌───────────────────────────────┐ │ ║
  Network   ║  │         │CtrldMachine  ╭─────→────╮     │ │ ║
 interface  ║  │         │    ┌─────────┴─┐      ┌─┴───┐ │ │ ║
  current ────→│┌──────┐ │    │           ├───←──┤ AVR │ │ │ ║
            ║  ││BusBar├o┼────┤ Sauer-Pai │      └─────┘ │ │ ║
  voltage ←────│└──────┘ │    │ Machine   │      ┌─────┐ │ │ ║
            ║  │         │    │           ├───←──┤ Gov │ │ │ ║
            ║  │         │    └─────────┬─┘      └─┬───┘ │ │ ║
            ║  │         │              ╰─────→────╯     │ │ ║
            ║  │         └───────────────────────────────┘ │ ║
            ║  └───────────────────────────────────────────┘ ║
            ╚════════════════════════════════════════════════╝
```
=#

@named ctrld_machine_bus_template = compile_bus(
    MTKBus(controlled_machine);
)
ctrld_machine_bus_template #hide #md

#=
### Generator + Load Bus (Controlled)
Buses with both controlled generators and loads

```asciiart
            ╔═════════════════════════════════════════════════╗
            ║ Ctrld Machine Load Bus (compiled)               ║
            ║  ┌────────────────────────────────────────────┐ ║
            ║  │MTKBus    ┌───────────────────────────────┐ │ ║
            ║  │          │CtrldMachine  ╭─────→────╮     │ │ ║
            ║  │          │    ┌─────────┴─┐      ┌─┴───┐ │ │ ║
            ║  │          │    │           ├───←──┤ AVR │ │ │ ║
  Network   ║  │        ┌─┼────┤ Sauer-Pai │      └─────┘ │ │ ║
 interface  ║  │        │ │    │ Machine   │      ┌─────┐ │ │ ║
  current ────→│┌──────┐│ │    │           ├───←──┤ Gov │ │ │ ║
            ║  ││BusBar├o │    └─────────┬─┘      └─┬───┘ │ │ ║
  voltage ←────│└──────┘│ │              ╰─────→────╯     │ │ ║
            ║  │        │ └───────────────────────────────┘ │ ║
            ║  │        │ ┌──────┐                          │ ║
            ║  │        └─┤ Load │                          │ ║
            ║  │          └──────┘                          │ ║
            ║  └────────────────────────────────────────────┘ ║
            ╚═════════════════════════════════════════════════╝
```
=#

@named ctrld_machine_load_bus_template = compile_bus(
    MTKBus(controlled_machine, load);
)
ctrld_machine_load_bus_template #hide #md

#=
### Generator + Load Bus (Uncontrolled)
Buses with uncontrolled generators and loads

```asciiart
            ╔════════════════════════════════╗
            ║ Unctr. Ma. Load Bus (compiled) ║
            ║  ┌────────────────────────┐    ║
  Network   ║  │MTKBus      ┌─────────┐ │    ║
 interface  ║  │          ┌─┤ Machine │ │    ║
  current ────→│ ┌──────┐ │ └─────────┘ │    ║
            ║  │ │BusBar├─o             │    ║
  voltage ←────│ └──────┘ │ ┌──────┐    │    ║
            ║  │          └─┤ Load │    │    ║
            ║  │            └──────┘    │    ║
            ║  └────────────────────────┘    ║
            ╚════════════════════════════════╝
```
=#

@named unctrld_machine_load_bus_template = compile_bus(
    MTKBus(uncontrolled_machine, load);
)
unctrld_machine_load_bus_template #hide #md

#=
## Bus Instantiation and Parameter Setting

Now we create the actual bus instances by copying templates and applying
specific parameters from the CSV data files.
=#

## Helper function to apply CSV parameters to bus components
function apply_csv_params!(bus, table, bus_index)
    row_idx = findfirst(table.bus .== bus_index)

    ## Apply all parameters except "bus" column
    row = table[row_idx, :]
    for col_name in names(table)
        col_name == "bus" && continue
        ## `V_b` in machine.csv duplicates `base_kv` in bus.csv. The voltage base is a
        ## property of the bus, not of the machine, so we take it from bus.csv below.
        col_name == "V_b" && continue
        set_default!(bus, Regex(col_name*"\$"), row[col_name])
    end
end
nothing #hide #md

#=
For each bus in the system, we:
1. Select the appropriate template based on its category
2. Create a bus instance with the correct vertex index and name
3. Apply component-specific parameters from CSV files
4. Set the power flow model (PQ, PV, or Slack)
=#

busses = []
for row in eachrow(bus_df)
    i = row.bus

    ## Select template based on bus category
    bus = if row.category == "junction"
        compile_bus(junction_bus_template; vidx=i, name=Symbol("bus$i"))
    elseif row.category == "load"
        compile_bus(load_bus_template; vidx=i, name=Symbol("bus$i"))
    elseif row.category == "ctrld_machine"
        compile_bus(ctrld_machine_bus_template; vidx=i, name=Symbol("bus$i"))
    elseif row.category == "ctrld_machine_load"
        compile_bus(ctrld_machine_load_bus_template; vidx=i, name=Symbol("bus$i"))
    elseif row.category == "unctrld_machine_load"
        compile_bus(unctrld_machine_load_bus_template; vidx=i, name=Symbol("bus$i"))
    end

    ## The voltage base is the one per-unit base that is genuinely per bus (this system mixes
    ## 16.5, 138, 230 and 345 kV). It only affects the SI observables (`u_kV`, `P_MW`, …) and
    ## is handed down to the incident lines during initialization.
    set_default!(bus, :busbar₊Vbase, row.base_kv)

    ## Apply component parameters from CSV files
    row.has_load && apply_csv_params!(bus, load_df, i)
    row.has_gen && apply_csv_params!(bus, machine_df, i)
    row.has_avr && apply_csv_params!(bus, avr_df, i)
    row.has_gov && apply_csv_params!(bus, gov_df, i)

    ## Set power flow model based on bus type
    pf_model = if row.bus_type == "PQ"
        pfPQ(P=row.P, Q=row.Q)  ## Load bus: fixed P and Q
    elseif row.bus_type == "PV"
        pfPV(P=row.P, V=row.V)  ## Generator bus: fixed P and V
    elseif row.bus_type == "Slack"
        pfSlack(V=row.V, δ=0)   ## Slack bus: fixed V and angle
    end
    set_pfmodel!(bus, pf_model)

    push!(busses, bus)
end

#=
## Transmission Line Creation

The IEEE 39-bus system includes both transmission lines and transformers,
all modeled using the π-line equivalent circuit model.

The model consists of several layers:
1. The `PiModel`, which satisfies the [Branch Interface](@ref) as it has two terminals
2. The [`MTKLine`](@ref) constructor, which creates a MTK model fulfilling the [MTKLine Interface](@ref)
3. The compiled `EdgeModel` created by calling the [`compile_line`](@ref) constructor
```
       ╔══════════════════════════════════╗
       ║ EdgeModel (compiled)             ║
   src ║ ┌──────────────────────────────┐ ║ dst
vertex ║ │MTKLine                       │ ║ vertex
   u ───→│┌───────┐ ┌────────┐ ┌───────┐│←─── u
       ║ ││LineEnd├o┤ PiLine ├o┤LineEnd││ ║
   i ←───│└───────┘ └────────┘ └───────┘│───→ i
       ║ └──────────────────────────────┘ ║
       ╚══════════════════════════════════╝
```
(We used the `PiLine_fault` model since we plan on simulating short circuits later.)

A line has two voltage bases, `src₊Vbase` and `dst₊Vbase`, one per `LineEnd` — a transformer
between the 16.5 kV generator buses and the 345 kV grid has genuinely different bases on its
two ends. We do not set them here: each `LineEnd` declares that it inherits its `Vbase` from
the bus it is attached to, and that inheritance is resolved during initialization, once the
line knows its neighbours. Part II shows the result.
=#

@named piline_template = compile_line(MTKLine(PiLine_fault(;name=:piline)))

#=
Each transmission element is created by:
1. Instantiating a line from the template with source and destination buses
2. Setting electrical parameters (resistance, reactance, susceptance) from CSV data
=#

branches = []
for row in eachrow(branch_df)
    ## Create line instance with topology
    line = compile_line(piline_template; src=row.src_bus, dst=row.dst_bus)

    ## Apply electrical parameters from CSV data
    for col_name in names(branch_df)
        if col_name ∉ ["src_bus", "dst_bus", "transformer"]
            set_default!(line, Regex(col_name*"\$"), row[col_name])
        end
    end

    push!(branches, line)
end

#=
## Network Assembly

Finally, we combine all buses and transmission lines into a complete network model.
This creates the IEEE 39-bus test system ready for initialization and simulation.
=#

nw = Network(busses, branches)

#=
All models are built, so we restore the global bases to their defaults. Calling the setters
without an argument does that, and it keeps this page's 60 Hz from leaking into whatever is
constructed next in the same session.
=#
set_Sbase!()
set_fbase!()

#=
The network `nw` now contains the complete IEEE 39-bus model structure.
In Part 2 of this tutorial series, we'll initialize this network by solving
the power flow and setting up the dynamic initial conditions.
=#
