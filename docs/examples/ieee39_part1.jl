#=
# [IEEE39 Bus Tutorial - Part I: Model Creation](@id ieee39-part1)

This is the first part of a three-part tutorial series for the IEEE 39-bus test system:

- **Part I: Model Creation** (this tutorial) - Build the network structure with buses, lines, and components
- **Part II: Initialization** - Perform power flow calculations and dynamic initialization
- **Part III: Dynamic Simulation** - Run time-domain simulations and analyze system behavior

In this first part, we'll construct the complete IEEE 39-bus network model using PowerDynamics.jl,
including generators, loads, transmission lines, and control systems.

This script can be downloaded as a normal Julia script [here](@__NAME__.jl). #md

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
using ModelingToolkit
using NetworkDynamics
using DataFrames
using CSV

DATA_DIR = joinpath(pkgdir(PowerDynamics), "docs", "examples", "ieee39data")
nothing #hide

#=
The system data is stored in CSV files containing:

!!! details "bus.csv - Bus Configuration Data"
    | Parameter | Description |
    |:----------|:------------|
    | `bus` | Bus number (unique identifier) |
    | `bus_type` | Power flow bus type: "PQ" (load), "PV" (generator), "Slack" (reference) |
    | `category` | Component category: "junction", "load", "ctrld\_machine", "ctrld\_machine\_load", "unctrld\_machine\_load" |
    | `P` | Active power injection [pu] (positive = generation, negative = load) |
    | `Q` | Reactive power injection [pu] (positive = generation, negative = load) |
    | `V` | Voltage magnitude [pu] (for PV and Slack buses) |
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
    | `V_b` | System voltage basis [kV] |
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
    | `T2` | Second transient time constant [s]  |
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
nothing #hide

#=
System base values follow the IEEE 39-bus standard:
=#

BASE_MVA = 100.0
BASE_FREQ = 60.0
nothing #hide


#=
## Subcomponent Definition
As stated above, our buses fall in 5 different categories.
We will define a "template" for each of those categories and then create the individual buses from those templates.
By doing so, we can reach substantial performance improvements, as we do not have repeatedly **compile** the same models (the symbolic simplification is quite costly).
Instead, we copy the templates and adjust parameters.

However, before we can define the bus templates, we need to define the individual subcomponents.
Those subcomponents are MTK models and not yet compiled node models. See [Modeling Concepts](@ref) and the [custom bus tutorial](@ref custom_bus).

### Load Model
We use the ZIP load model which represents loads. Those loads satisfy the [Injector Interface](@ref).
```
(t) ┌──────────┐
 o──┤ ZIP Load │
    └──────────┘
```
=#
load = ZIPLoad(;name=:ZIPLoad)
nothing #hide

#=
### Generator Models

For generators, we use the Sauer-Pai machine model, which is a 6th-order synchronous machine model.
We create two variants:

**Uncontrolled Machine**: No external control inputs for mechanical torque or field voltage.
This model satisfies the [Injector Interface](@ref) directly.
```
(t) ┌─────────┐
 o──┤ Machine │
    └─────────┘
```
=#

uncontrolled_machine = SauerPaiMachine(;
    τ_m_input=false,  ## No external mechanical torque input
    vf_input=false,   ## No external field voltage input
    name=:machine,
)
nothing #hide

#=
**Controlled Machine**: Includes automatic voltage regulator (AVR) and turbine governor controls.

The controlled machine is modeled as a **composite injector**, it consists
of 3 subcomponents: the machine, the AVR and the governor.
The AVR receives the voltage magnitude measurement from the terminal of the machine and sets the field voltage.
The governor receives the frequency measurement and sets the mechanical torque.
Together, they satisfy the [Injector Interface](@ref).

```
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
)
_avr = AVRTypeI(;
    name=:avr,
    ceiling_function=:quadratic,
)
_gov = TGOV1(; name=:gov,)

controlled_machine = CompositeInjector(
    [_machine, _avr, _gov],
    name=:ctrld_gen
)
nothing # hide

#=
## Bus Template Creation

Now we have all the components (i.e. the MTK models) so we can combine them into full bus models and compile the methods.

### Junction Bus
Pure transmission buses with no generation or load
```
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

@named junction_bus_template = Bus(MTKBus())
strip_defaults!(junction_bus_template)  ## Clear default parameters for manual setting
junction_bus_template #hide

#=
### Load Bus
Buses with only load components

```
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

@named load_bus_template = Bus(MTKBus(load))
strip_defaults!(load_bus_template)
load_bus_template #hide

#=
### Generator Bus (Controlled)
Buses with controlled generators (machine + AVR + governor)

```
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

@named ctrld_machine_bus_template = Bus(
    MTKBus(controlled_machine);
)
strip_defaults!(ctrld_machine_bus_template)
## Set system-wide base values for all generators
set_default!(ctrld_machine_bus_template, r"S_b$", BASE_MVA)
set_default!(ctrld_machine_bus_template, r"ω_b$", 2π*BASE_FREQ)
ctrld_machine_bus_template # hide

#=
### Generator + Load Bus (Controlled)
Buses with both controlled generators and loads

```
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

@named ctrld_machine_load_bus_template = Bus(
    MTKBus(controlled_machine, load);
)
strip_defaults!(ctrld_machine_load_bus_template)
set_default!(ctrld_machine_load_bus_template, r"S_b$", BASE_MVA)
set_default!(ctrld_machine_load_bus_template, r"ω_b$", 2π*BASE_FREQ)
ctrld_machine_load_bus_template #hide

#=
### Generator + Load Bus (Uncontrolled)
Buses with uncontrolled generators and loads

```
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

@named unctrld_machine_load_bus_template = Bus(
    MTKBus(uncontrolled_machine, load);
)
strip_defaults!(unctrld_machine_load_bus_template)
set_default!(unctrld_machine_load_bus_template, r"S_b$", BASE_MVA)
set_default!(unctrld_machine_load_bus_template, r"ω_b$", 2π*BASE_FREQ)
unctrld_machine_load_bus_template #hide

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
        if col_name != "bus"
            set_default!(bus, Regex(col_name*"\$"), row[col_name])
        end
    end
end
nothing #hide

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
        Bus(junction_bus_template; vidx=i, name=Symbol("bus$i"))
    elseif row.category == "load"
        Bus(load_bus_template; vidx=i, name=Symbol("bus$i"))
    elseif row.category == "ctrld_machine"
        Bus(ctrld_machine_bus_template; vidx=i, name=Symbol("bus$i"))
    elseif row.category == "ctrld_machine_load"
        Bus(ctrld_machine_load_bus_template; vidx=i, name=Symbol("bus$i"))
    elseif row.category == "unctrld_machine_load"
        Bus(unctrld_machine_load_bus_template; vidx=i, name=Symbol("bus$i"))
    end

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
    set_metadata!(bus, :pfmodel, pf_model)

    push!(busses, bus)
end

#=
## Transmission Line Creation

The IEEE 39-bus system includes both transmission lines and transformers,
all modeled using the π-line equivalent circuit model.

The model consists of several layers:
1. The `PiModel`, which fulfills the [Branch interface](@ref) as it has two terminals
2. The [`MTKLine`](@ref) constructor, which creates a MTK model fulfilling the [Line Interface](@ref)
3. The compiled `EdgeModel` created by calling the [`Line`](@ref) constructor
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
=#

@named piline_template = Line(MTKLine(PiLine(;name=:piline)))

#=
Each transmission element is created by:
1. Instantiating a line from the template with source and destination buses
2. Setting electrical parameters (resistance, reactance, susceptance) from CSV data
=#

branches = []
for row in eachrow(branch_df)
    ## Create line instance with topology
    line = Line(piline_template; src=row.src_bus, dst=row.dst_bus)

    ## Apply electrical parameters from CSV data
    for col_name in names(branch_df)
        if col_name ∉ ["src_bus", "dst_bus", "transformer"]
            set_default!(line, col_name, row[col_name])
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
The network `nw` now contains the complete IEEE 39-bus model structure.
In Part 2 of this tutorial series, we'll initialize this network by solving
the power flow and setting up the dynamic initial conditions.
=#
