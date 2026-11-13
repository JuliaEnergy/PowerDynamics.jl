# # IEEE39 Bus Tutorial - Part I: Model Creation
#
#
# This is the first part of a four-part tutorial series for the IEEE 39-bus test system:
#
# - **Part I: Model Creation** (this tutorial) - Build the network structure with buses, lines, and components
# - **Part II: Initialization** - Perform power flow calculations and dynamic initialization
# - **Part III: Dynamic Simulation** - Run time-domain simulations and analyze system behavior
# - **Part IV: Advanced Modeling & Parameter Optimization** - Create custom components and optimize system parameters
#
# In this first part, we'll construct the complete IEEE 39-bus network model using PowerDynamics.jl,
# including generators, loads, transmission lines, and control systems.
#
# ## System Structure
# The system consists of 39 buses (with 10 generators and 19 loads) and 46 branches (12 of which are transformers).
#
# The buses fall into the following categories:
# - Junction: pure transient buses without dynamic components
# - Load: buses with loads only
# - Controlled Machine: buses with controlled machines (generators with AVR and GOV)
# - Controlled Machine + Load: buses with controlled machines and loads
# - Uncontrolled Machine + Load: buses with uncontrolled machines and loads
#
# For the power flow solution, we have a slack bus, PV buses and PQ buses.
#
# For the dynamic simulation, we will use the following models:
# - ZIP Load for loads,
# - 6th Order Sauer-Pai Machine and
# - AVR Type I and TGOV1 for controlled machines.
#
# ## Setup and Data Loading
#
# > **No Standardized Data Import**
# >
# > As of now, PowerDynamics.jl does not support any advanced import mechanisms for power grids.
# > Therefore, this tutorial loads the data from some custom CSV files.
#
# First, we'll load the required packages and read the system data from CSV files.
# The IEEE 39-bus system data is organized into separate files for different components.

using PowerDynamics
using PowerDynamics.Library
using ModelingToolkitBase
using NetworkDynamics
using DataFrames
using CSV

DATA_DIR = joinpath(pkgdir(PowerDynamics), "docs", "examples", "ieee39data")

# The system data is stored in CSV files containing:
#
# > **bus.csv - Bus Configuration Data**
# >
# > | Parameter | Description |
# > |:----------|:------------|
# > | `bus` | Bus number (unique identifier) |
# > | `bus_type` | Power flow bus type: "PQ" (load), "PV" (generator), "Slack" (reference) |
# > | `category` | Component category: "junction", "load", "ctrld\_machine", "ctrld\_machine\_load", "unctrld\_machine\_load" |
# > | `P` | Active power injection [pu] \(positive = generation, negative = load\) |
# > | `Q` | Reactive power injection [pu] \(positive = generation, negative = load\) |
# > | `V` | Voltage magnitude [pu] \(for PV and Slack buses\) |
# > | `base_kv` | Base voltage level [kV] |
# > | `has_load` | Boolean flag indicating presence of load component |
# > | `has_gen` | Boolean flag indicating presence of generator component |
# > | `has_avr` | Boolean flag indicating presence of automatic voltage regulator |
# > | `has_gov` | Boolean flag indicating presence of turbine governor |
#
# > **branch.csv - Transmission Line and Transformer Data**
# >
# > | Parameter | Description |
# > |:----------|:------------|
# > | `src_bus` | Source bus number |
# > | `dst_bus` | Destination bus number |
# > | `transformer` | Transformer flag (0 = line, 1 = transformer) |
# > | `r_src` | Source end transformation ratio [pu] |
# > | `R` | Series resistance [pu] |
# > | `X` | Series reactance [pu] |
# > | `G_src` | Source end shunt conductance [pu] |
# > | `G_dst` | Destination end shunt conductance [pu] |
# > | `B_src` | Source end shunt susceptance [pu] |
# > | `B_dst` | Destination end shunt susceptance [pu] |
#
# > **load.csv - ZIP Load Model Parameters**
# >
# > | Parameter | Description |
# > |:----------|:------------|
# > | `bus` | Bus number where load is connected |
# > | `Pset` | Active power at operation point [pu] |
# > | `Qset` | Reactive power at operation point [pu] |
# > | `KpZ` | Active power constant impedance fraction |
# > | `KqZ` | Reactive power constant impedance fraction |
# > | `KpI` | Active power constant current fraction |
# > | `KqI` | Reactive power constant current fraction |
# > | `KpC` | Active power constant power fraction (1-KpZ-KpI) |
# > | `KqC` | Reactive power constant power fraction (1-KqZ-KqI) |
#
#     Note: ZIP loads combine constant impedance (Z), constant current (I), and constant power (P) components.
#
# > **machine.csv - Generator (Sauer-Pai Machine) Parameters**
# >
# > | Parameter | Description |
# > |:----------|:------------|
# > | `bus` | Bus number where generator is connected |
# > | `Sn` | Machine power rating [MVA] |
# > | `V_b` | Bus voltage base [kV] \(unused: duplicates `base_kv` in bus.csv\) |
# > | `Vn` | Machine voltage rating [kV] |
# > | `R_s` | Stator resistance [pu] |
# > | `X_ls` | Stator leakage reactance [pu] |
# > | `X_d` | d-axis synchronous reactance [pu] |
# > | `X_q` | q-axis synchronous reactance [pu] |
# > | `X′_d` | d-axis transient reactance [pu] |
# > | `X′_q` | q-axis transient reactance [pu] |
# > | `X″_d` | d-axis subtransient reactance [pu] |
# > | `X″_q` | q-axis subtransient reactance [pu] |
# > | `T′_d0` | d-axis transient time constant [s] |
# > | `T′_q0` | q-axis transient time constant [s] |
# > | `T″_d0` | d-axis subtransient time constant [s] |
# > | `T″_q0` | q-axis subtransient time constant [s] |
# > | `H` | Inertia constant [s] |
# > | `D` | Direct shaft damping coefficient |
#
# > **avr.csv - Automatic Voltage Regulator (AVR Type I) Parameters**
# >
# > | Parameter | Description |
# > |:----------|:------------|
# > | `bus` | Bus number where AVR-controlled generator is located |
# > | `Ka` | Amplifier gain |
# > | `Ke` | Field circuit integral deviation |
# > | `Kf` | Stabilizer gain |
# > | `Ta` | Amplifier time constant [s] |
# > | `Tf` | Stabilizer time constant [s] |
# > | `Te` | Field circuit time constant [s] |
# > | `Tr` | Measurement time constant [s] |
# > | `vr_min` | Minimum regulator voltage [pu] |
# > | `vr_max` | Maximum regulator voltage [pu] |
# > | `E1` | First ceiling voltage [pu] |
# > | `Se1` | First ceiling saturation factor |
# > | `E2` | Second ceiling voltage [pu] |
# > | `Se2` | Second ceiling saturation factor |
#
# > **gov.csv - Turbine Governor (TGOV1) Parameters**
# >
# > | Parameter | Description |
# > |:----------|:------------|
# > | `bus` | Bus number where governor-controlled generator is located |
# > | `V_min` | Minimum valve position [pu] |
# > | `V_max` | Maximum valve position [pu] |
# > | `R` | Governor droop [Machine PU] |
# > | `T1` | First transient time constant [s] |
# > | `T2` | Second transient time constant [s] |
# > | `T3` | Third transient time constant [s] |
# > | `DT` | Turbine damping coefficient |
# > | `ω_ref` | Frequency setpoint [pu] |

branch_df = CSV.read(joinpath(DATA_DIR, "branch.csv"), DataFrame)
bus_df = CSV.read(joinpath(DATA_DIR, "bus.csv"), DataFrame)
load_df = CSV.read(joinpath(DATA_DIR, "load.csv"), DataFrame)
machine_df = CSV.read(joinpath(DATA_DIR, "machine.csv"), DataFrame)
avr_df = CSV.read(joinpath(DATA_DIR, "avr.csv"), DataFrame)
gov_df = CSV.read(joinpath(DATA_DIR, "gov.csv"), DataFrame)

# System base values follow the IEEE 39-bus standard. `Sbase` and the frequency base are
# *global* in PowerDynamics: they are set once, before any model is constructed, and every
# component built afterwards bakes them in. Since they are read at construction time, changing
# them later has no effect on models that already exist — and by the same token, a script that
# changes them should restore the defaults when it is done, which we do at the very end of this
# page.

BASE_MVA = 100.0
BASE_FREQ = 60.0
set_Sbase!(BASE_MVA)
set_fbase!(BASE_FREQ)

# ## Device Models
# Each bus is a busbar plus whatever devices the data lists for it. Before we build the buses,
# we look at the device models we are going to use. Those are MTK models and not yet compiled
# node models. See Modeling Concepts and the custom bus tutorial.
#
# ### Load Model
# We use the `ZIPLoad` model to represent loads. This model satisfies the Injector Interface.
# ```asciiart
# (t) ┌──────────┐
#  o──┤ ZIP Load │
#     └──────────┘
# ```
#
# ### Generator Models
# For generators, we use the `SauerPaiMachine`, a 6th-order synchronous machine model.
# It appears in two variants.
#
# **Uncontrolled Machine**: No external control inputs for mechanical torque or field voltage.
# This model satisfies the Injector Interface directly.
# ```asciiart
# (t) ┌─────────┐
#  o──┤ Machine │
#     └─────────┘
# ```
#
# **Controlled Machine**: Includes an automatic voltage regulator (`AVRTypeI`) and a
# turbine governor (`TGOV1`).
#
# The controlled machine is modeled as a **composite injector**. It consists
# of 3 subcomponents: the machine, the AVR and the governor.
# The AVR receives the voltage magnitude measurement from the terminal of the machine and sets the field voltage.
# The governor receives the frequency measurement and sets the mechanical torque.
# Together, they satisfy the Injector Interface.
#
# ```asciiart
#       ┌───────────────────────────────┐
#       │ CtrldMachine  u_mag_meas      │
#       │              ╭─────→────╮     │
#       │    ┌─────────┴─┐      ┌─┴───┐ │
#   (t) │    │           ├───←──┤ AVR │ │
#    o──┼────┤ Sauer-Pai │ vf   └─────┘ │
#       │    │ Machine   │ τ_m  ┌─────┐ │
#       │    │           ├───←──┤ Gov │ │
#       │    └─────────┬─┘      └─┬───┘ │
#       │              ╰─────→────╯     │
#       │                 ω_meas        │
#       └───────────────────────────────┘
# ```
#
# ## Bus Structures
# Combining those devices, the buses in the data come in five structures.
# We don't need to set them up ahead of time; they emerge when we build each bus from the
# devices listed for it. Still, it helps to know what is coming.
#
# **Junction Bus**: pure transmission buses with no generation or load.
# ```asciiart
#            ╔══════════════════════╗
#            ║ Junction (compiled)  ║
#  Network   ║  ┌─────────────────┐ ║
# interface  ║  │MTKBus           │ ║
#  current ────→│┌──────┐         │ ║
#            ║  ││BusBar│(nothing)│ ║
#  voltage ←────│└──────┘         │ ║
#            ║  └─────────────────┘ ║
#            ╚══════════════════════╝
# ```
#
# **Load Bus**: buses with only a load.
# ```asciiart
#            ╔═════════════════════╗
#            ║ Load (compiled)     ║
#  Network   ║  ┌────────────────┐ ║
# interface  ║  │MTKBus          │ ║
#  current ────→│┌──────┐ ┌────┐ │ ║
#            ║  ││BusBar├o┤Load│ │ ║
#  voltage ←────│└──────┘ └────┘ │ ║
#            ║  └────────────────┘ ║
#            ╚═════════════════════╝
# ```
#
# **Generator Bus (Controlled)**: buses with a controlled generator (machine + AVR + governor).
# ```asciiart
#             ╔════════════════════════════════════════════════╗
#             ║ Ctrld Machine Bus (compiled)                   ║
#             ║  ┌───────────────────────────────────────────┐ ║
#             ║  │MTKBus                                     │ ║
#             ║  │         ┌───────────────────────────────┐ │ ║
#   Network   ║  │         │CtrldMachine  ╭─────→────╮     │ │ ║
#  interface  ║  │         │    ┌─────────┴─┐      ┌─┴───┐ │ │ ║
#   current ────→│┌──────┐ │    │           ├───←──┤ AVR │ │ │ ║
#             ║  ││BusBar├o┼────┤ Sauer-Pai │      └─────┘ │ │ ║
#   voltage ←────│└──────┘ │    │ Machine   │      ┌─────┐ │ │ ║
#             ║  │         │    │           ├───←──┤ Gov │ │ │ ║
#             ║  │         │    └─────────┬─┘      └─┬───┘ │ │ ║
#             ║  │         │              ╰─────→────╯     │ │ ║
#             ║  │         └───────────────────────────────┘ │ ║
#             ║  └───────────────────────────────────────────┘ ║
#             ╚════════════════════════════════════════════════╝
# ```
#
# **Generator + Load Bus (Controlled)**: buses with both a controlled generator and a load.
# ```asciiart
#             ╔═════════════════════════════════════════════════╗
#             ║ Ctrld Machine Load Bus (compiled)               ║
#             ║  ┌────────────────────────────────────────────┐ ║
#             ║  │MTKBus    ┌───────────────────────────────┐ │ ║
#             ║  │          │CtrldMachine  ╭─────→────╮     │ │ ║
#             ║  │          │    ┌─────────┴─┐      ┌─┴───┐ │ │ ║
#             ║  │          │    │           ├───←──┤ AVR │ │ │ ║
#   Network   ║  │        ┌─┼────┤ Sauer-Pai │      └─────┘ │ │ ║
#  interface  ║  │        │ │    │ Machine   │      ┌─────┐ │ │ ║
#   current ────→│┌──────┐│ │    │           ├───←──┤ Gov │ │ │ ║
#             ║  ││BusBar├o │    └─────────┬─┘      └─┬───┘ │ │ ║
#   voltage ←────│└──────┘│ │              ╰─────→────╯     │ │ ║
#             ║  │        │ └───────────────────────────────┘ │ ║
#             ║  │        │ ┌──────┐                          │ ║
#             ║  │        └─┤ Load │                          │ ║
#             ║  │          └──────┘                          │ ║
#             ║  └────────────────────────────────────────────┘ ║
#             ╚═════════════════════════════════════════════════╝
# ```
#
# **Generator + Load Bus (Uncontrolled)**: buses with an uncontrolled generator and a load.
# ```asciiart
#             ╔════════════════════════════════╗
#             ║ Unctr. Ma. Load Bus (compiled) ║
#             ║  ┌────────────────────────┐    ║
#   Network   ║  │MTKBus      ┌─────────┐ │    ║
#  interface  ║  │          ┌─┤ Machine │ │    ║
#   current ────→│ ┌──────┐ │ └─────────┘ │    ║
#             ║  │ │BusBar├─o             │    ║
#   voltage ←────│ └──────┘ │ ┌──────┐    │    ║
#             ║  │          └─┤ Load │    │    ║
#             ║  │            └──────┘    │    ║
#             ║  └────────────────────────┘    ║
#             ╚════════════════════════════════╝
# ```

# ## Model Cache
# Getting from model equations to a compiled bus takes three steps:
# 1. Construct the MTK models of the individual devices, such as the machine or the AVR. This is
#    surprisingly costly, because it creates a lot of symbolic terms.
# 2. Combine the models into composite systems, for example with `CompositeInjector` or
#    `MTKBus`.
# 3. Compile the bus with `compile_bus`. NetworkDynamics simplifies the equations of the
#    whole bus and generates Julia code for them.
#
# Most buses in a grid share their structure. The eight buses with a controlled machine in this
# system have the same equations and differ only in their parameter values. Repeating steps 1
# and 3 for each of them would be wasted work. There are two ways to avoid that.
#
# **Templates**: build one bus per structure up front and create every other bus as a copy of
# it with `compile_bus(template; vidx, name)`, adjusting the parameters afterwards. Copying is the
# fastest way to get a bus. The catch is that you have to know the structures before you build
# the grid: sort the buses into categories, build a template for each category and pick the
# right one for every bus.
#
# **Model cache**: build every bus from scratch, as if it were the only one, and let
# NetworkDynamics remember the expensive steps. Device models defined with `@cached_component`,
# like all library models used here, cover step 1: constructing a device again with the same
# options returns the stored model with the new parameter values. `compile_bus` covers step 3:
# for a bus with the same structure as an earlier one, it reuses the simplified equations and
# the generated code. This is a bit slower than copying a template, but the code simply follows
# the data.
#
# What counts as the same structure are the equations, variables and model options (such as
# whether the machine has a field voltage input). Parameter values don't matter. If you are
# curious what was reused, `NetworkDynamics.mtk_cache_stats()` lists the cache hits.
#
# This tutorial uses the cache. It is off by default, so we switch it on.

NetworkDynamics.set_mtk_cache!(:on)

# ## Bus Construction
# The device parameters are spread over the CSV files, one row per bus. A small helper collects
# the row of a bus into a `NamedTuple`, so we can splat it right into the model constructor.

function csv_params(table, bus)
    row = table[findfirst(==(bus), table.bus), :]
    # `V_b` in machine.csv duplicates `base_kv` in bus.csv. The voltage base is a
    # property of the bus, not of the machine, so we take it from bus.csv below.
    cols = filter(c -> c ∉ ("bus", "V_b"), names(table))
    NamedTuple(Symbol(c) => row[c] for c in cols)
end

# Now we go through `bus.csv` and build every bus from the devices the flags `has_gen`,
# `has_avr`, `has_gov` and `has_load` ask for:
# 1. Build the device models with their parameters from the CSV files
# 2. Combine them in an `MTKBus`
# 3. Compile the bus with `compile_bus`, passing the vertex index, the power flow model
#    (PQ, PV or Slack) and the voltage base
#
# The voltage base is the one per-unit base that is genuinely per bus (this system mixes 16.5,
# 138, 230 and 345 kV). It only affects the SI observables (`u_kV`, `P_MW`, …) and is handed
# down to the incident lines during initialization. `Sbase` and `ωbase` need no attention: they
# were set globally above and are picked up by every model we construct.
#
# Compiling a bus is the expensive step. Since the cache is on, only the first bus of each
# structure is actually compiled; all others reuse it.

busses = []
for row in eachrow(bus_df)
    i = row.bus
    injectors = []

    if row.has_gen
        machine = SauerPaiMachine(;
            name=:machine,
            vf_input=row.has_avr,   ## field voltage set by an AVR
            τ_m_input=row.has_gov,  ## mechanical torque set by a governor
            csv_params(machine_df, i)...,
        )
        if row.has_avr || row.has_gov
            # machine and controllers form one composite injector
            ctrld_gen = Any[machine]
            if row.has_avr
                avr = AVRTypeI(; name=:avr, ceiling_function=:quadratic, csv_params(avr_df, i)...)
                push!(ctrld_gen, avr)
            end
            if row.has_gov
                gov = TGOV1(; name=:gov, csv_params(gov_df, i)...)
                push!(ctrld_gen, gov)
            end
            push!(injectors, CompositeInjector(ctrld_gen; name=:ctrld_gen))
        else
            push!(injectors, machine)
        end
    end

    if row.has_load
        load = ZIPLoad(; name=:ZIPLoad, csv_params(load_df, i)...)
        push!(injectors, load)
    end

    pf_model = if row.bus_type == "PQ"
        pfPQ(P=row.P, Q=row.Q)  ## Load bus: fixed P and Q
    elseif row.bus_type == "PV"
        pfPV(P=row.P, V=row.V)  ## Generator bus: fixed P and V
    elseif row.bus_type == "Slack"
        pfSlack(V=row.V, δ=0)   ## Slack bus: fixed V and angle
    end

    bus = compile_bus(MTKBus(injectors...);
        vidx=i, name=Symbol("bus$i"), pf=pf_model, Vbase=row.base_kv)
    push!(busses, bus)
end

# ## Transmission Line Creation
#
# The IEEE 39-bus system includes both transmission lines and transformers,
# all modeled using the π-line equivalent circuit model.
#
# The model consists of several layers:
# 1. The `PiModel`, which satisfies the Branch Interface as it has two terminals
# 2. The `MTKLine` constructor, which creates a MTK model fulfilling the MTKLine Interface
# 3. The compiled `EdgeModel` created by calling the `compile_line` constructor
# ```
#        ╔══════════════════════════════════╗
#        ║ EdgeModel (compiled)             ║
#    src ║ ┌──────────────────────────────┐ ║ dst
# vertex ║ │MTKLine                       │ ║ vertex
#    u ───→│┌───────┐ ┌────────┐ ┌───────┐│←─── u
#        ║ ││LineEnd├o┤ PiLine ├o┤LineEnd││ ║
#    i ←───│└───────┘ └────────┘ └───────┘│───→ i
#        ║ └──────────────────────────────┘ ║
#        ╚══════════════════════════════════╝
# ```
# (We used the `PiLine_fault` model since we plan on simulating short circuits later.)
#
# A line has two voltage bases, `src₊Vbase` and `dst₊Vbase`, one per `LineEnd` — a transformer
# between the 16.5 kV generator buses and the 345 kV grid has genuinely different bases on its
# two ends. We do not set them here: each `LineEnd` declares that it inherits its `Vbase` from
# the bus it is attached to, and that inheritance is resolved during initialization, once the
# line knows its neighbours. Part II shows the result.

# Each transmission element is created by:
# 1. Compiling a line model with its source and destination buses
# 2. Setting electrical parameters (resistance, reactance, susceptance) from CSV data

branches = []
for row in eachrow(branch_df)
    # Create line instance with topology
    line = compile_line(MTKLine(PiLine_fault(;name=:piline)); src=row.src_bus, dst=row.dst_bus)

    # Apply electrical parameters from CSV data
    for col_name in names(branch_df)
        if col_name ∉ ["src_bus", "dst_bus", "transformer"]
            set_default!(line, Regex(col_name*"\$"), row[col_name])
        end
    end

    push!(branches, line)
end

# ## Network Assembly
#
# Finally, we combine all buses and transmission lines into a complete network model.
# This creates the IEEE 39-bus test system ready for initialization and simulation.

nw = Network(busses, branches)

# All models are built, so we restore the global bases to their defaults. Calling the setters
# without an argument does that, and it keeps this page's 60 Hz from leaking into whatever is
# constructed next in the same session.

set_Sbase!()
set_fbase!()

# The network `nw` now contains the complete IEEE 39-bus model structure.
# In Part 2 of this tutorial series, we'll initialize this network by solving
# the power flow and setting up the dynamic initial conditions.

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
