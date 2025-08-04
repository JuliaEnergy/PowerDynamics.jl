#=
IEEE39 Part I: Model Creation
=========================

=#

using PowerDynamics
using PowerDynamics.Library
using ModelingToolkit
using NetworkDynamics
using DataFrames
using CSV

# Define data directory
DATA_DIR = joinpath(pkgdir(PowerDynamics), "docs", "examples", "ieee39data")

# Load CSV files
branch_df = CSV.read(joinpath(DATA_DIR, "branch.csv"), DataFrame)
bus_df = CSV.read(joinpath(DATA_DIR, "bus.csv"), DataFrame)
load_df = CSV.read(joinpath(DATA_DIR, "load.csv"), DataFrame)
machine_df = CSV.read(joinpath(DATA_DIR, "machine.csv"), DataFrame)
avr_df = CSV.read(joinpath(DATA_DIR, "avr.csv"), DataFrame)
gov_df = CSV.read(joinpath(DATA_DIR, "gov.csv"), DataFrame)

# Constants from original validation file
BASE_MVA = 100.0
BASE_FREQ = 60.0
ω_b = BASE_FREQ * 2π

#=
Bus Template Models
===================
Create template Bus models for each category to avoid recompilation.
Later we'll copy these templates and update parameters for each specific bus.

## Components
### ZIPLoad Injector Model**
=#
load = ZIPLoad(;name=:ZIPLoad)

#=
### Uncontrolled Machine Injector Model
=#
uncontrolled_machine = SauerPaiMachine(;
    τ_m_input=false,
    vf_input=false,
    name=:machine,
)

#=
### Controlled Machine Injector Model
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

#=
## Bus Templates

### Junction Bus Template
=#
@named junction_bus_template = Bus(MTKBus())
strip_defaults!(junction_bus_template) # make sure we set everything manually

#=
### Load Bus Template
=#
@named load_bus_template = Bus(MTKBus(load))
strip_defaults!(load_bus_template) # make sure we set everything manually

#=
### Controlled Machine Bus Template
=#
@named ctrld_machine_bus_template = Bus(
    MTKBus(controlled_machine);
)
strip_defaults!(ctrld_machine_bus_template) # make sure we set everything manually
set_default!(ctrld_machine_bus_template, r"S_b$", BASE_MVA) # set global values
set_default!(ctrld_machine_bus_template, r"ω_b$", ω_b)


#=
### Controlled Machine + Load Bus Template
=#
@named ctrld_machine_load_bus_template = Bus(
    MTKBus(controlled_machine, load);
)
strip_defaults!(ctrld_machine_load_bus_template) # make sure we set everything manually
set_default!(ctrld_machine_load_bus_template, r"S_b$", BASE_MVA)
set_default!(ctrld_machine_load_bus_template, r"ω_b$", ω_b)

#=
### Uncontrolled Machine + Load Bus Template
=#
@named unctrld_machine_load_bus_template = Bus(
    MTKBus(uncontrolled_machine, load);
)
strip_defaults!(unctrld_machine_load_bus_template) # make sure we set everything manually
set_default!(unctrld_machine_load_bus_template, r"S_b$", BASE_MVA)
set_default!(unctrld_machine_load_bus_template, r"ω_b$", ω_b)

#=
Bus Categorization and Creation
===============================
=#

# Helper function to apply CSV parameters to bus components
function apply_csv_params!(bus, table, bus_index)
    row_idx = findfirst(table.bus .== bus_index)

    # Apply all parameters except "bus" column
    row = table[row_idx, :]
    for col_name in names(table)
        if col_name != "bus"
            set_default!(bus, Regex(col_name*"\$"), row[col_name])
        end
    end
end

# Create buses using templates
busses = []
for row in eachrow(bus_df)
    i = row.bus
    println("Bus $i: $(row.category)")
    
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
    
    # Apply dynamic component parameters from CSV data
    row.has_load && apply_csv_params!(bus, load_df, i)
    row.has_gen && apply_csv_params!(bus, machine_df, i)
    row.has_avr && apply_csv_params!(bus, avr_df, i)
    row.has_gov && apply_csv_params!(bus, gov_df, i)
    
    # Apply powerflow model
    pf_model = if row.bus_type == "PQ"
        pfPQ(P=row.P, Q=row.Q)
    elseif row.bus_type == "PV"
        pfPV(P=row.P, V=row.V)
    elseif row.bus_type == "Slack"
        pfSlack(V=row.V, δ=0)
    end
    set_metadata!(bus, :pfmodel, pf_model)

    push!(busses, bus)
end
break

#=
Branch Template and Creation
============================
=#

# Single branch template for both regular lines and transformers
@named piline_template = Line(MTKLine(PiLine(;name=:piline)))

# Create branches using template
branches = []
for row in eachrow(branch_df)
    println("Branch: $(row.src_bus) -> $(row.dst_bus)")
    
    # Create line with src/dst
    line = Line(piline_template; src=row.src_bus, dst=row.dst_bus)
    
    # Apply parameters directly in loop
    for col_name in names(branch_df)
        if col_name ∉ ["src_bus", "dst_bus", "transformer"]
            set_default!(line, col_name, row[col_name])
        end
    end
    
    push!(branches, line)
end

nw = Network(busses, branches)
