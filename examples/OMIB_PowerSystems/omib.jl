using PowerDynamics
using PowerSystems
using ModelingToolkit

raw_file = joinpath(@__DIR__, "OMIB.raw")
dyr_file = joinpath(@__DIR__, "OMIB.dyr")
system = System(raw_file, dyr_file, runchecks = true)

collect(get_components(Generator, system))[1]
collect(get_components(DynamicGenerator, system))[1]

collect(get_components(Branch, system))[1]
collect(get_components(Branch, system))[2]

collect(get_components(ElectricLoad, system))
collect(get_components(Storage, system))
collect(get_components(LoadZone, system))[1]


dat = system.data
comp = dat.components

collect(get_components(Bus, system))[1]
collect(get_components(Bus, system))[2]

collect(get_components(Arc, system))[1]

dyngen = collect(get_components(DynamicGenerator, system))[1]
metagen = MetaGenerator(dyngen; verbose=false);
busnode = BusNode(metagen; name=:bus1)

# see information
symbolsof(busnode)
busnode.block
busnode.parameter_names .=> busnode.parameters
busnode.block.system.eqs
