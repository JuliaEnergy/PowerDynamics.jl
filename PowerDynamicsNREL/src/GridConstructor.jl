function PowerDynamics.PowerGrid(sys::System; verbose=false)
    isempty(get_components(Source, sys)) || @warn "Ignore devices of type Source"
    isempty(get_components(StaticInjectionSubsystem, sys)) || @warn "Ignore devices of type StaticInjectionSubsystem"
    isempty(get_components(Storage, sys)) || @warn "Ignore devices of type StaticInjectionSubsystem"
    isempty(get_components(RegulationDevice, sys)) || @warn "Ignore devices of type StaticInjectionSubsystem"

    busiter = get_components(Bus, sys)
    N = length(get_components(Bus, sys))
    busnames = get_name.(busiter)

    busdevices = Dict{String, Vector{BlockPara}}()
    for name in busnames
        busdevices[name] = BlockPara[]
    end

    # generate IOBlocks for all the loads
    verbose && println("Loads: ")
    for load in get_components(ElectricLoad, sys)
        busname = get_name(get_bus(load))
        io_load = get_io_load(load)
        push!(busdevices[busname], io_load)
        verbose && println(" - added load $(get_name(load)) to $busname")
    end

    # generate IOBlocks for all the Generators
    verbose && println("Generators: ")
    for static_gen in get_components(Generator, sys)
        dyn_gen = static_gen.dynamic_injector
        dyn_gen === nothing && error("Whoops $(static_gen.name) has no dynamic model.")
        busname = get_name(get_bus(static_gen))
        io_inj = get_io_injection(dyn_gen)
        push!(busdevices[busname], io_inj)
        verbose && println(" - added dyn gen $(get_name(dyn_gen)) to $busname")
    end

    # get the ordered dict of IONodes for PD
    verbose && println("Busses: ")
    buses = OrderedDict{String, IONode}()
    for (name, devices) in busdevices
        verbose && println(" - $(name): ", (d->d.block.name).(devices))
        buses[name] = BusNode(devices...; name=Symbol(name))
    end
    return buses
end
