function PowerDynamics.PowerGrid(sys::System)
    busiter = get_components(Bus, sys)
    N = length(get_components(Bus, sys))
    busnames = get_name.(busiter)

    busdevices = Dict(busnames .=> Ref(Tuple{IOBlock, Dict}[]))
    # generate IOBlocks for all the loads
    for load in get_components(ElectricLoad, sys)
        busname = get_name(get_bus(load))
        ioload = get_io_load(load)
        push!(busdevices[busname], ioload)
    end

    # generate IOBlocks for all the dynamic injectors
    for injection in get_components(DynamicInjection, sys)
        busname = get_name(get_bus(injection))
        ioinj = get_io_injector(injection)
        push!(busdevices[busname], ioinj)
    end

    # get the ordered dict of IONodes for PD
    buses = OrderedDict{String, IONode}
    for (name, devices) in busdevices
        buses[name] = BusNOde(devices)
    end
end
