export convert_initial_conditions, get_grid_and_x0

function PowerDynamics.PowerGrid(sys::System; verbose=false)
    # warnings for currently unssuported components
    isempty(get_components(StaticInjectionSubsystem, sys)) || @warn "Ignore devices of type StaticInjectionSubsystem"
    isempty(get_components(Storage, sys)) || @warn "Ignore devices of type StaticInjectionSubsystem"
    isempty(get_components(RegulationDevice, sys)) || @warn "Ignore devices of type StaticInjectionSubsystem"

    ####
    #### buses
    ####
    busnames = sort!(get_name.(get_components(Bus, sys)))
    N = length(busnames)

    busdevices = OrderedDict{String, Vector{BlockPara}}()
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
    buses = OrderedDict{String, AbstractNode}()
    for (name, devices) in busdevices
        verbose && println(" - $(name): ", (d->d.block.name).(devices))
        buses[name] = BusNode(devices...; name=Symbol(name))
    end

    # the Source is an infinite/slack bus
    for source in get_components(Source, sys)
        source.dynamic_injector !== nothing && error("There is a dynamic injector attached to the source $(get_name(source))!")
        bus = get_name(get_bus(source))
        !isempty(busdevices[bus]) && @warn("There are other devices attached to the slack $bus: ", busdevices[bus])

        U = get_internal_voltage(source) * exp(im* get_internal_angle(source))

        verbose && println("Override bus $bus with SlackAlgebraic(U=$U)")
        buses[bus] = SlackAlgebraic(; U)
    end

    ####
    #### Branches
    ####
    branches = OrderedDict{String, PowerDynamics.AbstractLine}()
    verbose && println("Branches ")
    for arc in get_components(Arc, sys)
        lines = get_components(Branch, sys, x->x.arc==arc) |> collect
        lines = [l for l in lines]
        verbose && println(" - $(get_name(arc)): ", (l->repr(typeof(l))*": "*get_name(l)).(lines))
        pdline = get_pd_branch(arc, lines)
        branches[get_name(arc)] = pdline
    end

    return PowerGrid(buses, branches)
end

function convert_initial_conditions(pg::PowerGrid, sys::System, ic)
    nd_syms = syms_containing(rhs(pg), "")
    x0 = zeros(length(nd_syms))
    busnames = pg.nodes.keys

    buses = get_components(Bus, sys)
    ps_busnrs = Dict(get_name.(buses) .=> get_number.(buses))

    for (state_i, sym) in enumerate(nd_syms)
        str = String(sym)
        # in calse of voltage
        m = match(r"^u_([r,i])_([0-9]+)$", str)
        if m !== nothing
            nd_busnr = parse(Int, m.captures[2])
            busname = busnames[nd_busnr]

            if m.captures[1] == "r"
                v_key = "V_R"
            elseif m.captures[1] == "i"
                v_key = "V_I"
            else
                error("Can't handle state symbol $str")
            end

            ps_busnr = ps_busnrs[busname]
            x0[state_i] = ic[v_key][ps_busnr]
            continue
        end

        m = match(r"^(.+)₊(.+)_([0-9]+)$", str)
        if m !== nothing
            blockname = m.captures[1]
            varname = m.captures[2]

            if haskey(ic, blockname)
                x0[state_i] = ic[blockname][Symbol(varname)]
                continue
            elseif !isempty(get_components(PowerLoad, sys, x->get_name(x)==blockname))
                loads = get_components(PowerLoad, sys, x->get_name(x)==blockname)
                length(loads)==1 || error("There is more than on load $(get_name.(loads))")
                load = first(loads)
                # power draw from load
                S = load.active_power + im * load.reactive_power

                # get voltage of node
                nd_busnr = parse(Int, m.captures[3])
                busname = busnames[nd_busnr]
                ps_busnr = ps_busnrs[busname]
                U = ic["V_R"][ps_busnr] + im*ic["V_I"][ps_busnr]

                # calculate the current injection (negative current draw)
                Iconj = - S/U

                if varname == "i_r"
                    x0[state_i] = real(Iconj)
                elseif varname == "i_i"
                    x0[state_i] = -imag(Iconj)
                end

                continue
            end
        end

        @warn "Can't handle state symbol $str"
    end

    return x0
end

function get_grid_and_x0(sim::Simulation)
    ic = get_initial_conditions(sim)
    pg = PowerGrid(sim.sys)
    x0 = convert_initial_conditions(pg, sim.sys, ic)
    return (pg, x0)
end
