using Ipopt, PowerModels

function power_flow(power_grid :: PowerGrid)
    global ang_min, ang_max
    ang_min = deg2rad(60)
    ang_max = deg2rad(-60)

    function make_generator(dict :: Dict{String, Any}, key_b :: Int, node)
        key = length(dict["gen"])+1
        (dict["gen"])[string(key)] = Dict{String, Any}()
        ((dict["gen"])[string(key)])["mBase"] = 1
        ((dict["gen"])[string(key)])["gen_bus"] = key_b
        ((dict["gen"])[string(key)])["gen_status"] = 1
        ((dict["gen"])[string(key)])["source_id"] = Any["gen", key]
        ((dict["gen"])[string(key)])["index"] = key

        if isa(node, SlackAlgebraic)
            ((dict["gen"])[string(key)])["vg"] = node.U
            # ((dict["gen"])[string(key)])["pg"] = no
            ((dict["gen"])[string(key)])["pmin"] = 0
            # ((dict["gen"])[string(key)])["pmax"] = 1.1 * node.P
        else
            ((dict["gen"])[string(key)])["pg"] = node.P
            ((dict["gen"])[string(key)])["pmin"] = 0.9 * node.P
            ((dict["gen"])[string(key)])["pmax"] = 1.1 * node.P
            ((dict["gen"])[string(key)])["vg"] = node.V
        end

        ((dict["gen"])[string(key)])["qg"] = 0
        ((dict["gen"])[string(key)])["qmin"] = -1.5
        ((dict["gen"])[string(key)])["qmax"] = 1.5
    end

    function make_bus_ac(data :: Dict{String, Any}, node)
        key_e = length(data["bus"])+1
        (data["bus"])[string(key_e)] = Dict{String, Any}()
        dict = (data["bus"])[string(key_e)]
        dict["source_id"] = Any["bus", key_e]
        dict["index"] = key_e
        #dict["bus_i"] = key_e
        #dict["zone"] = 1
        #dict["area"] = 1
        dict["vmin"] = 0.9
        dict["vmax"] = 1.1
        dict["vm"] = 1
        dict["va"] = 0
        #dict["base_kv"] = 1

        if isa(node, PQAlgebraic)
            dict["bus_type"] = 1

            key_l = length(data["load"]) + 1
            (data["load"])[string(key_l)] = Dict{String, Any}()
            dict_l = (data["load"])[string(key_l)]
            dict_l["source_id"] = dict["source_id"]
            dict_l["load_bus"] = key_e
            dict_l["status"] = 1
            dict_l["pd"] = node.P
            dict_l["qd"] = node.Q
            dict_l["index"] = key_l
        elseif (isa(node, PVAlgebraic) || isa(node, SwingEqLVS))
            dict["bus_type"] = 2
            dict["vm"] = abs(node.V)  # assumed p.u.
            dict["vmin"] = 0.9 * dict["vm"]
            dict["vmax"] = 1.1 * dict["vm"]
            dict["va"] = angle(node.V)

            isa(node, SwingEqLVS) ? make_generator(data, key_e, node) : nothing
        elseif isa(node, SlackAlgebraic)
            dict["bus_type"] = 3
            dict["vm"] = abs(node.U)  # assumed p.u.
            dict["vmin"] = 0.9 * dict["vm"]
            dict["vmax"] = 1.1 * dict["vm"]
            dict["va"] = angle(node.U)
            make_generator(data, key_e, node)
        else
            dict["bus_type"] = 4
        end
    end

    function make_branch_ac(data :: Dict{String, Any}, line)
        key_e = length(data["branch"])+1
        (data["branch"])[string(key_e)] = Dict{String, Any}()
        data = (data["branch"])[string(key_e)]

        data["f_bus"] = line.from
        data["t_bus"] = line.to
        data["source_id"] = Any["branch", key_e]
        data["index"] = key_e
        data["rate_a"] = 1
        data["rate_b"] = 1
        data["rate_c"] = 1
        data["c_rating_a"] = 1
        data["br_status"] = 1
        data["angmin"] = ang_min
        data["angmax"] = ang_max

        # default values
        data["transformer"] = false
        data["tap"] = 1
        data["shift"] = 0
        data["g_fr"] = 0
        data["b_fr"] = 0
        data["g_to"] = 0
        data["b_to"] = 0

        if isa(line, StaticLine)
            data["br_r"] = real(1/line.Y)
            data["br_x"] = imag(1/line.Y)
        elseif isa(line, RLLine)
            data["br_r"] = real(line.R)
            data["br_x"] = imag(line.L*line.ω0)
        elseif isa(line, Transformer)
            data["transformer"] = true
            data["tap"] = line.t_ratio
            data["br_r"] = real(1/line.Y)
            data["br_x"] = imag(1/line.Y)
        elseif isa(line, PiModel)
            data["transformer"] = true
            data["tap"] = line.t_km / line.t_mk
            data["g_fr"] = real(line.y_shunt_km / data["tap"]^2)
            data["b_fr"] = imag(line.y_shunt_km / data["tap"]^2)
            data["br_r"] = real(1/line.Y)
            data["br_x"] = imag(1/line.Y)
            data["g_to"] = real(line.y_shunt_mk)
            data["b_to"] = imag(line.y_shunt_mk)
        elseif isa(line, PiModelLine)
            data["g_fr"] = real(line.y_shunt_km)
            data["b_fr"] = imag(line.y_shunt_km)
            data["br_r"] = real(1/line.Y)
            data["br_x"] = imag(1/line.Y)
            data["g_to"] = real(line.y_shunt_mk)
            data["b_to"] = imag(line.y_shunt_mk)
        else
            throw(ArgumentError("Line type $(typeof(line)) does not exist."))
        end
    end

    data = Dict{String, Any}()
    data["source_type"] = "matpower"
    data["name"] = "network"
    data["source_version"] = "0.0.0"
    data["per_unit"] = true
    data["baseMVA"] = 1
    keywords = ["bus"; "shunt"; "dcline";
                "storage"; "switch"; "load"; "branch"; "gen"]
    for keyword in keywords
        data[keyword] = Dict{String, Any}()
    end

    for line in power_grid.lines
        make_branch_ac(data, line)
    end
    for node in power_grid.nodes
        make_bus_ac(data, node)
    end

    result = run_pf(data, ACPPowerModel, Ipopt.Optimizer)

    return data, result
end

export power_flow
