using Ipopt: Ipopt
using PowerModelsACDC: run_acdcpf
using PowerModels

#=
function _make_generator_header(dict::Dict{String,Any}, key_b::Int)
    key = length(dict["gen"]) + 1
    dict["gen"][string(key)] = Dict{String,Any}()
    dict["gen"][string(key)]["mBase"] = 1
    dict["gen"][string(key)]["gen_bus"] = key_b
    dict["gen"][string(key)]["gen_status"] = 1
    dict["gen"][string(key)]["source_id"] = Any["gen", key]
    dict["gen"][string(key)]["index"] = key

    dict["gen"][string(key)]["qg"] = 0
    dict["gen"][string(key)]["qmin"] = -1.5
    dict["gen"][string(key)]["qmax"] = 1.5
end

function make_generator!(dict::Dict{String,Any}, key_b::Int, node::SlackAlgebraic)
    _make_generator_header(dict, key_b)
    key = length(dict["gen"])
    dict["gen"][string(key)]["vg"] = node.U
end

function make_generator!(dict::Dict{String,Any}, key_b::Int, node::NormalForm)
    make_generator_header(dict, key_b)
    key = length(dict["gen"])
    dict["gen"][string(key)]["pg"] = node.P
    dict["gen"][string(key)]["pmin"] = 0.9 * node.P
    dict["gen"][string(key)]["pmax"] = 1.1 * node.P
    dict["gen"][string(key)]["vg"] = node.V
end

"""
The entries that are present for all busses.
"""
function _make_bus_ac_header(data::Dict{String,Any})
    key_e = length(data["bus"]) + 1
    data["bus"][string(key_e)] = Dict{String,Any}()
    bus_dict = data["bus"][string(key_e)]
    bus_dict["source_id"] = Any["bus", key_e]
    bus_dict["index"] = key_e
    bus_dict["vmin"] = 0.9
    bus_dict["vmax"] = 1.1
    bus_dict["vm"] = 1
    bus_dict["va"] = 0
    return bus_dict
end

"""
The default case.
"""
function make_bus_ac!(data::Dict{String,Any}, node)
    bus_dict = _make_bus_ac_header(data)
    bus_dict["bus_type"] = 4
end

function make_bus_ac!(data::Dict{String,Any}, node::PQAlgebraic)
    bus_dict = _make_bus_ac_header(data)
    bus_dict["bus_type"] = 1

    key_l = length(data["load"]) + 1
    data["load"][string(key_l)] = Dict{String,Any}()
    dict_l = (data["load"])[string(key_l)]
    dict_l["source_id"] = bus_dict["source_id"]
    dict_l["load_bus"] = bus_dict["index"]
    dict_l["status"] = 1
    dict_l["pd"] = -node.P
    dict_l["qd"] = -node.Q
    dict_l["index"] = key_l
end

function make_bus_ac!(data::Dict{String,Any}, node::PVAlgebraic)
    bus_dict = _make_bus_ac_header(data)
    bus_dict["bus_type"] = 2
    bus_dict["vm"] = abs(node.V)  # assumed p.u.
    bus_dict["vmin"] = 0.9 * bus_dict["vm"]
    bus_dict["vmax"] = 1.1 * bus_dict["vm"]
    bus_dict["va"] = angle(node.V)
end

# function make_bus_ac!(data::Dict{String,Any}, node::NormalForm)
#     bus_dict = _make_bus_ac_header(data)
#     bus_dict["bus_type"] = 2
#     bus_dict["vm"] = abs(node.V)
#     bus_dict["vmin"] = 0.9 * bus_dict["vm"]
#     bus_dict["vmax"] = 1.1 * bus_dict["vm"]
#     bus_dict["va"] = angle(node.V)
#     make_generator!(data, bus_dict["index"], node)
# end

function make_bus_ac!(data::Dict{String,Any}, node::SlackAlgebraic)
    bus_dict = _make_bus_ac_header(data)
    bus_dict["bus_type"] = 3
    bus_dict["vm"] = abs(node.U)  # assumed p.u.
    bus_dict["vmin"] = 0.9 * bus_dict["vm"]
    bus_dict["vmax"] = 1.1 * bus_dict["vm"]
    bus_dict["va"] = angle(node.U)
    make_generator!(data, bus_dict["index"], node)
end

# function make_bus_ac!(data::Dict{String,Any}, node::RLCLoad)
#     bus_dict = _make_bus_ac_header(data)
#     bus_dict["bus_type"] = 3
#     bus_dict["vm"] = 1  # assumed p.u.
#     bus_dict["vmin"] = 0.9 * bus_dict["vm"]
#     bus_dict["vmax"] = 1.1 * bus_dict["vm"]
#     bus_dict["va"] = 0
#     make_shunt!(data, bus_dict["index"], 100π, node) # atm ω is set fixed 100π
# end

function _make_branch_ac_header(data::Dict{String,Any}, dict::Dict{Any, Int}, line)
    key_e = length(data["branch"]) + 1
    data["branch"][string(key_e)] = Dict{String,Any}()
    branch_dict = data["branch"][string(key_e)]
    branch_dict["source_id"] = Any["branch", key_e]
    branch_dict["index"] = key_e
    branch_dict["rate_a"] = 1
    branch_dict["rate_b"] = 1
    branch_dict["rate_c"] = 1
    branch_dict["c_rating_a"] = 1
    branch_dict["br_status"] = 1
    branch_dict["angmin"] = ang_min
    branch_dict["angmax"] = ang_max

    # default values
    branch_dict["transformer"] = false
    branch_dict["tap"] = 1
    branch_dict["shift"] = 0
    branch_dict["g_fr"] = 0
    branch_dict["b_fr"] = 0
    branch_dict["g_to"] = 0
    branch_dict["b_to"] = 0

    # addition for dictionary use
    if isempty(dict)
        branch_dict["f_bus"] = line.from
        branch_dict["t_bus"] = line.to
    else
        branch_dict["f_bus"] = dict[line.from]
        branch_dict["t_bus"] = dict[line.to]
    end

    return branch_dict
end

"""
The default, non-specialised method.
"""
function make_branch_ac!(data::Dict{String,Any}, dict::Dict{Any, Int}, line)
    throw(ArgumentError("Line type $(typeof(line)) does not exist."))
end

function make_branch_ac!(data::Dict{String,Any}, dict::Dict{Any, Int}, line::StaticLine)
    branch_dict = _make_branch_ac_header(data, dict, line)
    branch_dict["br_r"] = real(1 / line.Y)
    branch_dict["br_x"] = imag(1 / line.Y)
end

# function make_branch_ac!(data::Dict{String,Any}, dict::Dict{Any, Int}, line::RLLine)
#     branch_dict = _make_branch_ac_header(data, dict, line)
#     branch_dict["br_r"] = real(line.R)
#     branch_dict["br_x"] = imag(line.L * line.ω0)
# end

# function make_branch_ac!(data::Dict{String,Any}, dict::Dict{Any, Int}, line::Transformer)
#     branch_dict = _make_branch_ac_header(data, dict, line)
#     branch_dict["transformer"] = true
#     branch_dict["tap"] = line.t_ratio
#     branch_dict["br_r"] = real(1 / line.y)
#     branch_dict["br_x"] = imag(1 / line.y)
# end

# function make_branch_ac!(data::Dict{String,Any}, dict::Dict{Any, Int}, line::PiModelLine)
#     branch_dict = _make_branch_ac_header(data, dict, line)
#     branch_dict["g_fr"] = real(line.y_shunt_km)
#     branch_dict["b_fr"] = imag(line.y_shunt_km)
#     branch_dict["br_r"] = real(1 / line.y)
#     branch_dict["br_x"] = imag(1 / line.y)
#     branch_dict["g_to"] = real(line.y_shunt_mk)
#     branch_dict["b_to"] = imag(line.y_shunt_mk)
# end



export power_flow
function powermodels_pf(nodes, lines)
    # TODO write a wrapper for mapping dict strings to integer lists and remember it for mapping back the results to string bus names
    data = Dict{String,Any}()
    data["source_type"] = "PowerDynamics"
    data["name"] = "network"
    data["source_version"] = "2.3.2"
    data["per_unit"] = true
    data["baseMVA"] = 1
    keywords = [
        "bus"
        "shunt"
        "dcline"
        "storage"
        "switch"
        "load"
        "branch"
        "gen"
    ]

    for keyword in keywords
        data[keyword] = Dict{String,Any}()
    end

    # addition for dictionary
    dict = Dict{Any, Int}()
    for node in nodes
        make_bus_ac!(data, node)
    end

    for line in lines
        make_branch_ac!(data, dict, line)
    end

    s = Dict("output" => Dict("branch_flows" => true), "conv_losses_mp" => true)
    result = run_acdcpf(data, ACPPowerModel, Optimizer; setting = s)

    return data, result
end
=#

const ang_min = deg2rad(60)
const ang_max = deg2rad(-60)
function powermodels_pf(nodes, lines)
    data = Dict{String,Any}()
    data["per_unit"] = true
    data["baseMVA"] = 1

    # sub-dictionaries
    data["bus"] = Dict{String,Any}()
    data["gen"] = Dict{String,Any}()
    data["load"] = Dict{String,Any}()
    data["branch"] = Dict{String,Any}()
    data["shunt"] = Dict{String,Any}() # not used
    data["storage"] = Dict{String,Any}() # not used
    data["switch"] = Dict{String,Any}() # not used
    data["dcline"] = Dict{String,Any}() # not used

    for (i, node) in enumerate(nodes)
        busd = data["bus"][string(i)] = Dict{String,Any}()
        busd["index"] = i
        busd["vmin"] = 0.9
        busd["vmax"] = 1.1
        if node isa PQAlgebraic
            loadid = length(data["load"]) + 1
            load = data["load"][string(loadid)] = Dict{String,Any}()
            load["load_bus"] = i
            load["pd"] = node.P
            load["qd"] = node.Q
            load["index"] = loadid
        elseif node isa PVAlgebraic
            busd["bus_type"] = 2
            busd["vm"] = node.V
            genid = length(data["gen"]) + 1
            gen = data["gen"][string(genid)] = Dict{String,Any}()
            gen["gen_bus"] = i
            gen["gen_status"] = 1
            gen["pg"] = node.P
            gen["index"] = genid
        elseif node isa SlackAlgebraic
            busd["bus_type"] = 3
            busd["vm"] = node.U
            busd["va"] = 0
            genid = length(data["gen"]) + 1
            gen = data["gen"][string(genid)] = Dict{String,Any}()
            gen["gen_status"] = 1
            gen["gen_bus"] = i
            gen["index"] = genid
        else
            error("Uknown node type $(typeof(node))")
        end
    end
    for (i, line) in enumerate(lines)
        branch = data["branch"][string(i)] = Dict{String,Any}()
        branch["f_bus"] = line.from
        branch["t_bus"] = line.to
        branch["br_status"] = 1
        branch["index"] = i
        branch["angmin"] = PowerDynamicsPrototype.ang_min
        branch["angmax"] = PowerDynamicsPrototype.ang_max
        branch["tap"] = 1
        branch["shift"] = 0
        if line isa StaticLine
            branch["br_r"] = real(1 / line.Y)
            branch["br_x"] = imag(1 / line.Y)
            branch["g_fr"] = 0
            branch["b_fr"] = 0
            branch["g_to"] = 0
            branch["b_to"] = 0
        elseif line isa PiModelLine
            branch["br_r"] = real(1 / line.y)
            branch["br_x"] = imag(1 / line.y)
            branch["g_fr"] = real(line.y_shunt_km)
            branch["b_fr"] = imag(line.y_shunt_km)
            branch["g_to"] = real(line.y_shunt_mk)
            branch["b_to"] = imag(line.y_shunt_mk)
        else
            error("Uknown line type $(typeof(line))")
        end
    end

    # PowerModels.print_summary(data)
    res = PowerModels.solve_ac_pf(data, Ipopt.Optimizer)
    busses = res["solution"]["bus"]
    voltages = zeros(ComplexF64, length(nodes))
    for i in eachindex(voltages)
        va = busses[string(i)]["va"]
        vm = busses[string(i)]["vm"]
        voltages[i] = vm * exp(im * va)
    end
    voltages
end

export power_flow
function power_flow(nodes, lines)
    voltages = PowerDynamicsPrototype.powermodels_pf(nodes, lines)

    nw = Network(nodes, lines)
    s = NWState(nw)
    s.v[:, :u_r] .= real.(voltages)
    s.v[:, :u_i] .= imag.(voltages)

    du = [Inf for _ in 1:dim(nw)]
    nw(du, uflat(s), pflat(s), s.t)
    if maximum(abs.(du)) > 1e-6
        @warn "Power Models did not find a suitable power flow solution."
    end
    s
end
