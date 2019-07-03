using LightGraphs: SimpleGraph, add_edge!, add_vertices!, Edge
using JSON: parsefile
using Logging: @warn

function read_network_from_json(file)
    json = parsefile(file; dicttype=Dict, inttype=Int64, use_mmap=true)
    nodes = get(json, "nodes", []) |> convert_nodes
    lines = get(json, "lines", []) |> convert_lines
    PowerGrid(nodes, lines)
end

function convert_nodes(nodes)
    map(convert_node, nodes)
end

function convert_node(node)
    type = get(node, "type", nothing)
    params = get(node, "params", [])
    sym_params = Dict(Symbol(k) => _map_complex(v) for (k, v) in params)
    if type == "SwingEq"
        SwingEq(;sym_params...)
    elseif type == "SlackAlgebraic"
        SlackAlgebraic(;sym_params...)
    else
        throw(ArgumentError("Invalid type: $type"))
    end
end

function convert_lines(lines)
    map(convert_line, lines)
end

function convert_line(line)
    type = get(line, "type", nothing)
    params = get(line, "params", [])
    sym_params = Dict(Symbol(k) => _map_complex(v) for (k, v) in params)
    if type == "StaticLine"
        StaticLine(;sym_params...)
    elseif type == "PiModelLine"
        PiModelLine(;sym_params...)
    else
        throw(ArgumentError("Invalid type: $type"))
    end
end

function _map_complex(v)
    if isa(v, Dict) && haskey(v, "r") && haskey(v, "i")
        Complex(get(v, "r", nothing), get(v, "i", nothing))
    else
        v
    end
end

export read_network_from_json
