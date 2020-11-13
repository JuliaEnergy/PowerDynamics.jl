using LightGraphs: SimpleGraph, add_edge!, add_vertices!, Edge
using JSON: parsefile, print
using Logging: @warn


abstract type Json <: Format end

"""
    read_powergrid(file, Json)

Parses an existing model in JSON format into a [`PowerGrid`](@ref)
"""
function read_powergrid(file, ::Type{Json})
    json = parsefile(file; dicttype=Dict, inttype=Int64, use_mmap=true)
    nodes = get(json, "nodes", []) |> convert_nodes
    lines = get(json, "lines", []) |> convert_lines
    PowerGrid(nodes, lines)
end

"""
    write_powergrid(powergrid, file, Json)

Writes a [`PowerGrid`](@ref) model into a file as JSON.
"""
function write_powergrid(pg::PowerGrid, file, ::Type{Json})
    json_nodes = map(write_type, values(pg.nodes))
    json_lines = map(write_type, values(pg.lines))
    json_dict = Dict("version" => "1","nodes" => json_nodes, "lines" => json_lines)
    open(file,"w") do f
        print(f, json_dict, 2)
    end
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
    elseif type == "SwingEqLVS"
        SwingEqLVS(;sym_params...)
    elseif type == "FourthOrderEq"
        FourthOrderEq(;sym_params...)
    elseif type == "FourthOrderEqGovernorExciterAVR"
        FourthOrderEqGovernorExciterAVR(;sym_params...)
    elseif type == "SlackAlgebraic"
        SlackAlgebraic(;sym_params...)
    elseif type == "PQAlgebraic"
        PQAlgebraic(;sym_params...)
    elseif type == "PVAlgebraic"
        PVAlgebraic(;sym_params...)
    elseif type == "PVAlgebraic"
        PVAlgebraic(;sym_params...)
    elseif type == "VSIMinimal"
        VSIMinimal(;sym_params...)
    elseif type == "VSIVoltagePT1"
        VSIVoltagePT1(;sym_params...)
    elseif type == "CSIMinimal"
        CSIMinimal(;sym_params...)
    elseif type == "ExponentialRecoveryLoad"
        ExponentialRecoveryLoad(;sym_params...)
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
    if isa(v, Dict) && haskey(v, "re") && haskey(v, "im")
        Complex(get(v, "re", nothing), get(v, "im", nothing))
    else
        v
    end
end

function write_type(t)
    params = typedict(t)
    dict = Dict("type" => string(t |> typeof |> nameof), "params" => params)
end

typedict(x) = Dict(fn => getfield(x, fn) for fn ∈ fieldnames(typeof(x)))

export read_powergrid, write_powergrid, Json
