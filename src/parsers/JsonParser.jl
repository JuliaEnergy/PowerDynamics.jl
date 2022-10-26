using Graphs: SimpleGraph, add_edge!, add_vertices!, Edge
using JSON: parsefile, print
using Logging: @warn


abstract type Json <: Format end

"""
    read_powergrid(file, Json)

Parses an existing model in JSON format into a [`PowerGrid`](@ref)
"""
function read_powergrid(file, ::Type{Json})
    # only use mmap on non-windows systems
    # https://github.com/JuliaIO/JSON.jl/issues/112
    mmap = !Sys.iswindows()
    json = parsefile(file; dicttype=Dict, inttype=Int64, use_mmap=mmap)
    nodes = get(json, "nodes", []) |> convert_nodes |> component_format
    lines = get(json, "lines", []) |> convert_lines |> component_format
    PowerGrid(nodes, lines)
end

"""
    write_powergrid(powergrid, file, Json)

Writes a [`PowerGrid`](@ref) model into a file as JSON.
"""
function write_powergrid(pg::PowerGrid, file, ::Type{Json})
    json_nodes = write_types(pg.nodes)
    json_lines = write_types(pg.lines)
    json_dict = Dict("version" => "1","nodes" => json_nodes, "lines" => json_lines)
    open(file,"w") do f
        print(f, json_dict, 2)
    end
end

function component_format(components)
    components
end

function component_format(components::Array{<:Tuple})
    OrderedDict(components)
end

function convert_nodes(nodes)
    map(convert_node, nodes)
end

function convert_node(node)
    name = get(node, "name", nothing)
    type = get(node, "type", nothing)
    params = get(node, "params", [])
    sym_params = Dict(Symbol(k) => _map_complex(v) for (k, v) in params)
    if type == "SwingEq"
        node = SwingEq(;sym_params...)
    elseif type == "SwingEqLVS"
        node = SwingEqLVS(;sym_params...)
    elseif type == "FourthOrderEq"
        node = FourthOrderEq(;sym_params...)
    elseif type == "FourthOrderEqGovernorExciterAVR"
        node = FourthOrderEqGovernorExciterAVR(;sym_params...)
    elseif type == "SlackAlgebraic"
        node = SlackAlgebraic(;sym_params...)
    elseif type == "PQAlgebraic"
        node = PQAlgebraic(;sym_params...)
    elseif type == "PVAlgebraic"
        node = PVAlgebraic(;sym_params...)
    elseif type == "PVAlgebraic"
        node = PVAlgebraic(;sym_params...)
    elseif type == "VSIMinimal"
        node = VSIMinimal(;sym_params...)
    elseif type == "VSIVoltagePT1"
        node = VSIVoltagePT1(;sym_params...)
    elseif type == "CSIMinimal"
        node = CSIMinimal(;sym_params...)
    elseif type == "ExponentialRecoveryLoad"
        node = ExponentialRecoveryLoad(;sym_params...)
    elseif type == "VoltageDependentLoad"
        node = VoltageDependentLoad(;sym_params...)
    elseif type == "NormalForm"
        node = NormalForm(;sym_params...)
    else
        throw(ArgumentError("Invalid type: $type"))
    end
    if typeof(name) == Nothing
        node
    else
        (name,node)
    end
end

function convert_lines(lines)
    map(convert_line, lines)
end

function convert_line(line)
    name = get(line, "name", nothing)
    type = get(line, "type", nothing)
    params = get(line, "params", [])
    sym_params = Dict(Symbol(k) => _map_complex(v) for (k, v) in params)
    if type == "StaticLine"
        line = StaticLine(;sym_params...)
    elseif type == "PiModelLine"
        line = PiModelLine(;sym_params...)
    elseif type == "Transformer"
        line = Transformer(;sym_params...)
    else
        throw(ArgumentError("Invalid type: $type"))
    end
    if typeof(name) == Nothing
        line
    else
        (name,line)
    end
end

function _map_complex(v)
    if isa(v, Dict) && haskey(v, "re") && haskey(v, "im")
        Complex(get(v, "re", nothing), get(v, "im", nothing))
    else
        v
    end
end

function write_type(component)
    params = typedict(component)
    dict = Dict("type" => string(component |> typeof |> nameof), "params" => params)
end

function write_type(component,name)
    params = typedict(component)
    dict = Dict("type" => string(component |> typeof |> nameof), "params" => params, "name" => string(name))
end

function write_types(components)
    map(write_type, components)
end

function write_types(components::OrderedDict)
    map(write_type, values(components), keys(components))
end

typedict(x) = Dict(fn => getfield(x, fn) for fn ∈ fieldnames(typeof(x)))

export read_powergrid, write_powergrid, Json
