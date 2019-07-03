using CSV: getfield, CSV
using DataFrames: names!, getfield
using LightGraphs: SimpleGraph, add_edge!, add_vertices!, Edge
using Logging: @warn


function read_network_from_csv(bus_file, line_file)
    graph = SimpleGraph()
    nodes = _read_nodes_from_csv(bus_file)
    num_nodes = length(nodes)
    add_vertices!(graph, num_nodes)
    lines = _read_lines_from_csv(line_file, graph, num_nodes)
    PowerGrid(graph, nodes, lines)
end

function _read_nodes_from_csv(filename)
    busses_df = CSV.read(filename)[[2,5,6,7,8,9,10]]
    df_names = [:type, :P_gen, :Q_gen, :P_load, :Q_load, :intertia, :damping]
    names!(getfield(busses_df, :colindex), df_names)

    node_list = []
    for bus_index = 1:size(busses_df)[1]
        bus_data = busses_df[bus_index,:]
        if busses_df[bus_index,:type] == "S"
            append!(node_list, [SlackAlgebraic(U=1)])
        elseif busses_df[bus_index,:type] == "G"
            append!(node_list, [SwingEqLVS(
                H=busses_df[bus_index,:intertia]  ,
                P=(busses_df[bus_index,:P_gen] - busses_df[bus_index,:P_load]),
                D=busses_df[bus_index,:damping],
                Ω=50,
                Γ= 2,
                V=1
            )])
        elseif busses_df[bus_index,:type] == "L"
            append!(node_list, [PQAlgebraic(
                S= -busses_df[bus_index,:P_load] - im*busses_df[bus_index,:Q_load]
            )])
        end
    end
    node_list
end

function _read_lines_from_csv(filename, graph, num_nodes)
    lines_df = CSV.read(filename)
    df_names = [:from, :to, :R, :X, :charging, :tap_ratio]
    names!(getfield(lines_df, :colindex), df_names)

    line_list = []
    for line_index = 1:size(lines_df)[1]
        from = lines_df[line_index,:from]
        to = lines_df[line_index,:to]
        if (from > num_nodes) || (to > num_nodes)
            @warn("Skipping line $line_index from $from to $(to)!")
            continue
        end
        admittance = 1/(lines_df[line_index,:R] + im*lines_df[line_index,:X])
        line = StaticLine(from=from, to=to, Y=admittance)
        push!(line_list, line)
        add_edge!(graph, from, to);
    end
    line_list
end
