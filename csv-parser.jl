using CSV
using DataFrames
using LightGraphs
using PowerDynBase
using NetworkDynamics

function read_network_from_csv(bus_file, line_file)
    nodes = _read_nodes_from_csv(bus_file)
    network_graph, lines = _read_lines_from_csv(line_file, length(nodes))
    network_dynamics(map(construct_node_dynamics, nodes), lines, network_graph)
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
    return node_list
end

function _read_lines_from_csv(filename, num_nodes)
    lines_df = CSV.read(filename)
    df_names = [:from, :to, :R, :X, :charging, :tap_ratio]
    names!(getfield(lines_df, :colindex), df_names)

    g = SimpleGraph(num_nodes);
    edge_list = Array{StaticEdge, 1}()
    for line_index = 1:size(lines_df)[1]
        from = lines_df[line_index,:from]
        to = lines_df[line_index,:to]
        if (from > num_nodes) || (to > num_nodes)
            warn("Skipping line $line_index from $from to $(to)!")
            continue
        end
        admittance = 1/(lines_df[line_index,:R] + im*lines_df[line_index,:X])
        edge = construct_edge(StaticLine(Y=admittance))
        append!(edge_list, [edge])
        add_edge!(g, from, to);
    end
    g, edge_list
end
