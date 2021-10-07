# function flow_plot(op,t_idx)

#     g = op.grid.graph
#     dg = SimpleDiGraph(g)
#     arrow_size = []
#     rem_list = []

#     linename = [(e,string("branch",i)) for (i,e) in enumerate(edges(g))] |> Dict
#     lineset = linename |> keys |> Set
#     lines = op.grid.lines 

#     for i in 1:ne(dg)

#         e = collect(dg |> edges)[i]
#         src_bus = string("bus",e.src)
#         dst_bus = string("bus",e.dst)
#         u_src = op[src_bus,:u]
#         u_dst = op[dst_bus,:u]

#         if e in lineset
#             line = lines[linename[e]]
#         else
#             e_inv = Edge(e.dst,e.src)
#             line = lines[linename[e_inv]]
#         end

#         if typeof(line) == StaticLine
#             y = line.Y
#         else
#             y = line.y
#         end

#         i = y*(u_src - u_dst)
#         p = u_src*conj(i) |> real
#         b = y |> imag |> abs
#         p = b*(angle(u_src)-angle(u_dst))

#         if p < 0
#             push!(rem_list,e)
#         else
#             push!(arrow_size,10*p)
#         end
#     end
    
#     map(e->rem_edge!(dg,e),rem_list)

#     loads = [key for (key,value) in op.grid.nodes if typeof(value) == PQAlgebraic]

#     node_color = []
#     node_shape = []
#     cm = ColorSchemes.copper
#     colorval = get(cm,t_idx / maximum(t_idx))
#     colordict = Dict(loads .=> colorval)

#     for (bus,val) in buses

#         if typeof(val) == PQAlgebraic
#             push!(node_color,colordict[bus])
#             push!(node_shape,:circle)
#         elseif typeof(val) == FourthOrderEq
#             push!(node_color,:black)
#             push!(node_shape,'□')
#         else
#             push!(node_color,:black)
#             push!(node_shape,'□')
#         end
#     end

#     nlabels = repr.(1:nv(g))

#     return dg, nlabels, arrow_size, node_color, node_shape

# end

function flow_plot(op,t_idx)

    g = op.grid.graph
    nodes = op.grid.nodes
    dg = SimpleDiGraph(g)
    arrow_size = []
    rem_list = []

    linenumber = [(e,i) for (i,e) in enumerate(edges(g))] |> Dict
    lineset = edges(g) |> collect |> Set
    dg_edges = collect(dg |> edges)
    lines = op.grid.lines 

    for i in 1:ne(dg)

        e = dg_edges[i]
        src_bus = e.src
        dst_bus = e.dst
        u_src = op[src_bus,:u]
        u_dst = op[dst_bus,:u]

        if e in lineset
            line = lines[linenumber[e]]
        else
            e_inv = Edge(e.dst,e.src)
            line = lines[linenumber[e_inv]]
        end

        if typeof(line) == StaticLine
            y = line.Y
        else
            y = line.y
        end

        i = y*(u_src - u_dst)
        #p = u_src*conj(i) |> real
        b = y |> imag |> abs
        p = b*sin(angle(u_src)-angle(u_dst))

        if p < 0
            push!(rem_list,e)
        else
            push!(arrow_size,10*p)
        end
    end
    
    map(e->rem_edge!(dg,e),rem_list)

    loads = [idx for (idx,val) in enumerate(pg.nodes) if typeof(val) == PQAlgebraic];

    node_color = []
    node_shape = []
    cm = ColorSchemes.copper
    ti = (t_idx .- minimum(t_idx)); ti /= maximum(ti);
    colorval = get(cm,ti)
    colordict = Dict(loads .=> colorval)

    for (idx,bus) in enumerate(nodes)

        if typeof(bus) == PQAlgebraic
            push!(node_color,colordict[idx])
            push!(node_shape,:circle)
        elseif typeof(bus) == FourthOrderEq
            push!(node_color,:black)
            push!(node_shape,'□')
        else
            push!(node_color,:black)
            push!(node_shape,'□')
        end
    end

    nlabels = repr.(1:nv(g))

    return dg, nlabels, arrow_size, node_color, node_shape

end