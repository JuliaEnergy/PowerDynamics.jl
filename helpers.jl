import PowerDynamics.dimension
using Setfield

variable_index(nodes, n::Integer, s::Symbol) = begin
    first_idx = findfirst(ns -> ns == s, symbolsof(nodes[n]))
    if first_idx == nothing
        throw(StateError("Variable: $s not defined for node: $n"))
    else
        startindex(nodes, n) + first_idx
    end
end
@views startindex(nodes, n) = begin
    if n == 1
        0
    else
        sum(map(node -> dimension(node), nodes[1:n-1]))
    end
end

function determine_rocof_nadir(powergrid,sol,final_time)
    ω_indices = findall(n -> isa(n, SwingEqLVS), powergrid.nodes)
    rocof_max = zeros(length(ω_indices))
    nadir = zeros(length(ω_indices))
    j=0
    for i in ω_indices
        j+=1
        rocof= sol(range(0.,stop=final_time,length=10000),Val{1},idxs=variable_index(powergrid.nodes, i, :ω)).u
        rocof_max[j]=maximum(abs.(rocof))
        omega= sol(range(0.,stop=final_time,length=10000),Val{0},idxs=variable_index(powergrid.nodes, i, :ω)).u
        nadir[j]=maximum(abs.(omega))
    end
    (nadir,rocof_max)
end;

function change_inertia(scaling_factor,nodes,powergrid)
    nodes_copy = copy(powergrid.nodes)
    node_to_change = nodes_copy[nodes]
    nodes_changed= []
    for node in node_to_change
        node = @set node.H *= scaling_factor
        append!(nodes_changed,[node])
    end
    nodes_copy[nodes] = nodes_changed
    PowerGrid(nodes_copy, powergrid.lines)
end
