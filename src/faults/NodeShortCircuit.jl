"""
```Julia
    NodeShortCircuit(;node,Y,tspan_fault)
```
# Keyword Arguments
- `node`: number or name of the node
- `Y`: admittance of the short circuit
- `tspan_fault`: short circuit timespan

# Optional Keyword Arguments 
- `var`: symbol of the shunt admittance, default: `:Y_n`
"""
Base.@kwdef struct NodeShortCircuit <:AbstractPerturbation
    node
    Y
    tspan_fault
    var = :Y_n
end

function (nsc::NodeShortCircuit)(powergrid)
    typestable_node_field_update(powergrid, nsc.node, nsc.var, nsc.Y)
end

export NodeShortCircuit
