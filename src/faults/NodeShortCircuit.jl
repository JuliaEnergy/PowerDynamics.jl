"""
```Julia
    NodeShortCircuit(;node_number,Y,tspan_fault)
```
# Keyword Arguments
- `node`: number or name of the node
- `Y`: admittance of the short circuit
- `tspan_fault`: short circuit timespan
"""

Base.@kwdef struct NodeShortCircuit <:AbstractPerturbation
    node
    Y
    tspan_fault
    var = :Y_n
end

function (nsc::NodeShortCircuit)(powergrid)
    #mapField(powergrid, nsc, Y_n -> nsc.Y)
    typestable_field_update(powergrid, nsc.node, nsc.var, nsc.Y)
end

#NodeShortCircuit(node,Y,tspan_fault)=NodeShortCircuit(;node,Y,tspan_fault)


#simulate(nsc::NodeShortCircuit, op::State, timespan) = simulate(nsc, op.grid, op.vec, timespan)

export NodeShortCircuit
export NodeShortCircuitError
#export simulate
