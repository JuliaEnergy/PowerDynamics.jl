"""
```Julia
PowerPerturbation(;node_number,fault_power,tspan_fault,var)
```
# Keyword Arguments
- `node_number`: number or name of the node
- `fault_power`: value of power during fault event
- `tspan_fault`: PowerPerturbation timespan
- `var`: parameter symbol on the node that represents Power, default is :P
"""
Base.@kwdef struct PowerPerturbation <:AbstractPerturbation
    node
    fault_power
    tspan_fault
    var = :P
end

function (pd::PowerPerturbation)(powergrid)
    typestable_node_field_update(powergrid, pd.node, pd.var, pd.fault_power)
end

export PowerPerturbation
