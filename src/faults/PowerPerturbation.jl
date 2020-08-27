"""
```Julia
PowerPerturbation(;node_number,power_new,tspan_fault,var)
```
# Keyword Arguments
- `node_number`: number or name of the node
- `power_new`: value of power during fault event
- `tspan_fault`: PowerPerturbation timespan
- `var`: parameter symbol on the node that represents Power, default is :P
"""

Base.@kwdef struct PowerPerturbation <:AbstractNodePerturbation
    node
    power_new
    tspan_fault
    var = :P
end

function (pd::PowerPerturbation)(powergrid)
    #mapField(powergrid, pd, p -> p*0.0+convert(Float64, pd.power_new))
    typestable_field_update(powergrid, pd.node, pd.var, pd.power_new)
end

#PowerPerturbation(node,power_new,tspan_fault)=PowerPerturbation(;node,power_new,tspan_fault)

#simulate(pd::PowerPerturbation, op::State, timespan) = simulate(pd, op.grid, op.vec, timespan)

"Error to be thrown if something goes wrong during power perturbation"
struct PowerPerturbationError <: PowerDynamicsError
    msg::String
end

export PowerPerturbation
#export simulate
export PowerPerturbationError