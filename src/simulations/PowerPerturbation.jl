using OrdinaryDiffEq: ODEProblem, Rodas4, init, solve!, step!, reinit!, savevalues!, u_modified!
using Setfield

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

Base.@kwdef struct PowerPerturbation <:NodePerturbation
    node
    power_new
    tspan_fault
    var = :P
end

function (pd::PowerPerturbation)(powergrid)
    mapField(powergrid, pd, p -> p*0.0+convert(Float64, pd.power_new))
end

#PowerPerturbation(node,power_new,tspan_fault)=PowerPerturbation(;node,power_new,tspan_fault)

simulate(pd::PowerPerturbation, op::State, timespan) = simulate(pd, op.grid, op.vec, timespan)


export PowerPerturbation
export simulate