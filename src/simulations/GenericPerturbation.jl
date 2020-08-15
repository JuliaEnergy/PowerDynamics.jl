using OrdinaryDiffEq: ODEProblem, Rodas4, init, solve!, step!, reinit!, savevalues!, u_modified!
using Setfield

"""
```Julia
GemericPerturbation(;node_number,fraction,tspan_fault)
```
# Keyword Arguments
- `node`: number or name of the node
- `var_new`: value of power during fault event
- `tspan_fault`: PowerPerturbation timespan
- `var`:  symbol of parameter that is perturbed
"""

Base.@kwdef struct GenericPerturbation <:NodePerturbation
    node
    var_new
    tspan_fault
    var 
end

function (gp::GenericPerturbation)(powergrid)
    mapField(powergrid, gp, p -> p*0.0+convert(Float64, gp.var_new))
end

#GenericPerturbation(node,var_new,tspan_fault,var)=GenericPerturbation(;node,var_new,tspan_fault,var)

simulate(gp::GenericPerturbation, op::State, timespan) = simulate(gp, op.grid, op.vec, timespan)


export GenericPerturbation
export simulate
