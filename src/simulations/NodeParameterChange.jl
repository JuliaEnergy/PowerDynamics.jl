using OrdinaryDiffEq: ODEProblem, Rodas4, init, solve!, step!, reinit!, savevalues!, u_modified!
using Setfield

"""
```Julia
NodeParameterChange(;node_number,fraction,tspan_fault)
```
# Keyword Arguments
- `node`: number or name of the node
- `var_new`: value of power during fault event
- `tspan_fault`: PowerPerturbation timespan
- `var`:  symbol of parameter that is perturbed
"""

Base.@kwdef struct NodeParameterChange <:AbstractNodePerturbation
    node
    var_new
    tspan_fault
    var 
end

function (gp::NodeParameterChange)(powergrid)
    mapField(powergrid, gp, p -> p*0.0+convert(Float64, gp.var_new))
end


simulate(gp::NodeParameterChange, op::State, timespan) = simulate(gp, op.grid, op.vec, timespan)


export NodeParameterChange
export simulate
