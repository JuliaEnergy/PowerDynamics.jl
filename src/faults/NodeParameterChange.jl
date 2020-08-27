"""
```Julia
NodeParameterChange(;node_number,fraction,tspan_fault)
```
# Keyword Arguments
- `node`: number or name of the node
- `value`: value of power during fault event
- `tspan_fault`: PowerPerturbation timespan
- `var`:  symbol of parameter that is perturbed
"""

Base.@kwdef struct NodeParameterChange <:AbstractNodePerturbation
    node
    value
    tspan_fault::Tuple
    var::Symbol
end

function (gp::NodeParameterChange)(powergrid)
    #mapField(powergrid, gp, p -> p*0.0+convert(Float64, gp.value))
    typestable_field_update(powergrid, gp.node, gp.var, gp.value)
end


#simulate(gp::NodeParameterChange, op::State, timespan) = simulate(gp, op.grid, op.vec, timespan)

"Error to be thrown if something goes wrong during node parameter perturbation"
struct NodePerturbationError <: PowerDynamicsError
    msg::String
end

export NodeParameterChange
#export simulate
export NodePerturbationError
