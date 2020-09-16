@doc """
```Julia
NodeParameterChange(;node, value, tspan_fault, var)
```
# Keyword Arguments
- `node`: number or name of the node
- `value`: value of power during fault event
- `tspan_fault`: PowerPerturbation timespan
- `var`:  symbol of parameter that is perturbed
"""
Base.@kwdef struct NodeParameterChange <:AbstractPerturbation
    node
    value
    tspan_fault::Tuple
    var::Symbol
end

function (gp::NodeParameterChange)(powergrid)
    typestable_node_field_update(powergrid, gp.node, gp.var, gp.value)
end

"""
Error to be thrown if something goes wrong during node parameter perturbation.
"""
struct NodePerturbationError <: PowerDynamicsError
    msg::String
end

export NodeParameterChange
export NodePerturbationError
