using OrdinaryDiffEq: ODEProblem, Rodas4, init, solve!, step!, reinit!, savevalues!, u_modified!
using Setfield

"""
```Julia
NodePerturbation(;node_number,tspan_fault)
```
# Keyword Arguments
- `node_number`: number  of the node
- `tspan_fault`: PowerPerturbation timespan
"""
Base.@kwdef struct NodePerturbation
    node
    tspan_fault
end

"Error to be thrown if something goes wrong during power perturbation"
struct NodePerturbationError <: PowerDynamicsError
    msg::String
end
