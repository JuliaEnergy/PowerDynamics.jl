using BlockSystems

export IONode

"""
    IONode{T,M} <: AbstractNode
    IONode(bp::BlockPara)
    IONode(blk::IOBlock, parameters::Dict)

Create an `IONode` based on an `IOBlock` and a parameter dict.

The block needs to fulfil the folowing interface:
- inputs:  `i_r`, `i_i`
- outputs: `u_r`, `u_i`

The block gets the flow sum of the connected lines as an input. It should calculate
the resulting node voltage as an output (either per DGL or constraint).

The parameters should be provided as a `Dict{Symbol,Float}` such as

    p = Dict(:a => 1.0,
             :b => -π)
"""
struct IONode{T,M} <: AbstractNode
    block::IOBlock
    generated::T
    parameter_names::Vector{Symbol}
    parameters::Vector{Float64}
    mass_matrix::M
end

function IONode(bp::BlockPara)
    blk, para = bp.block, bp.para

    # BlockSpec: blk must by of type (i_r, i_i) ↦ (u_r, u_i)
    if !fulfills(blk, BlockSpec([:i_r, :i_i], [:u_r, :u_i]))
        throw(ArgumentError("Block has to follow PowerDynamics i/o conventions!"))
    end

    # create vectors from key, value pairs
    p_keys = collect(keys(para))
    p_vals = collect(Float64, values(para))

    gen = generate_io_function(blk,
                               f_states=[blk.u_r, blk.u_i],
                               f_inputs=[blk.i_r, blk.i_i],
                               f_params=p_keys, warn=false,
                               type=:ode);

    IONode(blk, gen, p_keys, p_vals, gen.massm)
end

IONode(blk::IOBlock, para::Dict) = IONode(BlockPara(blk, para))

# extend the necessary functions for the `AbstractNode` interface
function construct_vertex(ion::IONode)
    gen = ion.generated
    function rhs!(dx, x, edges, _p, t)
        i = total_current(edges)
        gen.f_ip(dx, x, (real(i), imag(i)), ion.parameters, t)
    end
    ODEVertex(f = rhs!, dim = length(gen.states), mass_matrix = gen.massm, sym = Symbol.(gen.states))
end

symbolsof(ionode::IONode) = Symbol.(ionode.generated.states)
dimension(ionode::IONode) = length(ionode.generated.states)
