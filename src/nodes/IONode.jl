using BlockSystems
using ModelingToolkit: getname, value

export IONode

struct IONode{T,M} <: AbstractNode
    block::IOBlock
    generated::T
    parameter_names::Vector{Symbol}
    parameters::Vector{Float64}
    mass_matrix::M
end

function IONode(blk::IOBlock, parameters::Dict)
    # BlockSpec: blk must by of type (i_r, i_i) ↦ (u_r, u_i)
    spec = BlockSpec([:i_r, :i_i], [:u_r, :u_i])
    @assert spec(blk) "Block has to follow PowerDynamics i/o conventions!"

    # parameters may be given as Num oder Symbol types
    p_keys = [k isa Symbol ? k : getname(value(k)) for k in keys(parameters)]
    p_vals = collect(Float64, values(parameters))

    gen = generate_io_function(blk,
                               f_states=[blk.u_r, blk.u_i],
                               f_inputs=[blk.i_r, blk.i_i],
                               f_params=p_keys, warn=false,
                               type=:ode);

    IONode(blk, gen, p_keys, p_vals, gen.massm)
end

function construct_vertex(ion::IONode)
    gen = ion.generated
    function rhs!(dx, x, edges, _p, t)
        i = total_current(edges)
        gen.f_ip(dx, x, (real(i), imag(i)), ion.parameters, t)
    end
    ODEVertex(f! = rhs!, dim = length(gen.states), mass_matrix = gen.massm, sym = Symbol.(gen.states))
end

symbolsof(ionode::IONode) = Symbol.(ionode.generated.states)
dimension(ionode::IONode) = length(ionode.generated.states)
