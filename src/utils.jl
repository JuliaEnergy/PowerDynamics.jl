"""
This type wraps the default constructor of a ModelingToolkit.Model and allows to attach metadata to it.
Two types of metadata is possible:
  - a "name" field which will become the "name" field of the system
  - a named tuple with aribrary fields which will become the "metadata" field of the system
"""
struct ModelMetadataConstructor{F, NT}
    original_constructor::F
    name::Union{Nothing, Symbol}
    metadata::NT
end
function ModelMetadataConstructor(model::ModelingToolkit.Model, metadata::NamedTuple)
    name, finalmeta = _split_name(metadata)
    ModelMetadataConstructor(model.f, name, finalmeta)
end
function ModelMetadataConstructor(model::ModelingToolkit.Model{<:ModelMetadataConstructor}, metadata::NamedTuple)
    mc = model.f
    oldmeta = _full_metadata(mc)
    newmeta = merge(oldmeta, metadata)
    name, finalmeta = _split_name(newmeta)
    ModelMetadataConstructor(mc.original_constructor, name, finalmeta)
end
function _split_name(metadata)
    if haskey(metadata, :name)
        name = metadata.name
        final_metadata = NamedTuple{filter(!isequal(:name), keys(metadata))}(metadata)
    else
        name = nothing
        final_metadata = metadata
    end
    return name, final_metadata
end
_full_metadata(mc::ModelMetadataConstructor) = (; name=mc.name, mc.metadata...)

function (mc::ModelMetadataConstructor)(args...; kwargs...)
    sys = if isnothing(mc.name)
        mc.original_constructor(args...; kwargs...)
    else
        mc.original_constructor(args...; name=mc.name, kwargs...)
    end
    if !isnothing(sys.metadata)
        @warn "Overwriting metadata for $(sys.name): $(sys.metadata). This should not happen, please report issue!"
    end
    @set sys.metadata = mc.metadata
end

"""
    @attach_metadata! Model metadata

Allows you to attach additonal metadata to a `Model` which was previously defined using `@mtkmodel`.
The metadata needs to be in the form of a named tuple `(; name=..., field1=..., field2=...)`.

If `name` is present in the metadata, it will be used as the default name of the system and stripped from the metadata.
The rest of the named tuple will be attachde to the `ODESystem`s metadata.
"""
macro attach_metadata!(model, metadata)
     quote
         $(esc(model)) = _attach_metadata($(esc(model)), $(esc(metadata)))
     end
end
function _attach_metadata(model::ModelingToolkit.Model, metadata::NamedTuple)
    mc = ModelMetadataConstructor(model, metadata)
    @set model.f = mc
end


function freep(sys::ODESystem)
    return filter(p -> !haskey(ModelingToolkit.defaults(sys), p), ModelingToolkit.parameters(sys))
end
function freep(cf::NetworkDynamics.ComponentModel)
    return filter(p -> !NetworkDynamics.has_default(cf, p), NetworkDynamics.psym(cf))
end
function freep(nw::Network)
    setdiff(NetworkDynamics.SII.parameter_symbols(nw), keys(NetworkDynamics.SII.default_values(nw)))
end
function freeu(nw::Network)
    setdiff(NetworkDynamics.SII.variable_symbols(nw), keys(NetworkDynamics.SII.default_values(nw)))
end
function freeu(cf::NetworkDynamics.ComponentModel)
    return filter(u -> !NetworkDynamics.has_default(cf, u), NetworkDynamics.sym(cf))
end

set_voltage!(cf::VertexModel; mag, arg) = set_voltage!(cf, mag * exp(im * arg))
function set_voltage!(cf::VertexModel, c::Complex)
    set_default!(cf, :busbar₊u_r, real(c))
    set_default!(cf, :busbar₊u_i, imag(c))
    c
end

get_voltage(cf::VertexModel) = get_default(cf, :busbar₊u_r) + im * get_default(cf, :busbar₊u_i)

function set_current!(cf::VertexModel; P, Q)
    @assert has_default(cf, :busbar₊u_r) && has_default(cf, :busbar₊u_i)
    u = get_default(cf, :busbar₊u_r) + im * get_default(cf, :busbar₊u_i)
    i = conj((P + im * Q) / u)
    set_current!(cf, -i)
end
# set_current!(cf::VertexModel; mag, arg) = set_current!(cf, mag * exp(im * arg))
function set_current!(cf::VertexModel, c::Complex)
    set_default!(cf, :busbar₊i_r, real(c))
    set_default!(cf, :busbar₊i_i, imag(c))
    c
end

get_current(cf::VertexModel) = -1 * (get_default(cf, :busbar₊i_r) + im * get_default(cf, :busbar₊i_i))
get_power(cf::VertexModel) = get_voltage(cf) * conj(get_current(cf))

function normalize_angle(θ)
    θ = mod2pi(θ)
    θ >= π ? θ - 2π : θ
end
