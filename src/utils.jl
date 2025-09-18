struct CustomMetadata end

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
    sys = ModelingToolkit.setmetadata(sys, CustomMetadata, mc.metadata)
    sys
end

function get_custom_metadata(sys::System)
    md = ModelingToolkit.getmetadata(sys, CustomMetadata, nothing)
    if isnothing(md)
        error("No custom metadata attached to the system. Use `@attach_metadata!` to attach metadata to the @mtkmodel")
    end
    md
end

get_custom_metadata(sys::System, key) = get_custom_metadata(sys)[key]

"""
    @attach_metadata! Model metadata

Allows you to attach additonal metadata to a `Model` which was previously defined using `@mtkmodel`.
The metadata needs to be in the form of a named tuple `(; name=..., field1=..., field2=...)`.

If `name` is present in the metadata, it will be used as the default name of the system and stripped from the metadata.
The rest of the named tuple will be attachde to the `System`s metadata.
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


function freep(sys::System)
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

"""
    load_pdtesting()

Internal function to dynamicially load the PowerDynamicsTesting code into Main.
If called in an interactive context and Revise.jl is available, it will use Revise.includet to load the code.
"""
function load_pdtesting(force=false)
    if !force && isdefined(Main, :PowerDynamicsTesting)
        println("PowerDynamicsTesting already loaded into Main")
        return
    end
    path = pdtesting_path()
    @eval Main begin
        if isinteractive() && isdefined(Main, :Revise)
            println("Loading PowerDynamicsTesting into Main with Revise")
            Revise.includet($path)
            dir = dirname($path)
            for file in readdir(joinpath(dir))
                if endswith(file, ".jl") && file != "PowerDynamicsTesting.jl"
                    path = PowerDynamics.pdtesting_path()
                    println(" -> track $file")
                    Revise.track(Main.PowerDynamicsTesting, joinpath(dir, file))
                end
            end
        else
            println("Loading PowerDynamicsTesting into Main")
            include($path)
        end
    end
end
pdtesting_path() = joinpath(pkgdir(@__MODULE__), "test", "PowerDynamicsTesting", "PowerDynamicsTesting.jl")
