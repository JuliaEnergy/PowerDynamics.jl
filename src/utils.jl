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
