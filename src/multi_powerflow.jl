abstract type PowerFlowType end

@kwdef struct Junction <: PowerFlowType
    subs::Vector{PowerFlowType} = Vector{PowerFlowType}()
    metadata::Dict{Symbol,Any} = Dict{Symbol,Any}()
end
@kwdef struct SlackType <: PowerFlowType
    V::Float64 = 1.0
    δ::Float64 = 0.0
    subs::Vector{PowerFlowType} = Vector{PowerFlowType}()
    metadata::Dict{Symbol,Any} = Dict{Symbol,Any}()
end
@kwdef struct PVType <: PowerFlowType
    P::Float64 = 0.0
    V::Float64 = 1.0
    subs::Vector{PowerFlowType} = Vector{PowerFlowType}()
    metadata::Dict{Symbol,Any} = Dict{Symbol,Any}()
end
@kwdef struct PQYType <: PowerFlowType
    P::Float64 = 0.0
    Q::Float64 = 0.0
    G::Float64 = 0.0
    B::Float64 = 0.0
    subs::Vector{PowerFlowType} = Vector{PowerFlowType}()
    metadata::Dict{Symbol,Any} = Dict{Symbol,Any}()
end

PQType(; kwargs...) = PQYType(; kwargs...)
YType(; kwargs...) = PQYType(; kwargs...)


# S + S
function combine(sA::SlackType, sB::SlackType)
    angle = sA.δ == sB.δ ? sA.δ : error("Incompatible voltage angles for slack merge: $(sA.δ) vs $(sB.δ)!")
    subs = combined_subs(sA, sB)
    metadata = combined_metadata(sA, sB)
    SlackType(combined_voltage(sA, sB), angle, subs, metadata)
end
# S + PV
function combine(s::SlackType, pv::PVType)
    subs = combined_subs(s, pv)
    metadata = combined_metadata(s, pv)
    SlackType(combined_voltage(s, pv), s.δ, subs, metadata)
end
combine(pv::PVType, s::SlackType) = combine(s, pv)
# S + PQY (Y component absorbed by fixed voltage)
function combine(s::SlackType, pqy::PQYType)
    subs = combined_subs(s, pqy)
    metadata = combined_metadata(s, pqy)
    SlackType(s.V, s.δ, subs, metadata)
end
combine(pqy::PQYType, s::SlackType) = combine(s, pqy)

# PV + PV
function combine(pvA::PVType, pvB::PVType)
    subs = combined_subs(pvA, pvB)
    metadata = combined_metadata(pvA, pvB)
    PVType(pvA.P + pvB.P, combined_voltage(pvA, pvB), subs, metadata)
end
# PV + PQY (Y component absorbed by fixed voltage)
function combine(pv::PVType, pqy::PQYType)
    subs = combined_subs(pv, pqy)
    metadata = combined_metadata(pv, pqy)
    # Y admittance contributes constant power at fixed voltage V
    P_total = pv.P + pqy.P + pqy.G * pv.V^2
    PVType(P_total, pv.V, subs, metadata)
end
combine(pqy::PQYType, pv::PVType) = combine(pv, pqy)

combined_voltage(v1::PowerFlowType, v2::PowerFlowType) = combined_voltage(v1.V, v2.V)
function combined_voltage(v1, v2)
    isnan(v1) && !isnan(v2) && return v2
    !isnan(v1) && isnan(v2) && return v1
    isapprox(v1, v2; rtol=1e-5, atol=1e-8) && return v1
    error("Incompatible voltage setpoints: $(str_significant(v1)) vs $(str_significant(v2))!")
end

# PQY + PQY
function combine(pqA::PQYType, pqB::PQYType)
    subs = combined_subs(pqA, pqB)
    metadata = combined_metadata(pqA, pqB)
    PQYType(pqA.P + pqB.P, pqA.Q + pqB.Q, pqA.G + pqB.G, pqA.B + pqB.B, subs, metadata)
end

combined_subs(A::PowerFlowType, B::PowerFlowType) = vcat(A.subs, B.subs)
combined_metadata(A::PowerFlowType, B::PowerFlowType) = merge(A.metadata, B.metadata)

const _PFTYPE_MODEL_CACHE = Dict{Symbol,VertexModel}()
function get_cached_vertex_model(type::Symbol)
    get!(_PFTYPE_MODEL_CACHE, type) do
        if type == :slack
            @named slack = Library.VδConstraint()
            compile_bus(MTKBus(slack))
        elseif type == :pv
            @named pv = Library.PVConstraint()
            compile_bus(MTKBus(pv))
        elseif type == :pq
            @named pq = Library.PQConstraint()
            compile_bus(MTKBus(pq))
        elseif type == :y
            @named shunt = Library.ConstantYLoad()
            compile_bus(MTKBus(shunt))
        elseif type == :pqy
            @named pq = Library.PQConstraint()
            @named shunt = Library.ConstantYLoad(allow_zero_conductance=true)
            compile_bus(MTKBus([pq, shunt]), assume_io_coupling=true)
        elseif type == :junction
            compile_bus(MTKBus())
        else
            error("Unknown vertex blueprint: $type")
        end
    end
end
function to_vertexmodel(pqy::PQYType, name)
    hasY = !(iszero(pqy.G) && iszero(pqy.B))
    hasPQ = !(iszero(pqy.P) && iszero(pqy.Q))
    !hasY && !hasPQ && error("PQYType must have at least one non-zero component!")
    # If no shunt admittance, use simple PQ model
    vm = if hasPQ && !hasY
        get_cached_vertex_model(:pq)
    elseif !hasPQ && hasY
        get_cached_vertex_model(:y)
    else
        get_cached_vertex_model(:pqy)
    end
    # set name
    vm = VertexModel(vm; name)

    hasPQ && set_default!(vm, :pq₊P, pqy.P)
    hasPQ && set_default!(vm, :pq₊Q, pqy.Q)
    hasY && set_default!(vm, :shunt₊G, pqy.G)
    hasY && set_default!(vm, :shunt₊B, pqy.B)

    return vm
end
function to_vertexmodel(pv::PVType, name)
    vm = get_cached_vertex_model(:pv)
    vm = VertexModel(vm; name)
    set_default!(vm, :pv₊P, pv.P)
    set_default!(vm, :pv₊V, pv.V)
    return vm
end
function to_vertexmodel(s::SlackType, name)
    vm = get_cached_vertex_model(:slack)
    vm = VertexModel(vm; name)
    set_default!(vm, :slack₊V, s.V)
    set_default!(vm, :slack₊δ, s.δ)
    return vm
end
function to_vertexmodel(junc::Junction, name)
    vm = get_cached_vertex_model(:junction)
    vm = VertexModel(vm; name)
    return vm
end

function to_vertexmodel(dict::AbstractDict, vm)
    overall = reduce(combine, values(dict))
    basemod = to_vertexmodel(overall, vm.name)
    VertexModel
    basemod.g
    basemod
end

# (out, u, ivec, p, t) -> function
#     _uvec = view(out, 1:2) # use obs output buffer to get voltage
#     NetworkDynamics.apply_compg!(fftype(basemod), basemod.g, (_uvec, ), u, (ivec, ), p, t)
#     i_r, i_i = ivec
#     u_r, u_i = _uvec

# end
