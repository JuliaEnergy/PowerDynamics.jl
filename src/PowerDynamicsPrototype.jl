module PowerDynamicsPrototype

using NetworkDynamics
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as Dt
using Graphs: SimpleGraph, add_edge!, nv

include("utils.jl")
include("library.jl")

####
#### Synthetic Power Grids Compat
####
abstract type NodeWrapper end
abstract type EdgeWrapper end

export tocomponent
_tocomplex(x) = x[1] + x[2]*im
_toreal(x) = (real(x), imag(x))

export StaticLine
@kwdef struct StaticLine <: EdgeWrapper
    from::Int
    to::Int
    Y::Complex{Float64}
end
function _staticlinef(e, src, dst, p, t)
    u_src = _tocomplex(src)
    u_dst = _tocomplex(dst)
    Y = _tocomplex(p)
    i_dst = Y * (u_dst - u_src)
    e[1], e[2] = _toreal(i_dst)
    nothing
end
function tocomponent(s::StaticLine)
    StaticEdge(_staticlinef; sym=[:i_r, :i_i], psym=[:Y_r=>real(s.Y), :Y_i=>imag(s.Y)],
        coupling=AntiSymmetric(), name=:StaticLine)
end
#=
s = StaticLine(from=1, to=2, Y=1.0-im)
tocomponent(s)
=#

export PiModelLine
@kwdef struct PiModelLine <: EdgeWrapper
    from::Int
    to::Int
    y::Complex{Float64}
    y_shunt_km::Complex{Float64}
    y_shunt_mk::Complex{Float64}
    t_km::Float64 = 1
    t_mk::Float64 = 1
end
function _pimodellinef(e, src, dst, p, t)
    u_src = _tocomplex(src)
    u_dst = _tocomplex(dst)

    y = _tocomplex(view(p,1:2))
    y_shunt_km = _tocomplex(view(p,3:4))
    y_shunt_mk = _tocomplex(view(p,5:6))
    t_km = p[7]
    t_mk = p[8]

    Π = SMatrix{2,2}(
        - abs2(t_km) * (y + y_shunt_km),
        - conj(t_mk) * t_km * y,
        conj(t_km) * t_mk * y,
        abs2(t_mk) * (y + y_shunt_mk)
    )

    i_src, i_dst = Π * SVector{2}(u_src, u_dst)

    e[1], e[2] = _toreal( i_dst)
    e[3], e[4] = _toreal(-i_src)
    nothing
end
function tocomponent(s::PiModelLine)
    StaticEdge(_pimodellinef; sym=[:i_dst_r, :i_dst_i, :i_src_r, :i_src_i],
        psym=[
            :Y_r=>real(s.y), :Y_i=>imag(s.y),
            :y_shunt_km_r=>real(s.y_shunt_km), :y_shunt_km_i=>imag(s.y_shunt_km),
            :y_shunt_mk_r=>real(s.y_shunt_mk), :y_shunt_mk_i=>imag(s.y_shunt_mk),
            :t_km=>s.t_km,
            :t_mk=>s.t_mk
        ],
        coupling=Fiducial(), name=:PiModelLine)
end
#=
s = PiModelLine(from=1, to=2, y=1.0-im, y_shunt_km=1.0-im, y_shunt_mk=1.0-im, t_km=1.0, t_mk=1.0)
tocomponent(s)
=#

export PVAlgebraic
@kwdef struct PVAlgebraic <: NodeWrapper
    P::Float64
    V::Float64
end
function _pvalgebraicf(dx, x, edges, (P,V), t)
    i = _tocomplex(edges)
    u = _tocomplex(x)
    v = abs(u)
    p = real(u*conj(i))
    du = (v-V) + im*(p-P)
    dx[1], dx[2] = _toreal(du)
    nothing
end
function tocomponent(n::PVAlgebraic)
    ODEVertex(_pvalgebraicf; sym=[:u_r=>1, :u_i=>0], psym=[:P=>n.P, :V=>n.V], mass_matrix=zeros(2),
        name=:PVAlgebraic)
end
#=
n = PVAlgebraic(P=1.0, V=1.0)
tocomponent(n)
=#

export PQAlgebraic
@kwdef struct PQAlgebraic <: NodeWrapper
    P::Float64
    Q::Float64
end
function _pqalgebraicf(dx, x, edges, (P,Q), t)
    i = _tocomplex(edges)
    u = _tocomplex(x)
    s = u*conj(i)
    du = complex(P, Q) - s
    dx[1], dx[2] = _toreal(du)
    nothing
end
function tocomponent(n::PQAlgebraic)
    ODEVertex(_pqalgebraicf; sym=[:u_r=>1, :u_i=>0], psym=[:P=>n.P, :Q=>n.Q], mass_matrix=zeros(2),
        name=:PQAlgebraic)
end
#=
n = PVAlgebraic(P=1.0, V=1.0)
tocomponent(n)
=#

export SlackAlgebraic
@kwdef struct SlackAlgebraic <: NodeWrapper
    U::Complex{Float64}
end
function _slackf(dx, x, edges, p, t)
    u = _tocomplex(x)
    Uref = _tocomplex(p)
    du = u - Uref
    dx[1], dx[2] = _toreal(du)
    nothing
end
function tocomponent(n::SlackAlgebraic)
    ODEVertex(_slackf; sym=[:u_r=>real(n.U), :u_i=>imag(n.U)],
        psym=[:u_ref_r=>real(n.U), :u_ref_i=>imag(n.U)], mass_matrix=zeros(2),
        name=:SlackAlgebraic)
end
#=
n = SlackAlgebraic(U=1.0-im)
tocomponent(n)
=#

include("normalform.jl")


function NetworkDynamics.Network(nodes, lines)
    g = SimpleGraph(length(nodes))
    for l in lines
        add_edge!(g, l.from, l.to)
    end
    Network(g, tocomponent.(nodes), tocomponent.(lines),
        execution=PolyesterExecution{false}(), aggregator=PolyesterAggregator(+))
end

include("power_flow.jl")
include("initialization.jl")


end
