export NormalFormSystem
using LinearAlgebra: dot, ⋅
using StaticArrays
using NetworkDynamics: subscript

export NormalForm
struct NormalForm{x_dims} <: NodeWrapper
    P::Float64
    Q::Float64
    V::Float64
    Bᵤ::SMatrix{1,x_dims,ComplexF64}
    Cᵤ::ComplexF64
    Gᵤ::ComplexF64
    Hᵤ::ComplexF64
    Bₓ::SMatrix{x_dims,x_dims,Float64}
    Cₓ::SVector{x_dims,Float64}
    Gₓ::SVector{x_dims,Float64}
    Hₓ::SVector{x_dims,Float64}
    Y_n::Complex{Float64}
end

NormalForm(; P, Q, V, Bᵤ=[], Cᵤ, Gᵤ, Hᵤ, Bₓ=[], Cₓ=[], Gₓ=[], Hₓ=[], Y_n = 0.0im) = NormalForm(P, Q, V, Bᵤ, Cᵤ, Gᵤ, Hᵤ, Bₓ, Cₓ, Gₓ, Hₓ, Y_n)

function NormalForm(P, Q, V, Bᵤ, Cᵤ, Gᵤ, Hᵤ, Bₓ, Cₓ, Gₓ, Hₓ, Y_n)
    x_dims = length(Cₓ)
    Bᵤ = SMatrix{1, x_dims, ComplexF64}(Bᵤ)
    Bₓ = SMatrix{x_dims, x_dims, Float64}(Bₓ)
    Cₓ = SVector{x_dims, Float64}(Cₓ)
    Gₓ = SVector{x_dims, Float64}(Gₓ)
    Hₓ = SVector{x_dims, Float64}(Hₓ)
    NormalForm(P, Q, V, Bᵤ, Cᵤ, Gᵤ, Hᵤ, Bₓ, Cₓ, Gₓ, Hₓ, Y_n)
end

struct NFRhs{N} end

idx_P(N) = 1
idx_Q(N) = 2
idx_V(N) = 3
idx_Bᵤ_re(N) = 4 : 4 + N - 1
idx_Bᵤ_im(N) = 4 + N : 4 + 2*N - 1
idx_Cᵤ(N) = 4 + 2*N : 4 + 2*N + 1
idx_Gᵤ(N) = 4 + 2*N + 2 : 4 + 2*N + 3
idx_Hᵤ(N) = 4 + 2*N + 4 : 4 + 2*N + 5
idx_Bₓ(N) = 4 + 2*N + 6 : 4 + 2*N + 6 + N^2 - 1
idx_Cₓ(N) = 4 + 2*N + 6 + N^2 : 4 + 2*N + 6 + N^2 + N - 1
idx_Gₓ(N) = 4 + 2*N + 6 + N^2 + N : 4 + 2*N + 6 + N^2 + N + N - 1
idx_Hₓ(N) = 4 + 2*N + 6 + N^2 + N + N : 4 + 2*N + 6 + N^2 + N + N + N - 1
idx_Yn(N) = 4 + 2*N + 6 + N^2 + N + N + N : 4 + 2*N + 6 + N^2 + N + N + N + 1
_getpdim(N) = 3 + 4*2 + 5*N + N^2

#=
using OpPoDyn: idx_P , idx_Q , idx_V , idx_Bᵤ_re, idx_Bᵤ_im, idx_Cᵤ, idx_Gᵤ, idx_Hᵤ, idx_Bₓ, idx_Cₓ, idx_Gₓ, idx_Hₓ, idx_Yn

N = 3
idxs = Int[]
append!(idxs, idx_P(N))
append!(idxs, idx_Q(N))
append!(idxs, idx_V(N))
append!(idxs, idx_Bᵤ_re(N))
append!(idxs, idx_Bᵤ_im(N))
append!(idxs, idx_Cᵤ(N))
append!(idxs, idx_Gᵤ(N))
append!(idxs, idx_Hᵤ(N))
append!(idxs, idx_Bₓ(N))
append!(idxs, idx_Cₓ(N))
append!(idxs, idx_Gₓ(N))
append!(idxs, idx_Hₓ(N))
append!(idxs, idx_Yn(N))
@assert idxs == 1:(3 + 4*2 + 5*N + N^2)
=#

function (::NFRhs{N})(dz, z, esum, p, t) where {N}
    P   = p[1]
    Q   = p[2]
    V   = p[3]
    Bᵤ  = SMatrix{1,N}(view(p, idx_Bᵤ_re(N))) + im * SMatrix{1,N}(view(p, idx_Bᵤ_re(N)))
    Cᵤ  = _tocomplex(view(p, idx_Cᵤ(N)))
    Gᵤ  = _tocomplex(view(p, idx_Gᵤ(N)))
    Hᵤ  = _tocomplex(view(p, idx_Hᵤ(N)))
    Bₓ  = SMatrix{N, N}(view(p, idx_Bₓ(N)))
    Cₓ  = SVector{N}(view(p, idx_Cₓ(N)))
    Gₓ  = SVector{N}(view(p, idx_Gₓ(N)))
    Hₓ  = SVector{N}(view(p, idx_Hₓ(N)))
    Y_n = _tocomplex(view(p, idx_Yn(N)))

    i = _tocomplex(esum) + Y_n * (z[1] + z[2] * 1im)
    u = z[1] + z[2] * im
    x = view(z, 3:3+N-1)
    s = u * conj(i)

    δp = real(s) - P
    δq = imag(s) - Q
    δv2 = abs2(u) - V^2

    du = (conj(Bᵤ) ⋅ x + Cᵤ * δv2 + Gᵤ * δp + Hᵤ * δq) * u  # conj(Bᵤ) because Julia's dot-product conjugates the first vector/matrix
    dx = Bₓ * x + Cₓ * δv2 + Gₓ * δp + Hₓ * δq

    dz[1] = real(du)
    dz[2] = imag(du)
    dz[3:end] .= real(dx)

    return nothing
end

function tocomponent(nf::NormalForm{N}) where {N}
    sym = [:u_r, :u_i]
    append!(sym, [Symbol("x"*subscript(i)) for i in 1:N])

    psym = [:P => nf.P, :Q => nf.Q, :V => nf.V]
    append!(psym, [Symbol("Bᵤ_r"*subscript(i)) => real(nf.Bᵤ[i]) for i in 1:N])
    append!(psym, [Symbol("Bᵤ_i"*subscript(i)) => real(nf.Bᵤ[i]) for i in 1:N])
    append!(psym, [:Cᵤ_r => real(nf.Cᵤ), :Cᵤ_i => imag(nf.Cᵤ)])
    append!(psym, [:Gᵤ_r => real(nf.Gᵤ), :Gᵤ_i => imag(nf.Gᵤ)])
    append!(psym, [:Hᵤ_r => real(nf.Hᵤ), :Hᵤ_i => imag(nf.Hᵤ)])
    append!(psym, [Symbol("Bₓ"*subscript(i)*subscript(j)) => nf.Bₓ[i,j] for i in 1:N, j in 1:N])
    append!(psym, [Symbol("Cₓ"*subscript(i)) => nf.Cₓ[i] for i in 1:N])
    append!(psym, [Symbol("Gₓ"*subscript(i)) => nf.Gₓ[i] for i in 1:N])
    append!(psym, [Symbol("Hₓ"*subscript(i)) => nf.Hₓ[i] for i in 1:N])
    append!(psym, [:Y_N_r => real(nf.Y_n), :Y_N_i => imag(nf.Y_n)])
    @assert length(psym) == _getpdim(N)

    ODEVertex(NFRhs{N}(); sym, psym, name=Symbol("NormalForm{$N}"))
end

# @mtkmodel NormalFormSystem begin
#     @extend DQBus()
#     @variables begin
#         ΔP(t), [description = "active power deviation"]
#         ΔQ(t), [description = "active power deviation"]
#         ΔV2(t), [description = "squared voltage deviation"]
#         i_int_r(t), [description = "internal d-current"]
#         i_int_i(t), [description = "internal q-current"]
#     end
#     @structural_parameters begin
#         N
#     end
#     @parameters begin
#         P_ref, [description = "Active power reference"]
#         Q_ref, [description = "Reactive power reference"]
#         V_ref = 1.0, [description = "Voltage reference"]
#         Y_r, [description = "real part of self inductance?"]
#         Y_i, [description = "imag part of self inductance?"]
#         Bᵤ[1:1,1:N]::Complex
#         # Cᵤ::Complex
#         # Gᵤ::Complex
#         # Hᵤ::Complex
#         # Bₓ[1:N,1:N]
#         # Cₓ[1:N]
#         # Gₓ[1:N]
#         # Hₓ[1:N]
#     end
#     @equations begin
#         i_int_r ~ i_r - Y_i*u_i + Y_r*u_r
#         i_int_i ~ i_r + Y_i*u_r + Y_r*u_i
#         ΔP ~ ( i_int_i*u_i + i_int_r*u_r) - P_ref
#         ΔQ ~ (-i_int_i*u_r + i_int_r*u_i) - Q_ref
#         ΔV2 ~ (u_r^2 + u_i^2) - V_ref^2
#         # Dt(u_r) ~ real( dot(conj(Bᵤ), x) + Cᵤ * ΔV2 + Gᵤ * ΔP + Hᵤ * ΔQ) * (u_r + im*u_i) )  # conj(Bᵤ) because Julia's dot-product conjugates the first vector/matrix
#         # Dt(u_i) ~ imag( dot(conj(Bᵤ), x) + Cᵤ * ΔV2 + Gᵤ * ΔP + Hᵤ * ΔQ) * (u_r + im*u_i) )  # conj(Bᵤ) because Julia's dot-product conjugates the first vector/matrix
#         # Dt(x) ~ Bₓ * x + Cₓ * ΔV2 + Gₓ * ΔP + Hₓ * ΔQ
#     end
# end
#= generate the equations
@variables t u_r(t) u_i(t) i_r(t) i_i(t) i_int_r(t) i_int_i(t) P(t)  Q(t)
@parameters Y_r Y_i

i_int_r ~ simplify(i_r + real((Y_r+im*Y_i)*(u_r+im*u_i)))
i_int_i ~ simplify(i_r + imag((Y_r+im*Y_i)*(u_r+im*u_i)))
P ~ simplify( real( conj(i_int_r+im*i_int_i)*(u_r+im*u_i) ) )
Q ~ simplify( imag( conj(i_int_r+im*i_int_i)*(u_r+im*u_i) ) )
=#
