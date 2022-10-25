using StaticArrays

@doc doc"""
```Julia
NormalForm(P, Q, V, Bᵤ, Cᵤ, Gᵤ, Hᵤ, Bₓ, Cₓ, Gₓ, Hₓ, Y_n, x_dims)
```

An implementation of the normal form for grid-forming actors following the publication:

    Kogler, R., Plietzsch, A., Schultz, P., & Hellmann, F. (2022).
    Normal Form for Grid-Forming Power Grid Actors. PRX Energy, 1(1), 013008.
    https://doi.org/10.1103/PRXEnergy.1.013008

# Mathematical Representation
The dynamics of the normal form is given by
```math
    \frac{\dot{u}}{u} = B_u x + C_u \delta v^2 + G_u \delta p + H_u \delta q\\
    \dot{x} = B_x x + + C_x \delta v^2 + G_x \delta p + H_x \delta q
```
where ``u`` is the complex voltage and ``x`` is a vector of internal variables
and ``\delta v^2``, ``\delta p`` and ``\delta q`` denote the deviations of 
squared voltage amplitude, active and reactive power from their setpoint values.
The internal variables ``x`` are defined such that their setpoint values are zero.

# Keyword Arguments
- `B_u`: complex parameter (matrix) 
- `C_u`: complex parameter
- `G_u`: complex parameter
- `H_u`: complex parameter
- `B_x`: real parameter (matrix) 
- `C_x`: real parameter (vector)
- `G_x`: real parameter (vector)
- `H_x`: real parameter (vector)

Note: The parameters `A_x` and `A_u` in the paper are always zero in the implementation.
(``A_x=0`` as ``x=0`` in the setpoint and ``A_u=j\omega_s`` with ``\omega_s=0`` in the co-rotating system.)

"""
struct NormalForm{x_dims} <: AbstractNode
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

function construct_vertex(nf::NormalForm)

    sym = symbolsof(nf)
    dim = dimension(nf)
    mass_matrix = ones(Int64, dim, dim) |> Diagonal

    P = nf.P
    Q = nf.Q
    V = nf.V
    
    Bᵤ = nf.Bᵤ
    Cᵤ = nf.Cᵤ
    Gᵤ = nf.Gᵤ
    Hᵤ = nf.Hᵤ
    Bₓ = nf.Bₓ
    Cₓ = nf.Cₓ
    Gₓ = nf.Gₓ
    Hₓ = nf.Hₓ

    Y_n = nf.Y_n

    rhs! = function (dz, z, edges, p, t)
        i = total_current(edges) + Y_n * (z[1] + z[2] * 1im)
        u = z[1] + z[2] * im
        @views x = z[3:end]  # @views is needed to avoid allocations
        s = u * conj(i)

        δp = real(s) - P
        δq = imag(s) - Q
        δv2 = abs2(u) - V^2

        du = (conj(Bᵤ) ⋅ x + Cᵤ * δv2 + Gᵤ * δp + Hᵤ * δq) * u  # conj(Bᵤ) because Julia's dot-product conjugates the first vector/matrix
        dx = (Bₓ * x + Cₓ * δv2 + Gₓ * δp + Hₓ * δq)
        
        dz[1] = real(du)  
        dz[2] = imag(du)
        dz[3:end] = real(dx)

        return nothing
    end

    ODEVertex(rhs!, dim, mass_matrix, sym)
end

function symbolsof(nf::NormalForm)
    x_dims = length(nf.Cₓ)
    symbols = [:u_r, :u_i]
    append!(symbols, [Symbol("x_$i") for i in 1:x_dims])
    return symbols
end

function dimension(nf::NormalForm)
    x_dims = length(nf.Cₓ)
    return x_dims + 2
end

export NormalForm