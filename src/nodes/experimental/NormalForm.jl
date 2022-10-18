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
    \frac{\dot{u}}{u} = B_u x + + C_u \delta v^2 + G_u \delta p + H_u \delta q\\
    \dot{x} = B_x x + + C_x \delta v^2 + G_x \delta p + H_x \delta q
```
where ``u`` is the complex voltage and ``x`` is a vector of internal variables
and ``\delta v^2``, ``\delta p`` and ``\delta q`` denote the deviations of 
squared voltage amplitude, active and reactive power from their setpoint values.

# Keyword Arguments
- `B_u`: complex parameter (matrix) 
- `C_u`: complex parameter
- `G_u`: complex parameter
- `H_u`: complex parameter
- `B_x`: real parameter (matrix) 
- `C_u`: real parameter (vector)
- `G_u`: real parameter (vector)
- `H_u`: real parameter (vector)
- `x_dims`: number of internal variables

"""
struct NormalForm <: AbstractNode
    P
    Q
    V
    Bᵤ
    Cᵤ
    Gᵤ
    Hᵤ
    Bₓ
    Cₓ
    Gₓ
    Hₓ
    Y_n
    x_dims
end

NormalForm(; P, Q, V, Bᵤ, Cᵤ, Gᵤ, Hᵤ, Bₓ, Cₓ, Gₓ, Hₓ, Y_n = 0, x_dims) = NormalForm(P, Q, V, Bᵤ, Cᵤ, Gᵤ, Hᵤ, Bₓ, Cₓ, Gₓ, Hₓ, Y_n, x_dims)

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

    if x_dims == 0  # Case 1: no internal variable

        rhs! = function (dz, z, edges, p, t)
            i = total_current(edges) + Y_n * (z[1] + z[2] * 1im)
            u = z[1] + z[2] * im 
            s = u * conj(i)

            δp = real(s) - P
            δq = imag(s) - Q
            δv2 = abs2(u) - V^2

            du = (Cᵤ * δv2 + Gᵤ * δp + Hᵤ * δq) * u
            
            dz[1] = real(du)  
            dz[2] = imag(du)

            return nothing
        end

    elseif x_dims == 1  # Case 2: one internal variable

        @assert typeof(Bₓ) <: Real "Bₓ must be a real number."
        @assert typeof(Cₓ) <: Real "Cₓ must be a real number."
        @assert typeof(Gₓ) <: Real "Gₓ must be a real number."
        @assert typeof(Hₓ) <: Real "Hₓ must be a real number."

        rhs! = function (dz, z, edges, p, t)
            i = total_current(edges) + Y_n * (z[1] + z[2] * 1im)
            u = z[1] + z[2] * im
            x = z[3]
            s = u * conj(i)

            δp = real(s) - P
            δq = imag(s) - Q
            δv2 = abs2(u) - V^2

            du = (Bᵤ * x + Cᵤ * δv2 + Gᵤ * δp + Hᵤ * δq) * u
            dx = (Bₓ * x + Cₓ * δv2 + Gₓ * δp + Hₓ * δq)
            
            dz[1] = real(du)  
            dz[2] = imag(du)
            dz[3] = real(dx)          

            return nothing
        end

    elseif x_dims > 1 # Case of multiple internal variables

        @assert typeof(Bₓ) <: Matrix "Bₓ must be a matrix."
        @assert all(imag(Bₓ) .== 0) "Bₓ must be real."
        @assert size(Bₓ) == (x_dims,x_dims) "Bₓ parameters have the wrong dimension."
        @assert typeof(Cₓ) <: Array "Cₓ must be an array."
        @assert all(imag(Cₓ) .== 0) "Cₓ must be real."
        @assert length(Cₓ) == x_dims "Cₓ parameters have the wrong dimension."
        @assert typeof(Gₓ) <: Array "Gₓ must be an array."
        @assert all(imag(Gₓ) .== 0) "Gₓ must be real."
        @assert length(Gₓ) == x_dims "Gₓ parameters have the wrong dimension."
        @assert typeof(Hₓ) <: Array "Hₓ must be an array."
        @assert all(imag(Hₓ) .== 0) "Hₓ must be real."
        @assert length(Hₓ) == x_dims "Hₓ parameters have the wrong dimension."
        @assert length(Bᵤ) == x_dims "Bᵤ parameters have the wrong dimension."

        rhs! = function (dz, z, edges, p, t)
            i = total_current(edges) + Y_n * (z[1] + z[2] * 1im)
            u = z[1] + z[2] * im
            x = [z[j] for j in 3:lastindex(z)]
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
    end

    ODEVertex(rhs!, dim, mass_matrix, sym)
end

function symbolsof(nf::NormalForm)
    @assert typeof(nf.x_dims) <: Int "The dimension of x must be an integer!"
    @assert nf.x_dims >= 0 "The dimension of x cannot be negative!"

    symbols = [:u_r, :u_i]
    append!(symbols, [Symbol("x_$i") for i in 1:x_dims])
    return symbols
end

function dimension(nf::NormalForm)
    @assert typeof(nf.x_dims) <: Int "The dimension of x must be an integer!"
    @assert nf.x_dims >= 0 "The dimension of x cannot be negative!"

    return 2 + nf.x_dims
end

export NormalForm