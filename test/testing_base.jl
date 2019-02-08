begin
        using Test
        using PowerDynBase
        # using SymEngine
        using SymPy
        Core.eval(PowerDynBase, :(using SymPy))
        Core.eval(PowerDynBase, :(pi = PI))
        using MacroTools
        using Test
        using Random
        random_seed = 1234
        @show random_seed
        Random.seed!(random_seed)
        using LinearAlgebra
end

begin
    # within construct_network_rhs, the array of real values (necessary for DifferentialEqautions.jl)
    # is converted to an array of complex values, mostly for the complex voltage.
    # As this is a small hack, we need to prepare here that SymPy Arrays can do the same.
    # In the following, we define something that emulates a reinterpreted vector.
    # Then we overload complex_view to use this emulated version instead
    # of the real reinterpret as done in GridDynamics.jl.
    using Base.Iterators: flatten, partition
    import Base: setindex!, size, getindex, IndexStyle
    import PowerDynBase: complexview

    struct ComplexTestVector{T} <: AbstractVector{T}
        vec::AbstractVector{T}
        function ComplexTestVector(vec::AbstractVector{T}) where T
            @assert length(vec) % 2 == 0 "need an even length to interpret it as an array of complex numbers"
            new{T}(vec)
        end
    end
    realindex(i::Integer) = 2*(i-1)+1
    imagindex(i::Integer) = realindex(i)+1
    function Base.setindex!(A::ComplexTestVector{T}, v, I::Vararg{Int}) where {T}
        (i, ) = I
        A.vec[[realindex(i), imagindex(i)]] = [real(v), imag(v)]
     end
    function getindex(v::ComplexTestVector{T}, i::Integer) where T
        i_real = 2*(i-1)+1
        i_imag = i_real + 1
        v.vec[i_real]+ v.vec[i_imag]*im
    end
    size(v::ComplexTestVector{T}) where T = (div(length(v.vec), 2),)
    Base.IndexStyle(::Type{<:ComplexTestVector{T}}) where T = IndexLinear()
    function complexview(vec::Array{Sym,1}, begin_index, num_complex_numbers)
        begin_index_original = 2*(begin_index-1)+1
        vec_view = view(vec, begin_index_original:begin_index_original + 2*(num_complex_numbers-1)+1)
        ComplexTestVector(vec_view)
    end

    # and two helper functions
    "Split an array of (complex) symbols into their real and complex parts."
    splitcomplexsymbols(vec::AbstractVector{Sym}) =
        collect(Base.Iterators.flatten([func(el) for func in (real, imag), el in vec]))
    "Merge an array of (real) symbols into (half the number of) complex symbols."
    mergecomplexsymbols(vec::AbstractVector{Sym}) =
        [vec[i]+im*vec[i+1] for i in 1:2:length(vec)]
    "Create a matrix of SymPy symbols."
    SymbolMatrix(var::AbstractString, n, m; kwargs...) = [symbols("$(var)_$i$j"; kwargs...) for i in 1:n, j in 1:m]
    SymbolMatrix(var::AbstractString, n; kwargs...) = SymbolMatrix(var, n, n; kwargs...)
end
