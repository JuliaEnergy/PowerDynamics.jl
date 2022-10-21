function smoketest_rhs(vertex; int_x=[], int_dx=[])
    du,u,i = rand(ComplexF64, 3)
    t = rand(0.1:10.0)

    dx = [real(du), imag(du)]
    x = [real(u), imag(u)]
    append!(x, int_x)
    append!(dx, int_dx)
    edges = [[real(i), imag(i), real(i), imag(i)]]
    p = nothing
    vertex.f(dx, x, edges, p, t) # first call for precompilation

    allocations = @allocated vertex.f(dx, x, edges, p, t)
    if allocations != 0
        @warn "RHS shows non-zero ($allocations) allocations!"
    end
end

function rand_real()
    rand(1.0:10.0)
end

function rand_real(n)
    rand(1.0:10.0, n)
end

function rand_positive()
    rand(1:10)
end

function rand_positive(n)
    rand(1:10,n)
end

function rand_negative()
    rand(-10:-1)
end

function rand_negative(n)
    rand(-10:-1,n)
end
