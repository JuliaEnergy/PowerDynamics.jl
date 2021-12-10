function smoketest_rhs(edge; int_x=[], int_dx=[])
    vs, vd, u, di, i, di_r, i_r = rand(ComplexF64, 7)
    t = rand(0.1:10.0)

    dx = [real(di), imag(di), real(di_r), imag(di_r)]
    x = [real(i), imag(i), real(i_r), imag(i_r)]
    append!(x, int_x)
    append!(dx, int_dx)
    v_s = [real(vs), imag(vs)]
    v_d = [real(vd), imag(vd)]
    p = []

    edge.f(dx, x, v_s, v_d, p, t)
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
