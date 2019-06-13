using SymPy: @syms

@syms du u i
@syms t real=true

function smoketest_rhs(vertex; int_x=[], int_dx=[])
    dx = [real(du), imag(du)]
    x = [real(u), imag(u)]
    append!(x, int_x)
    append!(dx, int_dx)
    e_s = []
    e_d = [[real(i), imag(i), real(i), imag(i)]]
    p = []
    vertex.f!(dx, x, e_s, e_d, p, t)
end
