using PowerDynamics
using DifferentialEquations
using Plots

mutable struct Pars
    τ
end


nodes = [PQAlgebraic(P=-1.0,Q=-im*0.5),
        #SwingEqLVS(H=1., P=-1.0, D=1, Ω=50, Γ=20, V=1),
        SwingEqLVS(H=1., P=1.0, D=1, Ω=50, Γ=20, V=1)]
lines = [StaticLine(Y=-50*im, from=1, to=2)]
pg = PowerGrid(nodes, lines)

f! = rhs(pg)
const out = zeros(systemsize(pg)) # Define a cache variable
h(out, p, t) = (out.=1.0)

function wrapper!(dx, x, h, p, t)
    h(out, p, t-p.τ) # updates out to be the correct history function
    f!.f(dx, x, p, t)
    dx[end] += out[end] # example: additive delay in last index
    return nothing
end

delay_pg = DDEFunction(wrapper!, mass_matrix=f!.mass_matrix, syms=f!.syms)

u0 = 0.1rand(systemsize(pg))
#du0 = similar(u0)
#wrapper!(du0, u0, h, p, 0)

p = Pars(0.1)

prob = DDEProblem(delay_pg, u0, h, (0., 5.), p; constant_lags=p.τ)

sol = solve(prob, MethodOfSteps(Rosenbrock23(autodiff=false)))

plot(sol)
