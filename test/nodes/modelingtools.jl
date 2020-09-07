using ModelingToolkit
using OrdinaryDiffEq
using Plots

t = Variable(:t; known = true)()
Ki = Variable(:Ki; known = true)()
Kp = Variable(:Kp; known = true)()
Kd = Variable(:Kd; known = true)()
x_ic = Variable(:x_ic; known = true)()
x = Variable(:x)(t)
u = Variable(:u)(t)
D = Differential(t)

de = ODESystem([D(x) ~ -Kp * x], t, [x], [Kp])
ode = ODEProblem(ODEFunction(de, [x], [Kp]), rand(1), (0., 10.), [Kp])

sol = solve(ode, Tsit5())
plot(sol)


function system(x, u; A, B)
        return D(x) ~ A * x + B * u
end

function observation(x; C)
        return C * x
end

function additive_error(y; y_ref)
        return y_ref - y
end

function P(e; Kp)
        return Kp * e
end

function PI(u, e; Kp, Ki, e_ic)
        return D(u) ~ Ki * (e - e_ic) + Kp
end

function closed_P(x)
        @variables y(t), e(t) #u(t)

        y = observation(x; C=1.)
        e = additive_error(y; y_ref=0.)
        u = P(e; Kp=Kp)
        # return type needs to be Vector{Equation}
        return [system(x, u; A=0., B=1.),]
end

function closed_PI(state)
        x, u = state
        @variables y(t), e(t) #u(t)

        y = observation(x; C=1.)
        e = additive_error(y; y_ref=0.)

        e_ic = additive_error(observation(x_ic; C=1.); y_ref=0.)
        # return type needs to be Vector{Equation}
        return [system(x, u; A=0., B=1.),
                PI(u, e; Kp=Kp, Ki=Ki, e_ic=e_ic)] # return type needs to be
end

function closed_PID(state)
        x, u = state
        @variables y(t), e(t) #u(t)

        A = 1.
        B = 1.
        C = 1.
        # return type needs to be Vector{Equation}
        return [D(x) ~ u
                D(u) ~ ((B*Kp*C-A) * u - B*Ki*C * (x - x_ic)) / (1-B*Kd*C) ] # return type needs to be
end


system_eqns = closed_P(x)
de = ODESystem(system_eqns, t, [x], [Kp])
f = ODEFunction(de)
ode_P = ODEProblem(f, rand(1), (0., 10.), ones(1))

sol = solve(ode_P, Tsit5())
plot(sol)

system_eqns = closed_PI([x, u])
de = ODESystem(system_eqns, t, [x, u], [Kp, Ki, x_ic])
f = ODEFunction(de)
x0 = rand(2)
ode_PI = ODEProblem(f, x0, (0., 10.), [1., 1., x0[1]])

sol = solve(ode_PI, Tsit5())
plot(sol)

system_eqns = closed_PID([x, u])
de = ODESystem(system_eqns, t, [x, u], [Kp, Ki, Kd, x_ic])
f = ODEFunction(de)
x0 = rand(2)
ode_PI = ODEProblem(f, x0, (0., 100.), [0.1, 0.1, 0.1, x0[1]])

sol = solve(ode_PI, Tsit5())
plot(sol)
