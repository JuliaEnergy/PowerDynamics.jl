using PowerDynamics
using PowerDynamics.Library
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as Dt
using OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqNonlinearSolve
using CairoMakie
#=
G1 (~)─╂──╮J1 ╭───╂─▷ L2
         ╺┷━┯━┷╸
            │ Breaker
         ╺┯━┷━┯╸
G2 (~)─╂──╯J2 ╰───╂─▷ L2
=#

@named swing = Swing()
G1 = compile_bus(MTKBus(swing); name=:G1, vidx=1, pf=pfSlack(V=1.), swing₊M=0.1, swing₊D=0.1)
G2 = compile_bus(MTKBus(swing); name=:G2, vidx=2, pf=pfPV(V=0.95, P=1.0), swing₊D=0.1)

@named load = ConstantYLoad()
L1 = compile_bus(MTKBus(load); name=:L1, vidx=3, pf=pfPQ(P=-1.0, Q=-0.1))
L2 = compile_bus(MTKBus(load); name=:L2, vidx=4, pf=pfPQ(P=-1.0, Q=-0.1))

J1 = compile_bus(MTKBus(); name=:J1, vidx=5)
J2 = compile_bus(MTKBus(); name=:J2, vidx=6)


# @named pibranch = PiLine(;R=0.01, X=0.085, B_src=0.088, B_dst=0.088)
# @named pibranch = PiLine(;R=0.01, X=0.085, B_src=0., B_dst=0.)
@named pibranch = PiLine(;R=0.01, X=0.005, B_src=0., B_dst=0.)
line_template = compile_line(MTKLine(pibranch))
lines = [
    EdgeModel(line_template; src=:G1, dst=:J1)
    EdgeModel(line_template; src=:J1, dst=:L1)
    EdgeModel(line_template; src=:G2, dst=:J2)
    EdgeModel(line_template; src=:J2, dst=:L2)
]

@mtkmodel Breaker begin
    @parameters begin
        closed=1, [description="Breaker closed (1) or open (0)"]
    end
    @components begin
        src = Terminal()
        dst = Terminal()
    end
    @variables begin
        i_r(t)=0, [description="Current real part through breaker"]
        i_i(t)=0, [description="Current imaginary part through breaker"]
    end
    @equations begin
        # 0 ~ dst.u_r - src.u_r + implicit_output(i_r)
        # 0 ~ dst.u_i - src.u_i + implicit_output(i_i)
        0 ~ ifelse(closed == 1, dst.u_r - src.u_r + implicit_output(i_r), i_r)
        0 ~ ifelse(closed == 1, dst.u_i - src.u_i + implicit_output(i_i), i_i)
        dst.i_r ~ i_r
        dst.i_i ~ i_i
        src.i_r ~ -dst.i_r
        src.i_i ~ -dst.i_i
    end
end
breaker_mod = Breaker(name=:breaker)
breaker = compile_line(MTKLine(breaker_mod), src=:J1, dst=:J2, name=:breaker)
copyi = @pfinitformula begin
    :breaker₊i_r = @pf(:breaker₊i_r)
    :breaker₊i_i = @pf(:breaker₊i_i)
end
set_pfinitformula!(breaker, copyi)

nw = Network([G1, G2, L1, L2, J1, J2], [lines... , breaker]; warn_order=false)

s0 = initialize_from_pf!(nw; subverbose=true)

toggle_breaker = ComponentAffect([],[:breaker₊closed]) do u, p, ctx
    current = p[:breaker₊closed]
    next = current == 1 ? 0 : 1
    println("Toggling breaker state from ", current, " to ", next, " at t=", ctx.t)
    p[:breaker₊closed] = next
end
open_breaker = PresetTimeComponentCallback(0.1, toggle_breaker)

close_cond = ComponentCondition([:src₊u_r, :src₊u_i, :dst₊u_r, :dst₊u_i], [:breaker₊closed]) do u, p, ctx
    p[:breaker₊closed] == 1 && return Inf
    # ctx.t < 0.2 && return Inf

    src_arg = atan(u[:src₊u_i], u[:src₊u_r])
    dst_arg = atan(u[:dst₊u_i], u[:dst₊u_r])
    println(src_arg - dst_arg)

    src_arg - dst_arg
end
close_breaker = ContinuousComponentCallback(close_cond, toggle_breaker)

prob = ODEProblem(nw, s0, (0,5); add_comp_cb=[EIndex(:breaker) => (open_breaker, close_breaker)])
sol = solve(prob, Rodas5P());
break

let
    fig = Figure(size=(600,600))
    ax = Axis(fig[1, 1]; title="Voltage Magnitude", xlabel="Time [s]", ylabel="Voltage [pu]")
    ts = refine_timeseries(sol.t)
    lines!(ax, ts, sol(ts; idxs=VIndex(:J1, :busbar₊u_mag)).u, label="J1")
    lines!(ax, ts, sol(ts; idxs=VIndex(:J2, :busbar₊u_mag)).u, label="J2")
    axislegend(ax)
    ax = Axis(fig[2, 1]; title="Voltage Angle", xlabel="Time [s]", ylabel="Angle [rad]")
    lines!(ax, ts, sol(ts; idxs=VIndex(:J1, :busbar₊u_arg)).u, label="J1")
    lines!(ax, ts, sol(ts; idxs=VIndex(:J2, :busbar₊u_arg)).u, label="J2")
    axislegend(ax)
    ax = Axis(fig[3, 1]; title="Swing Frequency", xlabel="Time [s]", ylabel="Frequency [pu]")
    lines!(ax, ts, sol(ts; idxs=VIndex(:G1, :swing₊ω)).u, label="G1")
    lines!(ax, ts, sol(ts; idxs=VIndex(:G2, :swing₊ω)).u, label="G2")
    axislegend(ax)
    fig
end
