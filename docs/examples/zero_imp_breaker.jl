#=
# [Zero-Impedance Circuit Breaker](@id zero_imp_breaker)

This tutorial can be downloaded as a normal Julia script [here](@__NAME__.jl). #md

This example demonstrates how to model a **zero-impedance circuit breaker** in PowerDynamics.jl.
It is primarily a **pedagogical example**, the system presented in this tutorial itself
is not meant as a realistic physical simulation scenario, but rather as a teaching tool for
understanding how to handle such components in PowerDynamics.

The example features a simple two-generator, two-load system with a breaker connecting two
junction buses. We'll open the breaker at a preset time and then automatically reclose it
when the voltage angles across the breaker are synchronized.

=#

#=
## Imports

We use the standard PowerDynamics.jl stack along with ModelingToolkit for component definitions
and OrdinaryDiffEq for solving the resulting differential-algebraic equations.
=#

using PowerDynamics
using PowerDynamics.Library
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as Dt
using OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqNonlinearSolve
using CairoMakie

#=
## System Topology

The system consists of two generators (G1, G2) supplying two loads (L1, L2) through
junction buses (J1, J2) connected by a breaker:

```asciiart
G1 (~)─╂──╮J1 ╭───╂─▷ L1
         ╺┷━┯━┷╸
            │ Breaker
         ╺┯━┷━┯╸
G2 (~)─╂──╯J2 ╰───╂─▷ L2
```

The breaker can be opened and closed, allowing us to study transient behavior during
switching operations.


## Bus Definitions

We define two swing generators with different power specifications:
- **G1**: Slack bus (reference) with V=1.0 pu
- **G2**: PV bus with V=0.95 pu and P=1.0 pu

The system also has two constant admittance loads and two junction buses that serve
as connection points for the breaker.
=#

@named swing = Swing()
G1 = compile_bus(MTKBus(swing); name=:G1, vidx=1, pf=pfSlack(V=1.), swing₊M=0.1, swing₊D=0.1)
G2 = compile_bus(MTKBus(swing); name=:G2, vidx=2, pf=pfPV(V=0.95, P=1.0), swing₊D=0.1)

## Constant admittance loads
@named load = ConstantYLoad()
L1 = compile_bus(MTKBus(load); name=:L1, vidx=3, pf=pfPQ(P=-1.0, Q=-0.1))
L2 = compile_bus(MTKBus(load); name=:L2, vidx=4, pf=pfPQ(P=-1.0, Q=-0.1))

## Junction buses (no dynamics, just connection points)
J1 = compile_bus(MTKBus(); name=:J1, vidx=5)
J2 = compile_bus(MTKBus(); name=:J2, vidx=6)

nothing #hide #md

#=
## Line Model

We use a simple π-line model with low reactance and shunt capacitance. The lines connect
each generator to its respective junction bus, and each junction bus to its load.

!!! tip "Copy-Constructor of components"
    Here we use the "copy-constructor" of `EdgeModel` to create a new `EdgeModel` based
    on a previous one. This is useful to create multiple lines and is far more efficient
    than going through the full compilation process each time.
=#

@named pibranch = PiLine(;R=0.01, X=0.005, B_src=0.08, B_dst=0.08)
line_template = compile_line(MTKLine(pibranch))
lines = [
    EdgeModel(line_template; src=:G1, dst=:J1)
    EdgeModel(line_template; src=:J1, dst=:L1)
    EdgeModel(line_template; src=:G2, dst=:J2)
    EdgeModel(line_template; src=:J2, dst=:L2)
]

nothing #hide #md

#=
## Breaker Component Definition

The breaker is a switchable component that can connect or disconnect two buses. When closed,
it enforces zero voltage difference across its terminals (zero impedance). When open, it
allows no current to flow.

The breaker has a parameter `closed` that controls its state:
- `closed = 1`: Breaker is closed (zero impedance, voltages equal)
- `closed = 0`: Breaker is open (zero current)

We use `ifelse` to switch between these two operating modes. When closed, the equations
enforce `u_dst = u_src` with the current as an implicit output. When open, the equations
enforce `i = 0`.

!!! note "Usage of `implicit_output`""
    [`implicit_output`](@ref NetworkDynamics.implicit_output) evaluates to zero.
    Including it is just a trick to convice MTK that the current variables `i_r` and `i_i`
    are in some sense part of this constraint, because by chaging the current the solver
    can satisfy the voltage equality.
    This is necessary, because MTK does not know about the explicit feedback loop `u_src = f(i_src)`.
=#

@mtkmodel Breaker begin
    @parameters begin
        closed=1, [description="Breaker closed (1) or open (0)"]
    end
    @components begin
        src = Terminal()
        dst = Terminal()
    end
    @variables begin
        i_r(t), [guess=0, description="Current real part through breaker"]
        i_i(t), [guess=0, description="Current imaginary part through breaker"]
    end
    @equations begin
        ## When closed: enforce u_dst = u_src (with current as implicit output)
        ## When open: enforce i = 0
        0 ~ ifelse(closed == 1, dst.u_r - src.u_r + implicit_output(i_r), i_r)
        0 ~ ifelse(closed == 1, dst.u_i - src.u_i + implicit_output(i_i), i_i)
        ## Current flow equations (standard for all edge components)
        dst.i_r ~ i_r
        dst.i_i ~ i_i
        src.i_r ~ -dst.i_r
        src.i_i ~ -dst.i_i
    end
end

breaker_mod = Breaker(name=:breaker)
breaker = compile_line(MTKLine(breaker_mod), src=:J1, dst=:J2, name=:breaker)

#=
## Network Assembly

We assemble the full network with all buses, lines, and the breaker.
=#
nw = Network([G1, G2, L1, L2, J1, J2], [lines... , breaker]; warn_order=false)

#=
## Power Flow Initialization

We initialize the system from a power flow solution. The breaker starts in the closed state,
so the two sides of the system are initially connected and synchronized.
=#

s0 = initialize_from_pf!(nw)

nothing #hide #md

#=
## Breaker Control Callbacks

We define two callbacks to control the breaker:

1. **Opening callback**: A preset-time callback that opens the breaker at t=0.1s
2. **Closing callback**: A continuous callback that automatically recloses the breaker
   when the voltage angles across the breaker are synchronized (angle difference crosses zero)

The `toggle_breaker` affect function switches the breaker state by flipping the `closed` parameter.

Once the breaker is opend, we expect power imbalances and thus frequency
deviations in both subsystem 1 and subsystem 2.
At some point, the voltage angles on both sides will align agin, triggering the reclosing
callback. This is obviously not a realistic breaker control strategy, but it serves to
demostrate both opening and closing of the breaker in a single simulation.
=#

toggle_breaker = ComponentAffect([],[:breaker₊closed]) do u, p, ctx
    current = p[:breaker₊closed]
    next = current == 1 ? 0 : 1
    println("Toggling breaker state from ", Int(current), " to ", Int(next), " at t=", ctx.t)
    p[:breaker₊closed] = next
end

## Open the breaker at t=0.1s
open_breaker = PresetTimeComponentCallback(0.1, toggle_breaker)

## Close the breaker when voltage angles are synchronized
close_cond = ComponentCondition([:src₊u_r, :src₊u_i, :dst₊u_r, :dst₊u_i], [:breaker₊closed]) do u, p, ctx
    ## Don't trigger if already/still closed
    p[:breaker₊closed] == 1 && return Inf

    ## Calculate voltage angles on both sides
    src_arg = atan(u[:src₊u_i], u[:src₊u_r])
    dst_arg = atan(u[:dst₊u_i], u[:dst₊u_r])

    ## Return angle difference (callback triggers when this crosses zero)
    src_arg - dst_arg
end
close_breaker = ContinuousComponentCallback(close_cond, toggle_breaker)

nothing #hide #md

#=
!!! note "Angle Wrapping"
    Be carfull: due to the wrapping behavior of `atan`, the angle difference can jump.
    Therefore, this simple closing conditions won't work in many scenarios. Implementing a
    robust synchronization detection algorithm is beyond the scope of this tutorial.

## Simulation

We create an ODE problem and solve it using the Rodas5P solver (a stiff DAE solver).
The callbacks are attached to the breaker component using the `add_comp_cb` keyword argument.
=#

prob = ODEProblem(nw, s0, (0,2); add_comp_cb=[EIndex(:breaker) => (open_breaker, close_breaker)])
sol = solve(prob, Rodas5P())

nothing #hide #md

#=
## Results Visualization

We see that the voltages of both junction buses diverge after the breaker is opened at t=0.1s.
When the breaker recloses (when the voltage angles synchronize), the voltages realign and stay aligned.
=#

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
