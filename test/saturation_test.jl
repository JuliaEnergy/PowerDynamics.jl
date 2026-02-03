using PowerDynamics
using PowerDynamics.Library
PowerDynamics.load_pdtesting()
using Main.PowerDynamicsTesting
using ModelingToolkit
using NetworkDynamics
using Graphs
using OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqNonlinearSolve
using ModelingToolkit: ModelingToolkit as MTK
using ModelingToolkit: t_nounits as t, D_nounits as Dt
# using CairoMakie
using Test
using ScopedValues
using DiffEqCallbacks: PresetTimeCallback
using SciMLBase: auto_dt_reset!


@testset "Test configuration stuff" begin
    # Test SaturationConfiguration
    set_saturation_config!(SaturationConfig(:rhs_soft, regularization=1e-6))
    @test get_saturation_config().method == :rhs_soft
    set_saturation_config!(:callback)
    @test get_saturation_config().method == :callback

    with_saturation_config(:rhs_soft) do
       @test get_saturation_config().method == :rhs_soft
    end
    ScopedValues.with(SaturationConfiguration => SaturationConfig(:rhs_hard)) do
       @test get_saturation_config().method == :rhs_hard
    end
    @test get_saturation_config().method == :callback

    # Test CallbackVerbosity
    set_callback_verbosity!(false)
    @test get_callback_verbosity() == false
    set_callback_verbosity!(true)
    @test get_callback_verbosity() == true

    with_callback_verbosity(false) do
        @test get_callback_verbosity() == false
    end
    ScopedValues.with(CallbackVerbosity => false) do
        @test get_callback_verbosity() == false
    end
    @test get_callback_verbosity() == true
end

@testset "clamp function tests" begin
    using PowerDynamics.Library: _soft_clamped_rhs, _hard_clamped_rhs

    @testset "_hard_clamped_rhs" begin
        # Test within range - should pass through unchanged
        @test _hard_clamped_rhs(1.0, 0.5, 0.0, 1.0) == 1.0
        @test _hard_clamped_rhs(-1.0, 0.5, 0.0, 1.0) == -1.0

        # Test at upper limit with positive forcing - should clamp to zero
        @test _hard_clamped_rhs(1.0, 1.0, 0.0, 1.0) == 0.0
        @test _hard_clamped_rhs(0.5, 1.0, 0.0, 1.0) == 0.0

        # Test at upper limit with negative forcing - should pass through
        @test _hard_clamped_rhs(-1.0, 1.0, 0.0, 1.0) == -1.0
        @test _hard_clamped_rhs(-0.5, 1.0, 0.0, 1.0) == -0.5

        # Test at lower limit with negative forcing - should clamp to zero
        @test _hard_clamped_rhs(-1.0, 0.0, 0.0, 1.0) == 0.0
        @test _hard_clamped_rhs(-0.5, 0.0, 0.0, 1.0) == 0.0

        # Test at lower limit with positive forcing - should pass through
        @test _hard_clamped_rhs(1.0, 0.0, 0.0, 1.0) == 1.0
        @test _hard_clamped_rhs(0.5, 0.0, 0.0, 1.0) == 0.5

        # Test beyond limits
        @test _hard_clamped_rhs(2.0, 1.5, 0.0, 1.0) == 0.0  # above max, positive forcing
        @test _hard_clamped_rhs(-1.0, 1.5, 0.0, 1.0) == -1.0  # above max, negative forcing
        @test _hard_clamped_rhs(-2.0, -0.5, 0.0, 1.0) == 0.0  # below min, negative forcing
        @test _hard_clamped_rhs(1.0, -0.5, 0.0, 1.0) == 1.0  # below min, positive forcing
    end

    @testset "_soft_clamped_rhs" begin
        ε = 0.01

        # Test well within range - should be approximately equal to u
        @test _soft_clamped_rhs(1.0, 0.5, 0.0, 1.0, ε) ≈ 1.0
        @test _soft_clamped_rhs(-1.0, 0.5, 0.0, 1.0, ε) ≈ -1.0

        # Test at center of range - maximum pass-through
        u_mid = _soft_clamped_rhs(2.0, 0.5, 0.0, 1.0, ε)
        @test u_mid ≈ 2.0 atol=0.01

        # Test near upper limit - should be attenuated
        u_near_max = _soft_clamped_rhs(1.0, 0.95, 0.0, 1.0, ε)
        @test u_near_max < 1.0
        @test u_near_max > 0.0

        # Test at upper limit - should be significantly attenuated
        u_at_max = _soft_clamped_rhs(1.0, 1.0, 0.0, 1.0, ε)
        @test u_at_max <= 0.5
        @test u_at_max > 0.0

        # Test beyond upper limit - should be heavily attenuated
        u_past_max = _soft_clamped_rhs(1.0, 1.2, 0.0, 1.0, ε)
        @test u_past_max < u_at_max
        @test u_past_max ≈ 0.0 atol=0.01

        # Test near lower limit - should be attenuated
        u_near_min = _soft_clamped_rhs(-1.0, 0.05, 0.0, 1.0, ε)
        @test abs(u_near_min) < 1.0
        @test abs(u_near_min) > 0.0

        # Test at lower limit - should be significantly attenuated
        u_at_min = _soft_clamped_rhs(-1.0, 0.0, 0.0, 1.0, ε)
        @test abs(u_at_min) <= 0.5
        @test abs(u_at_min) > 0.0

        # Test beyond lower limit - should be heavily attenuated
        u_past_min = _soft_clamped_rhs(-1.0, -0.2, 0.0, 1.0, ε)
        @test abs(u_past_min) < abs(u_at_min)
        @test u_past_min ≈ 0.0 atol=0.01

        # Test symmetry - clamping should be symmetric for symmetric limits
        @test _soft_clamped_rhs(1.0, 0.5, -1.0, 1.0, ε) ≈ _soft_clamped_rhs(-1.0, -0.5, -1.0, 1.0, ε) * -1.0
    end
end

# Integration test for saturation limiters
using PowerDynamics.Library: SimpleLagLim, LimIntegrator, SaturationConfig, SaturationConfiguration
using LinearAlgebra: norm

"""
SaturationTestInjector - Test node with both SimpleLagLim and LimIntegrator

This injector node accepts voltage as input and outputs current (loopback pattern).
It contains both a lag limiter and an integrator limiter to test saturation behavior.
"""
@mtkmodel SaturationTestInjector begin
    @structural_parameters begin
        min = -0.5
        max = 1.5
    end
    @components begin
        terminal = Terminal()
        # Lag limiter with asymmetric limits
        # Uses SaturationConfiguration scoped value
        lag = SimpleLagLim(;
            K = 1.0,
            T = 1.0,
            outMin = min,
            outMax = max,
            guess = 0.5
        )

        # Integrator with different limits
        integrator = LimIntegrator(;
            K = 1.0,
            T = 1.0,
            outMin = min,
            outMax = max,
            guess = 0.0,
        )
    end
    @equations begin
        # Drive limiters with voltage deviation from 1.0 pu
        lag.in ~ terminal.u_r
        integrator.in ~ terminal.u_i

        # Output small dummy current proportional to limiter outputs
        terminal.i_r ~ lag.out
        terminal.i_i ~ integrator.out
    end
end

function build_test_network(config=get_saturation_config(); min=-0.5, max=1.0)
    hub_vertex = pfSlack(; u_r=0, u_i=0, name=:slack_hub)

    sat_vertex = with(SaturationConfiguration => config) do
        @named injector = SaturationTestInjector(; min, max)

        # Create vertex model: voltage in, current out
        VertexModel(
            MTKBus(injector),
            [:busbar₊u_r, :busbar₊u_i],  # inputs (voltage)
            [:busbar₊i_r, :busbar₊i_i],  # outputs (current)
            ff_to_constraint=false,  # required for injector
            name = Symbol(:sat_, config.method)
        )
    end
    vertices = [
        hub_vertex,
        sat_vertex
    ]

    edge = LoopbackConnection(
        potential = [:busbar₊u_r, :busbar₊u_i],
        flow = [:terminal₊i_r, :terminal₊i_i],
        src = sat_vertex.name,
        dst = :slack_hub,
        name = :loopback_callback
    )

    Network(vertices, edge; warn_order=false)
end

function sim_network(nw)
    tstops = [1.0, 4.0, 7.0]

    affect = function(integrator)
        u = NWState(integrator)
        t = integrator.t

        if t == 1.0
            u.v[1, :slack₊u_r] = 2
            u.v[1, :slack₊u_i] = 2
        elseif t == 4.0
            u.v[1, :slack₊u_r] = -1
            u.v[1, :slack₊u_i] = -1
        elseif t == 7.0
            u.v[1, :slack₊u_r] = 0.1
            u.v[1, :slack₊u_i] = -1
        end

        auto_dt_reset!(integrator)
        save_parameters!(integrator)
    end

    cb = PresetTimeCallback(tstops, affect)
    s0 = NWState(nw, ufill=0)
    prob = ODEProblem(nw, s0, (0.0, 10.0); add_nw_cb=cb)
    sol = solve(prob, Rodas5P());
end

@testset "Different saturation limiters" begin
    configs = Dict(
        :callback => SaturationConfig(:callback),
        :complementary => SaturationConfig(:complementary; regularization=1e-8),
        :rhs_hard => SaturationConfig(:rhs_hard;),
        :rhs_soft => SaturationConfig(:rhs_soft; regularization=1e-8),
    )

    sols = Dict(k => sim_network(build_test_network(v)) for (k, v) in configs);

    ts = range(sols[:callback].t[begin], sols[:callback].t[end], length=1000)
    int_outs = Dict(k => sol(ts, idxs=VIndex(2, :injector₊integrator₊out)).u for (k, sol) in sols)
    @test maximum(map(row -> abs(-(extrema(row)...)), eachrow(hcat(values(int_outs)...)))) < 0.01
    lag_outs = Dict(k => sol(ts, idxs=VIndex(2, :injector₊lag₊out)).u for (k, sol) in sols)
    @test maximum(map(row -> abs(-(extrema(row)...)), eachrow(hcat(values(lag_outs)...)))) < 0.03

    #=
    fig = let
        fig = Figure(size=(600,800))
        for (row, (type, sol)) in enumerate(sols)
            Label(fig[row, 1], string(type); rotation=π/2, tellheight=false)
            # saturated integrator
            if row == length(sols)
                ax = Axis(fig[row, 2]; xlabel="Time (s)")
            elseif row == 1
                ax = Axis(fig[row, 2]; title="Saturated Integrator")
            else
                ax = Axis(fig[row, 2])
            end
            lines!(ax, sol; idxs=VIndex(2, :injector₊integrator₊forcing), color=:lightgray)
            lines!(ax, sol, idxs=VIndex(2, :injector₊integrator₊out))
            lines!(ax, sol; idxs=VIndex(2, :injector₊integrator₊min), color=:red, linestyle=:dash)
            lines!(ax, sol; idxs=VIndex(2, :injector₊integrator₊max), color=:red, linestyle=:dash)
            ylims!(ax, -2.2, 2.2)
            # ylims!(ax, -0.102, -0.099)
            # saturated lag
            if row == length(sols)
                ax = Axis(fig[row, 3]; xlabel="Time (s)")
            elseif row == 1
                ax = Axis(fig[row, 3]; title="Saturated PT1 Lag")
            else
                ax = Axis(fig[row, 3])
            end
            lines!(ax, sol; idxs=VIndex(2, :injector₊lag₊forcing), color=:lightgray)
            lines!(ax, sol, idxs=VIndex(2, :injector₊lag₊out))
            lines!(ax, sol; idxs=VIndex(2, :injector₊lag₊min), color=:red, linestyle=:dash)
            lines!(ax, sol; idxs=VIndex(2, :injector₊lag₊max), color=:red, linestyle=:dash)
            ylims!(ax, -2.2, 2.2)
            # ylims!(ax, -0.102, -0.099)
        end
        fig
    end
    =#
end

@testset "initialiation" begin
    configs = Dict(
        :callback => SaturationConfig(:callback),
        :complementary => SaturationConfig(:complementary; regularization=1e-8),
        :rhs_hard => SaturationConfig(:rhs_hard;),
        :rhs_soft => SaturationConfig(:rhs_soft; regularization=1e-8),
    )
    vms = Dict(k => build_test_network(v)[VIndex(2)] for (k, v) in configs)

    function test_init(vm, int, lag; tol=1e-10, verbose=false)
        set_default!(vm, :busbar₊u_i, int.in)
        set_default!(vm, :busbar₊u_r, lag.in)
        set_default!(vm, :busbar₊i_i, -int.out)
        set_default!(vm, :busbar₊i_r, -lag.out)

        lagt = false
        intt = false
        err = nothing
        try
            verbose && printstyled("\nInitializing vertex $(vm.name)... \n", color=:blue)
            cs = initialize_component(vm; tol, verbose)
            if verbose
                printstyled("Residuals:\n", color=:blue)
                init_residual(vm, cs; verbose=true)
            end
            lagt = isapprox(cs[:injector₊lag₊x], lag.out; atol=tol)
            intt = isapprox(cs[:injector₊integrator₊x], int.out; atol=tol)
        catch e
            err = e
        end

        if !isnothing(err) || !lagt || !intt
            initialize_component!(vm, tol=Inf)
            printstyled("\nInitialization failed for test_init at vertex $(vm.name)\n", color=:blue)
            dump_initial_state(vm)
            printstyled("\nGot following Residuals:\n", color=:blue)
            init_residual(vm; verbose=true)
            isnothing(err) || throw(err)
        end
        lagt && intt
    end

    # both in bounds
    int = (; in=0.0, out=0.3, res=0.3)
    lag = (; in=0.0, out=0.0, res=0.0)
    @test test_init(vms[:callback], int, lag, verbose=false)
    @test test_init(vms[:complementary], int, lag, tol=1e-7)
    @test test_init(vms[:rhs_hard], int, lag)
    @test test_init(vms[:rhs_soft], int, lag) # soft rhs might lead to reduce accuracy

    # lag above bound
    int = (; in=0.0, out=0.0, res=0.0)
    lag = (; in=1.5, out=1.0, res=1.0)
    @test_broken test_init(vms[:callback], int, lag, verbose=false) # callbacks do not respect bounds
    @test_broken test_init(vms[:complementary], int, lag, tol=1e-7)
    @test test_init(vms[:rhs_hard], int, lag)
    @test test_init(vms[:rhs_soft], int, lag, tol=1e-7) # soft rhs might lead to reduce accuracy

    # int above bound
    int = (; in=1, out=1.0, res=1.0)
    lag = (; in=0, out=0.0, res=0.0)
    @test_broken test_init(vms[:callback], int, lag, verbose=false) # callbacks do not respect bounds
    @test_broken test_init(vms[:complementary], int, lag, tol=1e-7)
    @test test_init(vms[:rhs_hard], int, lag)
    @test test_init(vms[:rhs_soft], int, lag, tol=1e-6) # soft rhs might lead to reduce accuracy

    # both below bound
    int = (; in=-1.0, out=-0.5, res=-0.5)
    lag = (; in=-1.0, out=-0.5, res=-0.5)
    @test_broken test_init(vms[:callback], int, lag, verbose=false) # callbacks do not respect bounds
    @test_broken test_init(vms[:complementary], int, lag, tol=1e-7)
    @test test_init(vms[:rhs_hard], int, lag)
    @test test_init(vms[:rhs_soft], int, lag, tol=1e-6) # soft rhs might lead to reduce accuracy
end
