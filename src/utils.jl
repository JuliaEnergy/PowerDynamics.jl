using Chairmarks
using CUDA

export compare_execution_styles
function compare_execution_styles(prob)
    styles = [
            # KAExecution{true}(),
              # KAExecution{false}(),
              SequentialExecution{true}(),
              SequentialExecution{false}(),
              PolyesterExecution{true}(),
              PolyesterExecution{false}(),
              # ThreadedExecution{true}(),
              # ThreadedExecution{false}()
              ]

    aggregators = [
        # KAAggregator,
        SequentialAggregator,
        PolyesterAggregator,
        # ThreadedAggregator,
        # SparseAggregator
    ]

    @assert prob isa ODEProblem "test_execution_styles only works for ODEProblems"

    u = copy(prob.u0)
    du = zeros(eltype(u), length(u))
    t = 0.0
    p = copy(prob.p)
    nw = prob.f.f
    nw(du, u, p, t)
    @assert u==prob.u0
    @assert p==prob.p

    exsaggs = [(ex, agg) for ex in styles for agg in aggregators]

    results = Dict()

    for (execution, aggregator) in exsaggs
        _nw = Network(nw; execution, aggregator=aggregator(nw.layer.aggregator.f))
        _du = zeros(eltype(u), length(u))
        try
            println("$execution\t $aggregator")
            b = @b $_nw($_du, $u, $p, $t) seconds=1
            show(b)
            results[(execution,aggregator)] = b
        catch e
            # XXX: fix for https://github.com/JuliaLang/julia/issues/55075
            if e isa MethodError && e.f == Base.elsize
                continue
            end
            println("Error in $execution with $aggregator: $e")
            @assert false
            continue
        end
        issame = isapprox(_du, du; atol=1e-10)
        if !issame
            println("$execution with $aggregator lead to different results: extrema(Δ) = $(extrema(_du - du))")
        end
        @assert issame
    end

    if CUDA.functional()
        to = CuArray
        u_d = adapt(to, u)
        p_d = adapt(to, p)

        for (execution, aggregator) in exsaggs
            (iscudacompatible(execution) && iscudacompatible(aggregator)) || continue

            _nw = Network(nw; execution, aggregator=aggregator(nw.layer.aggregator.f))
            _nw_d = adapt(to, _nw)
            _du_d = adapt(to, zeros(eltype(u), length(u)))


            println("$execution\t $aggregator")
            b = @b $_nw_d($_du_d, $u_d, $p_d, $t)
            show(b)
            # b =@b_nw_d(_du_d, u_d, p_d, t)
            issame = isapprox(Vector(_du_d), du; atol=1e-10)
            if !issame
                println("CUDA execution lead to different results: extrema(Δ) = $(extrema(_du - du))")
            end
            @assert issame
        end
    end
    return results
end
