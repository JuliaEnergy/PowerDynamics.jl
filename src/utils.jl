using Chairmarks
using CUDA
using NetworkDynamics.Adapt
using NetworkDynamics: iscudacompatible
import ArgCheck: @argcheck
import Symbolics: Symbolics, Symbolic, fixpoint_sub, simplify, isterm
import ModelingToolkit: getname, unwrap, rename

export get_PQV
function get_PQV(s::NWState)
    us = s[VIndex(:, :u_r)] .+ im.*s[VIndex(:, :u_i)]
    ubuf, aggbuf = NetworkDynamics.get_ustacked_buf(s)
    P = zeros(Float64, length(us))
    Q = zeros(Float64, length(us))
    for i in 1:nv(s.nw)
        aggr = s.nw.im.v_aggr[i]
        current = _tocomplex(@views aggbuf[aggr])
        S = us[i] * conj(current)
        P[i], Q[i] = _toreal(S)
    end
    V = abs.(us)
    P, Q, V
end

export get_VΘ
function get_VΘ(s::NWState)
    us = s[VIndex(:, :u_r)] .+ im.*s[VIndex(:, :u_i)]
    abs.(us), angle.(us)
end

export compare_execution_styles
function compare_execution_styles(prob)
    styles = [
              KAExecution{true}(),
              # KAExecution{false}(),
              SequentialExecution{true}(),
              # SequentialExecution{false}(),
              PolyesterExecution{true}(),
              PolyesterExecution{false}(),
              # ThreadedExecution{true}(),
              # ThreadedExecution{false}()
              ]

    aggregators = [
        KAAggregator,
        SequentialAggregator,
        PolyesterAggregator,
        # ThreadedAggregator,
        SparseAggregator
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

    println("CPU benchmarks")
    for (execution, aggregator) in exsaggs
        _nw = Network(nw; execution, aggregator=aggregator(nw.layer.aggregator.f))
        _du = zeros(eltype(u), length(u))
        try
            println("$execution\t $aggregator")
            b = @b $_nw($_du, $u, $p, $t) seconds=1
            # show(b)
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

    gpuresults = Dict()
    if CUDA.functional()
        println("GPU benchmarks")
        to = CuArray
        u_d = adapt(to, u)
        p_d = adapt(to, p)

        for (execution, aggregator) in exsaggs
            (iscudacompatible(execution) && iscudacompatible(aggregator)) || continue

            _nw = Network(nw; execution, aggregator=aggregator(nw.layer.aggregator.f))
            _nw_d = adapt(to, _nw)
            _du_d = adapt(to, zeros(eltype(u), length(u)))


            println("$execution\t $aggregator")
            b = @b $_nw_d($_du_d, $u_d, $p_d, $t) seconds=1
            # show(b)
            gpuresults[(execution,aggregator)] = b
            # b =@b_nw_d(_du_d, u_d, p_d, t)
            issame = isapprox(Vector(_du_d), du; atol=1e-10)
            if !issame
                println("CUDA execution lead to different results: extrema(Δ) = $(extrema(_du - du))")
            end
            @assert issame
        end
    end
    return results, gpuresults
end

export pinparameters

function pinparameters(sys::ODESystem, sub::Pair)
    pinparameters(sys, Dict(sub.first => sub.second))
end
function pinparameters(sys::ODESystem, subs::Dict)
    _dict = Dict(_resolve_to_namespaced_symbolic(sys, k) => v for (k,v) in subs)

    if any(isterm, values(_dict))
        @info "Parameter map for pinparameters cannot include any terms in the rhs. Try to resolve..."
        __dict = Dict()
        for (k,v) in _dict
            _fullname = string(getname(k))
            newname = replace(_fullname, r"^"*string(sys.name)*"₊" => "")
            __dict[rename(k, Symbol(newname))] = v
        end
        for (k,v) in _dict
            _dict[k] = simplify(fixpoint_sub(v, __dict))
            # _dict[k] = fixpoint_sub(v, __dict)
        end
        if any(Base.Fix2(isa, Symbolic), values(_dict))
            error("Could not resolve all dependent parameters in the parameter map.")
        end
    end

    _pinparameters(sys, _dict)
end

function _pinparameters(sys::ODESystem, _subs::Dict{<:Symbolic, <:Any})
    isempty(_subs) && return sys

    subs = Dict{Symbolic, Any}()
    for (k,v) in _subs
        _fullname = string(getname(k))
        contains(_fullname, r"^"*string(sys.name)*"₊") || continue
        newname = replace(_fullname, r"^"*string(sys.name)*"₊" => "")
        subs[rename(k, Symbol(newname))] = v
    end

    applicable = Dict{Symbolic, Any}()
    applicable_subs = intersect(keys(subs), parameters(sys))
    for k in applicable_subs
        applicable[k] = subs[k]
    end

    @info "$(sys.name): substitute following parameters $(keys(applicable))"
    if isempty(applicable)
        _eqs = unwrap(sys.eqs)
        _observed = unwrap(sys.observed)
        # _ps = sys.ps
        # _var_to_name = sys.var_to_name
        _defaults = sys.defaults
    else
        _eqs = [fixpoint_sub(eq, applicable) for eq in unwrap(sys.eqs)]

        # FIXME: check something drops the metadata
        # missingmetadata = check_metadata(_eqs)
        # if !isempty(missingmetadata)
        #     @warn "Metadata was dropped from the equations: $missingmetadata"
        #     _eqs = try_fix_metadata(_eqs, unwrap(sys.eqs))
        # end
        @info check_metadata(_eqs)
        fix_metadata!(_eqs, sys)


        _observed = [fixpoint_sub(eq, applicable) for eq in unwrap(sys.observed)]
        # missingmetadata = check_metadata(_observed)
        # if !isempty(missingmetadata)
        #     @warn "Metadata was dropped from the observed equations: $missingmetadata"
        #     _observed = try_fix_metadata(_observed, unwrap(sys.observed))
        # end
        fix_metadata!(_observed, sys)

        @argcheck unwrap(isempty(sys.ctrls)) "pin_parameters does not handle control variables"
        _defaults = copy(sys.defaults)
        for k in keys(applicable)
            delete!(_defaults, k)
        end
        @argcheck isnothing(sys.initializesystem) "pin_parameters does not handle initializesystem"
        @argcheck unwrap(isempty(sys.continuous_events)) "pin_parameters does not handle continuous events"
        @argcheck unwrap(isempty(sys.discrete_events)) "pin_parameters does not handle discrete events"
        @argcheck unwrap(isempty(sys.parameter_dependencies)) "pin_parameters does not parameter_dependencies"
        @argcheck isnothing(sys.discrete_subsystems) "pin_parameters does not handle discrete subsystems"
        @argcheck isnothing(sys.solved_unknowns) "pin_parameters does not handle solved unknowns"
    end

    _systems = [_pinparameters(s, subs) for s in sys.systems]

    newsys = ODESystem(_eqs, sys.iv;
            # controls = Num[],
            observed = _observed,
            systems = _systems,
            tspan = sys.tspan,
            name = sys.name,
            defaults = _defaults,
            guesses = sys.guesses,
            # initializesystem = nothing,
            # initialization_eqs = Equation[],
            schedule = sys.schedule,
            connector_type = sys.connector_type,
            preface = sys.preface,
            # continuous_events = nothing,
            # discrete_events = nothing,
            # parameter_dependencies = nothing,
            checks = true,
            metadata = sys.metadata,
            gui_metadata = sys.gui_metadata)
    @info "after" check_metadata(_eqs) check_metadata(newsys.eqs) check_metadata(full_equations(newsys))
    newsys
end
function _resolve_to_namespaced_symbolic(sys, var)
    ns = string(getname(sys))
    varname = string(getname(var))
    varname_nons = replace(varname, r"^"*ns*"₊" => "")
    parts = split(varname_nons, "₊")
    r = getproperty(sys, Symbol(parts[1]); namespace=true)
    for part in parts[2:end]
        r = getproperty(r, Symbol(part); namespace=true)
    end
    unwrap(r)
end

function check_metadata(exprs)
    nometadata = []
    for ex in exprs
        if ex isa Equation
            _check_metadata!(nometadata, ex.rhs)
            _check_metadata!(nometadata, ex.lhs)
        else
            _check_metadata!(nometadata, ex)
        end
    end
    return unique!(nometadata)
end
function _check_metadata!(nometadata, expr)
    vars = Symbolics.get_variables(expr)
    for v in vars
        isnothing(Symbolics.metadata(v)) && push!(nometadata, v)
    end
end

function fix_metadata!(invalid_eqs, sys)
    missingmetadata = check_metadata(invalid_eqs)
    if isempty(missingmetadata)
        return invalid_eqs
    end

    metadatasubs = Dict()
    allsyms = ModelingToolkit.all_symbols(sys)
    allnames = string.(ModelingToolkit.getname.(allsyms))
    for invalids in missingmetadata
        invalidname = getname(invalids)
        valid = if hasproperty(sys, getname(invalidname))
            getproperty(sys, getname(invalidname); namespace=false)
        else
            idxs = findall(contains(string(invalidname)), allnames)
            if length(idxs) == 1
                allsyms[only(idxs)]
            else
                @warn "Could not resolve invalid symbol $invalidname, options are $(allsyms[idxs])"
            end
        end
        metadatasubs[invalids] = valid
    end

    fixedeqs = [Symbolics.fast_substitute(eq, metadatasubs) for eq in invalid_eqs]
    if !isempty(check_metadata(fixedeqs))
        @warn "Some transformation droped metadata ($missingmetadata)! Could NOT be fixed. $(check_metadata(fixedeqs))"
    else
        @warn "Some transformation droped metadata ($missingmetadata)! Could be fixed."
    end
    invalid_eqs .= fixedeqs
end

export structural_simplify_bus
function structural_simplify_bus(sys, busbar=:busbar)
    io = ([getproperty(sys, busbar; namespace=false).i_r, getproperty(sys, busbar; namespace=false).i_i], [])
    structural_simplify(sys, io)[1]
end
