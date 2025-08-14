pin_parameters(sys::System) = pin_parameters(sys, defaults(sys))
function pin_parameters(sys::System, sub::Pair)
    pin_parameters(sys, Dict(sub.first => sub.second))
end
function pin_parameters(sys::System, subs::Dict)
    _dict = Dict(_resolve_to_namespaced_symbolic(sys, k) => v for (k,v) in subs)

    if any(iscall, values(_dict))
        @info "Parameter map for pin_parameters cannot include any terms in the rhs. Try to resolve..."
        __dict = Dict()
        for (k,v) in _dict
            _fullname = string(getname(k))
            newname = replace(_fullname, r"^"*string(get_name(sys))*"₊" => "")
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

function _pinparameters(sys::System, _subs::Dict)
    isempty(_subs) && return sys
    iscomplete(sys) && ArgumentError("pin_parameters does not handle complete/simplified systems")

    subs = Dict{Symbolic, Any}()
    for (k,v) in _subs
        _fullname = string(getname(k))
        contains(_fullname, r"^"*string(get_name(sys))*"₊") || continue
        newname = replace(_fullname, r"^"*string(get_name(sys))*"₊" => "")
        subs[rename(k, Symbol(newname))] = v
    end

    applicable = Dict{Symbolic, Any}()
    applicable_subs = intersect(keys(subs), parameters(sys))
    for k in applicable_subs
        applicable[k] = subs[k]
    end

    if isempty(applicable)
        _eqs = unwrap(get_eqs(sys))
        _observed = unwrap(get_observed(sys))
        _defaults = get_defaults(sys)
    else
        _eqs = [fixpoint_sub(eq, applicable) for eq in unwrap(get_eqs(sys))]
        # FIXME: check something drops the metadata
        fix_metadata!(_eqs, sys)

        _observed = [fixpoint_sub(eq, applicable) for eq in unwrap(get_observed(sys))]
        fix_metadata!(_observed, sys)

        @argcheck unwrap(isempty(get_ctrls(sys))) "pin_parameters does not handle control variables"
        _defaults = copy(get_defaults(sys))
        for k in keys(applicable)
            delete!(_defaults, k)
        end
        @argcheck isnothing(get_initializesystem(sys)) "pin_parameters does not handle initializesystem"
        @argcheck unwrap(isempty(get_continuous_events(sys))) "pin_parameters does not handle continuous events"
        @argcheck unwrap(isempty(get_discrete_events(sys))) "pin_parameters does not handle discrete events"
        @argcheck unwrap(isempty(get_parameter_dependencies(sys))) "pin_parameters does not parameter_dependencies"
        @argcheck isnothing(get_discrete_subsystems(sys)) "pin_parameters does not handle discrete subsystems"
        @argcheck isnothing(get_solved_unknowns(sys)) "pin_parameters does not handle solved unknowns"
    end

    _systems = [_pinparameters(s, subs) for s in get_systems(sys)]

    newsys = System(_eqs, get_iv(sys);
            # controls = Num[],
            observed = _observed,
            systems = _systems,
            tspan = get_tspan(sys),
            name = get_name(sys),
            defaults = _defaults,
            guesses = get_guesses(sys),
            # initializesystem = nothing,
            # initialization_eqs = Equation[],
            schedule = get_schedule(sys),
            connector_type = get_connector_type(sys),
            preface = get_preface(sys),
            # continuous_events = nothing,
            # discrete_events = nothing,
            # parameter_dependencies = nothing,
            checks = true,
            metadata = ModelingToolkit.get_metadata(sys),
            gui_metadata = get_gui_metadata(sys))
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
