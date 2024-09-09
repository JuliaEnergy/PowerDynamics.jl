function isterminal(sys::ODESystem)
    vars = Set(getname.(unknowns(sys))) == Set([:u_r, :u_i, :i_r, :i_i])
    inputs = isempty(ModelingToolkit.unbound_inputs(sys))
    outputs = isempty(ModelingToolkit.unbound_outputs(sys))
    vars && inputs && outputs
end
function iscomponentmodel(sys::ODESystem)
    # HACK: use try catch instead of hasproperty https://github.com/SciML/ModelingToolkit.jl/issues/3016
    try
        return isterminal(sys.terminal)
    catch e
        if e isa ArgumentError && contains(e.msg, "variable terminal does not exist")
            return false
        else
            rethrow(e)
        end
    end
end


function isbusbar(sys::ODESystem)
    vars = getname.(unknowns(sys)) ⊇ Set([:u_r, :u_i, :i_r, :i_i, :P, :Q])
    inputs = Set(getname.(ModelingToolkit.unbound_inputs(sys))) == Set([:i_r, :i_i])
    outputs = Set(getname.(ModelingToolkit.unbound_outputs(sys))) == Set([:u_r, :u_i])
    vars && inputs && outputs
end
function isbusmodel(sys::ODESystem)
    # HACK: use try catch instead of hasproperty https://github.com/SciML/ModelingToolkit.jl/issues/3016
    try
        return isbusbar(sys.busbar)
    catch e
        if e isa ArgumentError && contains(e.msg, "variable busbar does not exist")
            return false
        else
            rethrow(e)
        end
    end
end


function islineend(sys::ODESystem)
    vars = getname.(unknowns(sys)) ⊇ Set([:u_r, :u_i, :i_r, :i_i])
    inputs = Set(getname.(ModelingToolkit.unbound_inputs(sys))) == Set([:u_r, :u_i])
    outputs = Set(getname.(ModelingToolkit.unbound_outputs(sys))) == Set([:i_r, :i_i])
    vars && inputs && outputs
end
function islinemodel(sys::ODESystem)
    # HACK: use try catch instead of hasproperty https://github.com/SciML/ModelingToolkit.jl/issues/3016
    try
        return islineend(sys.src)
    catch e
        if e isa ArgumentError && contains(e.msg, "variable src does not exist")
            return false
        else
            rethrow(e)
        end
    end
    try
        return islineend(sys.dst)
    catch e
        if e isa ArgumentError && contains(e.msg, "variable dst does not exist")
            return false
        else
            rethrow(e)
        end
    end
end

function isbranchmodel(sys::ODESystem)
    # HACK: use try catch instead of hasproperty https://github.com/SciML/ModelingToolkit.jl/issues/3016
    try
        return isterminal(sys.src)
    catch e
        if e isa ArgumentError && contains(e.msg, "variable src does not exist")
            return false
        else
            rethrow(e)
        end
    end
    try
        return isterminal(sys.dst)
    catch e
        if e isa ArgumentError && contains(e.msg, "variable dst does not exist")
            return false
        else
            rethrow(e)
        end
    end
end
