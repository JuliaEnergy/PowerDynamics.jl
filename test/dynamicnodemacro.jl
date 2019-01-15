
include("testing_base.jl")

begin
    _removelinereferences(b) = b
    function _removelinereferences(b::Expr)
        if ~(b.head in (:block, :macrocall))
            return b
        end
        b.args = [ el for el in b.args if ~(el isa LineNumberNode) ]
        b
    end
    removelinereferences(q::Expr) = MacroTools.postwalk( _removelinereferences , q)
    rlr = removelinereferences
    _removesingleblocks(b) = b
    function _removesingleblocks(b::Expr)
        if ~(b.head in (:block, )) || length(b.args) != 1
            return b
        end
        b.args[1]
    end
    removesingleblocks(q::Expr) = MacroTools.postwalk( _removesingleblocks , q)
    rsb = removesingleblocks
end

struct_def = PowerDynBase.buildparameterstruct(:MyParams, [:p1, :p2])
struct_def_tests = :(struct MyParams <: AbstractNodeParameters
        p1
        p2
        MyParams(; p1, p2) = new(p1, p2)
        end)
@test struct_def |> rlr == struct_def_tests |> rlr |> rsb
internals = PowerDynBase.getinternalvars(Val{:OrdinaryNodeDynamics}, :([[a, da], [b, db]]))
@test internals == (vars=[:a, :b], dvars=[:da, :db])
@test PowerDynBase.getinternalvars(Val{:OrdinaryNodeDynamicsWithMass}, :([[a, da], [b, db]])) == internals
@test_throws NodeDynamicsError PowerDynBase.getinternalvars(Val{:unknown}, :([[a, da], [b, db]]))
dynamicscall = :(OrdinaryNodeDynamicsWithMass(a=1, b=3))
funcbody = quote
    funcbodyline1
    funcbodyline2
    end
prep = quote
        prepline1
        prepline2
    end
cndfunction = quote
        function construct_node_dynamics(par::MyParams)
                p1 = par.p1
                p2 = par.p2
                prepline1
                prepline2
        end
end |> rlr |> rsb
PowerDynBase.cndfunction_builder!(Val{:OrdinaryNodeDynamicsWithMass},
        internals,
        dynamicscall,
        funcbody,
        cndfunction)
cnd_function_test = quote
        function construct_node_dynamics(par::MyParams)
                p1 = par.p1
                p2 = par.p2
                prepline1
                prepline2
                function rhs!(dint::AbstractVector, u, i, int::AbstractVector, t)
                        a = int[1]
                        b = int[2]
                        funcbodyline1
                        funcbodyline2
                        try
                                dint[1] = da
                                dint[2] = db
                                return du
                        catch e
                                if typeof(e) === UndefVarError
                                        throw(NodeDynamicsError("you need to provide $(e.var)"))
                                else
                                        throw(e)
                                end
                        end
                end
                OrdinaryNodeDynamicsWithMass(a=1, b=3, rhs=rhs!, n_int=2, parameters=par)
        end
end |> rlr |> rsb
@test cndfunction |> rlr |> rsb == cnd_function_test
symbolsof_fct = PowerDynBase.generate_symbolsof_fct(Val{:OrdinaryNodeDynamicsWithMass}, :MyParams, internals) |> rlr |> rsb
symbolsof_fct_test = :(symbolsof(::Type{MyParams}) = ODENodeSymbols([:a, :b], [:da, :db])) |> rlr |> rsb
@test symbolsof_fct == symbolsof_fct_test
full_macro_return = PowerDynBase.DynamicNode(:(MyParams(p1, p2) <: OrdinaryNodeDynamicsWithMass(a=1, b=3)), prep, :([[a, da], [b, db]]), funcbody ) |> rlr |> rsb
full_macro_return_test = quote
        @__doc__ $(struct_def_tests)
        $(cnd_function_test)
        $(symbolsof_fct_test)
end |> rlr |> rsb
@test full_macro_return == full_macro_return_test
