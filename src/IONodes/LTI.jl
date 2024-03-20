export identify_lti, balresid

using ControlSystems
using BlockSystems: @check
using ModelingToolkit: value
using ModelingToolkit.SymbolicUtils: Symbolic
using LinearAlgebra

function _state_matrix(expr, vars)
    if eltype(expr) <: Equation
        expr = getproperty.(expr, :rhs)
    end

    A = Matrix(undef, length(expr), length(vars))
    fill!(A, NaN)
    for (ieq, ex) in enumerate(expr)
        # ex = simplify(ex, expand=true)
        for (istate, state) in enumerate(vars)
            coeff = (ex - substitute(ex, state => 0))/state |> simplify
            if !isempty(Set(get_variables(coeff)) ∩ Set(vars))
                error("Could not properly extract $state from $(ex), got $coeff")
            end
            A[ieq, istate] = coeff
        end
    end
    return _narrow_type(A)
end

"""
    A,B,C,D,x,y,u=identify_lti(blk::IOBlock)

Identify the matrices A,B,C and D from IOBlock.
"""
function identify_lti(blk::IOBlock)
    inputs  = blk.inputs
    outputs = blk.outputs
    output_rhs = similar(outputs, Any)
    odestates = Symbolic[]
    odeidx    = Int[]
    algstates = Symbolic[]
    algidx    = Int[]
    for (i, eq) in enumerate(equations(blk))
        (type, var) = BlockSystems.eq_type(eq)
        if type == :explicit_diffeq
            push!(odestates, var)
            push!(odeidx, i)
        elseif type == :explicit_algebraic
            @assert var ∈ Set(outputs) "Explicit algebraic equations musst represent outputs!"
            push!(algstates, var)
            push!(algidx, i)
        else
            error("Only pure LTI systems. Can not handle $type-type equations!")
        end
        #check if variable is an output
        outidx = findfirst(isequal(var), outputs)
        if !isnothing(outidx) # is output variable
            if type == :explicit_diffeq
                output_rhs[outidx] = var
            elseif type == :explicit_algebraic
                output_rhs[outidx] = eq.rhs
            else
                error()
            end
        end
    end

    for eq in equations(blk)
        @assert isempty(Set(get_variables(eq.rhs)) ∩ algstates) "The rhs of the equations must not contain algebraic states!"
    end

    A = _state_matrix(equations(blk)[odeidx], odestates)
    B = _state_matrix(equations(blk)[odeidx], inputs)
    C = _state_matrix(output_rhs, odestates)
    D = _state_matrix(output_rhs, inputs)

    Avars = isempty(A) ? Set() : Set(mapreduce(get_variables, vcat, A))
    Bvars = isempty(B) ? Set() : Set(mapreduce(get_variables, vcat, B))
    Cvars = isempty(C) ? Set() : Set(mapreduce(get_variables, vcat, C))
    Dvars = isempty(D) ? Set() : Set(mapreduce(get_variables, vcat, D))

    @assert Avars ⊆ Set(blk.iparams) "Matrix A contains non-parameters. Thats an error! $Avars"
    @assert Bvars ⊆ Set(blk.iparams) "Matrix B contains non-parameters. Thats an error! $Bvars"
    @assert Cvars ⊆ Set(blk.iparams) "Matrix C contains non-parameters. Thats an error! $Cvars"
    @assert Dvars ⊆ Set(blk.iparams) "Matrix D contains non-parameters. Thats an error! $Dvars"

    if eltype(A) == Any
        # @warn "There appear to be `Any` terms in A"
        A = fixdiv(A)
    end
    if eltype(B) == Any
        # @warn "There appear to be `Any` terms in B"
        B = fixdiv(B)
    end
    if eltype(C) == Any
        # @warn "There appear to be `Any` terms in C"
        C = fixdiv(C)
    end
    if eltype(D) == Any
        # @warn "There appear to be `Any` terms in D"
        D = fixdiv(D)
    end

    return A, B, C, D, odestates, outputs, inputs
end

function fixdiv(A)
    for idx in eachindex(A)
        if A[idx] isa SymbolicUtils.Div
            try
                A[idx] = eval(Meta.parse(repr(A[idx])))
            catch e
                if !isa(e, UndefVarError)
                    rethrow(e)
                end
            end
        end
    end
    A = _narrow_type(A)
end

function IOBlock(A::Matrix,B::Matrix,C::Matrix,D::Matrix,x,y,u; name=gensym(:LTI), warn=BlockSystems.WARN[], rem_eqs=Equation[])
    @assert size(A)[1] == size(A)[2] == size(B)[1] == size(C)[2]  == length(x)
    @assert size(B)[2] == size(D)[2] == length(u)
    @assert size(C)[1] == size(D)[1] == length(y)
    iv_candidate = unique(vcat(Symbolics.arguments.(value.(x)),
                               Symbolics.arguments.(value.(y)),
                               Symbolics.arguments.(value.(u))))
    @assert length(iv_candidate) == 1 && length(iv_candidate[begin]) == 1
    iv = iv_candidate[1][1]
    dt = Differential(iv)
    eqs = vcat(dt.(x) .~ A*x + B*u,
               y .~ C*x + D*u)
               @show
    filter!(eq -> !isequal(eq.lhs, eq.rhs), eqs)

    IOBlock(eqs, u, y; name, warn, rem_eqs)
end

"""
    balresid(blk::IOBlock, order; warn=BlockSystems.WARN[], verbose=false, reconstruct=false)

This function creates a model reduced version of `blk` of given order using the balanced residualization method.
Only works for linear systems. Internaly, it uses the `baltrunc` function from `ControlSystems.jl`
to apply the method to an LTI. The LTI system matrices A,B,C and D are determined symbolicially from the
`IOBlock` object.

If `reconstruct=true` the resulting system will include `removed_eqs` to reconstruct the original states from the projected states `z_i`.

"""
function balresid(blk::IOBlock, order; warn=BlockSystems.WARN[], verbose=false, reconstruct=false, getT=false)
    @check isempty(blk.iparams) "In order to balresid a system it can not have any internal parameters!"
    A, B, C, D, x, y, u = identify_lti(blk)
    verbose && @info "Block is LTI of order $(length(x))"

    if length(x) <= order
        throw(ArgumentError("Cannot reduce a system of order $(length(x)) to order $order."))
    end

    ss = StateSpace(A,B,C,D)
    ss_bal, G, _ = ControlSystems.baltrunc(ss; residual=true, n=order)
    _, _, T = ControlSystems.balreal(ss)

    t = get_iv(blk)
    z = Num[]
    for i in 1:length(x)
        zs = subscript(:z, i)
        append!(z, @variables $zs(t))
    end

    z_r = z[begin:order]
    z_t = z[order+1:end]
    A_r = ss_bal.A
    B_r = ss_bal.B
    C_r = ss_bal.C
    D_r = ss_bal.D

    ## check results
    An = T*A*inv(T)
    Bn = T*B
    Cn = C*inv(T)
    Dn = D

    A11 = An[1:order, 1:order]
    A12 = An[1:order, order+1:end]
    A21 = An[order+1:end, 1:order]
    A22 = An[order+1:end, order+1:end]
    @assert [A11 A12; A21 A22] == An

    B1 = Bn[1:order, :]
    B2 = Bn[order+1:end, :]
    @assert [B1; B2] == Bn

    C1 = Cn[:, 1:order]
    C2 = Cn[:, order+1:end]
    @assert [C1 C2] == Cn

    # get matrices for reduced model
    @assert A_r ≈ A11 - A12*inv(A22)*A21
    @assert B_r ≈ B1 - A12*inv(A22)*B2
    @assert C_r ≈ C1 - C2*inv(A22)*A21
    @assert D_r ≈ Dn - C2*inv(A22)*B2
    ##

    verbose && @info "Truncated $(length(x)-order) states!"

    # in the sigular pertubation we assume dot(z_t)=0 but z_t is given by constraint
    # we calculate z_t based on z_r and u
    balA = T*A*inv(T)
    A21 = balA[order+1:end, 1:order]
    A22 = balA[order+1:end, order+1:end]
    B2  = (T*B)[order+1:end, :]
    z_t = -inv(A22)*A21*z_r -inv(A22)*B2*u

    outidx = findall(x -> x ∈ Set(y), x)
    state_estimates = x .~ inv(T)*vcat(z_r, z_t)
    deleteat!(state_estimates, outidx) # outputs are still present in the full system

    if reconstruct
        substitutions = Dict(eq.lhs => eq.rhs for eq in state_estimates)
        old_rem_eqs = map(eq->BlockSystems.eqsubstitute(eq, substitutions), blk.removed_eqs)
        rem_eqs = vcat(state_estimates, old_rem_eqs)
    else
        rem_eqs = Equation[]
    end

    name = Symbol(string(blk.name) * "_trunc_$order")

    blk = IOBlock(A_r, B_r, C_r, D_r, z_r, y, u; rem_eqs, name, warn)
    return getT ? (blk, T) : blk
end

function ControlSystems.balreal(blk::IOBlock)
    A, B, C, D, x, y, u = identify_lti(blk)
    ss = StateSpace(A,B,C,D)
    ss_bal, G, T = ControlSystems.balreal(ss)

    return G, T, x
end

function ControlSystems.StateSpace(iob::IOBlock)
    A,B,C,D = identify_lti(iob)
    StateSpace(A,B,C,D)
end

function _narrow_type(A::AbstractArray)
    isempty(A) && return A
    elt = mapreduce(typeof, promote_type, A)
    convert.(elt, A)
end

function subscript(s, i, add...)
    Symbol(s, _subscript(i), add...)
end

function _subscript(i::Int)
    dig = reverse(digits(i))
    String(map(i -> Char(0x02080 + i), dig))
end
