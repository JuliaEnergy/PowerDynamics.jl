using NonlinearSolve
using NetworkDynamics: LazyBufferCache
using LinearAlgebra: norm

export UnsuccessfullError
struct UnsuccessfullError <: Exception
   msg::String
end
function Base.showerror(io::IO, err::UnsuccessfullError)
    print(io, "UnsuccesffulError: ")
    print(io, err.msg)
end


export find_operationpoint

function find_operationpoint(nw, guess)
    # nw = safehouse.pg
    # guess = safehouse.ic_guess
    init = deepcopy(guess)
    if any(s -> isnan(s) || ismissing(s), uflat(guess))

        println("Network is not fully initialized from powerflow, trying to initialize it component wise...")
        @assert all(NetworkDynamics.isstatic.(nw.im.edgef))
        _, aggbuf = NetworkDynamics.get_ustacked_buf(guess)
        @assert !any(isnan.(aggbuf))
        for i in 1:nv(nw)
            vertexf = nw.im.vertexf[i]
            _u = guess.v[i, :]
            _p = guess.p.v[i, :]
            _e = aggbuf[nw.im.v_aggr[i]]
            uinit, success = initialize_vertex(vertexf, _u, _e, _p, guess.t)
            init.v[i, :] .= uinit
        end
    end

    du = [Inf for i in 1:dim(nw)]
    nw(du, uflat(init), pflat(init), init.t)
    residual = maximum(abs.(du))
    if residual < 1e-8
        println("Component wise initialization reached NW steady state. Residual: $residual")
        return init
    else
        println("Component wise initialisation did not reach state state Residual: $residual. Attempt dynamic relaxation...")
        relax = find_fixpoint(nw, init)
        du = [Inf for i in 1:dim(nw)]
        nw(du, uflat(relax), pflat(relax), relax.t)
        relaxresidual = maximum(abs.(du))
        if relaxresidual > 1e-8
            throw(UnsuccessfullError("Dynamic relaxation did not converge. Residual: $relaxresidual"))
        end
        println("Dynamic relaxation reached NW steady state. Residual: $residual")
        return relax
    end
end


function initialize_vertex(v::ODEVertex, u, esum, p, t)
    initmask = map(s -> isnan(s) || ismissing(s), u)

    _f = let u=u, esum=esum, cachepool = LazyBufferCache()
        function(_du, _u, _p)
            _tmp = cachepool[_du, length(u)]
            _tmp .= u
            _tmp[initmask] .= _u
            v.f(_du, _tmp, esum, _p, t)
        end
    end
    u0 = zeros(sum(initmask))
    prob = NonlinearLeastSquaresProblem(
        NonlinearFunction(_f, resid_prototype = similar(u)), u0, p)
    sol = solve(prob)
    success = SciMLBase.successful_retcode(sol)
    if !success
        @warn "Unsuccesfull retcode $(sol.retcode)"
    end
    unew = copy(u)
    unew[initmask] .= sol.u

    du = [Inf for i in 1:length(u)]
    v.f(du, unew, esum, p, t)
    err = maximum(abs.(du))
    if err > 1e-6
        @warn "Initialisation of node $(v.name) lead to residual of $err"
    end

    unew, success
end



# function _initialize_vertex(v::ODEVertex, p; u0=nothing, vmag)
#     g = path_graph(2)
#     slack = ODEVertex(Slack())
#     line = StaticEdge(RMSRLLine(;ω0=2π*50, R=0, L=1e-4));
#     # line = ODEEdge(EMTRLLine(;ω0=2π*50, R=0, L=1e-4));
#     nd = network_dynamics([slack, v], line, g)
#     if u0==nothing
#         # need to scale down initial guess to match the voltage magnitude
#         u0 = vmag .* u0guess(nd; noise=false)
#         nd.syms .=> u0
#         u0[1] = vmag
#         u0[2] = 0.0
#     else
#         u0 = copy(u0)
#         prepend!(u0, [vmag, 0.0])
#         @assert length(u0) == length(nd.syms)
#     end
#     pnd = (p, nothing)

#     prob = SteadyStateProblem(nd, u0, pnd)
#     sol = solve(prob, DynamicSS(Rodas5P()))

#     # prob = ODEProblem(nd, u0, (0,10.01),pnd)
#     # sol = solve(prob, Rodas5P())
#     # Plots.plot(sol)
#     # nd.syms .=> round.(sol[end]; digits=1)
#     # nd.syms .=> round.(sol[1]; digits=1)
#     # return sokl

#     du = fill(NaN, length(sol.u))
#     nd(du, sol.u, pnd, NaN)
#     if maximum(abs.(du)) > 1e-6
#         @warn "Vertex initialization did not converge! Max res = $(maximum(abs.(du)))"
#         @info "Residuals" nd.syms .=> du
#     end
#     return sol.u[3:2+length(v.sym)]
# end

# function _initialize_edge(e::ODEEdge, src, dst, p)
#     u0 = zeros(length(e.sym))
#     prob = NonlinearProblem((du,u,p)->e.f(du,u,src,dst,p,NaN), u0, p)
#     sol = solve(prob, NewtonRaphson())

#     if !SciMLBase.successful_retcode(sol)
#         @warn "Steady state search for edge did not converge!" sol.retcode
#     end
#     sol.u
# end
