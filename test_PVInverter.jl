using Pkg
Pkg.instantiate()
using PowerDynamics

node_list=[]
    #append!(node_list,[SlackAlgebraic(U=1.)])
    append!(node_list, [SwingEqLVS(H=0.5, P=-1, D=0.1, Ω=50,Γ=0.1,V=1.)])
    append!(node_list,[PVInverterWithFrequencyControl(I_n=1.,k_PLL=1,f=50,f_s=50.02,T_m=0.1,k_P=0.2)])
    line_list=[]
    append!(line_list,[StaticLine(from=1,to=2,Y=-1im/0.02)])
    #append!(line_list,[StaticLine(from=2,to=3,Y=-1im/0.02)])


powergrid = PowerGrid(node_list,line_list)

#vertices = map(construct_vertex, powergrid.nodes)
#println(v.mass_matrix for v in vertices)

operationpoint = find_operationpoint(powergrid)

#using DiffEqCallbacks: SavedValues, SavingCallback
using OrdinaryDiffEq: get_du, ODEProblem, Rodas4, solve
p = Perturbation(1, :ω, Inc(0.5))
    x1 = p(operationpoint)
    problem = ODEProblem{true}(rhs(powergrid),x1.vec,(0.0,3.))
    result2 = solve(problem,Rodas4(autodiff=false),saveat=1e-4)#,callback=cb)

#saved_values = SavedValues(Float64, Vector{Float64})
#cb = SavingCallback((u,t,integrator)->get_du(integrator), saved_values, saveat = 1e-6:1e-6:1.)
#solution = solve(problem, Rodas4(autodiff=true),saveat=1e-6, callback=cb)#reltol=1e-8,abstol=1e-8,
#PowerGridSolution(solution, pg)
result = simulate(Perturbation(1, :ω, Inc(0.5)), powergrid, operationpoint, timespan = (0.0,3.))

#import DiffEqBase: solve
#Susing DiffEqCallbacks, OrdinaryDiffEq, LinearAlgebra
#saved_values = SavedValues(Float64, Vector{Float64})
#cb = SavingCallback((u,t,integrator)->(get_du(integrator), saved_values))


include("plotting.jl")
plot_res(result,powergrid,2)

p1 = plot(result, [2], :θ_PLL)
    p2 = plot(result, [1], :ω)
    p3 = plot(result2(range(0.,stop=3.,length=10000),Val{1},idxs=[6]))
    #p4 = plot(test[:,2])
    p5 = plot(result2(range(0.,stop=3.,length=10000),Val{1},idxs=[3]))
    #p6 = plot(test[:,3])
    plot(
    p1,p2,p3,p5;
    layout=(4,1),
    size = (500, 500))



#using DiffEqCallbacks, OrdinaryDiffEq, LinearAlgebra, ODEInterfaceDiffEq
#prob = ODEProblem{true}(rhs(powergrid), x1.vec, (0.0,1.0))
#saved_values = SavedValues(Float64, Array{Float64})
#cb = SavingCallback((u,t,integrator)-> integrator(t,Val{1}), saved_values,saveat=1e-6:1e-6:1.)
#sol = solve(prob, Ros4LStab(autodiff=false),saveat=1e-7,callback=cb)
#test=hcat(saved_values.saveval...)'
#plot(test[:,6])
