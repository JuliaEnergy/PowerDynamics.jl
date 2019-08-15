Pkg.instantiate()
cd(@__DIR__)
using PowerDynamics

begin
    node_list = []
    append!(node_list, [SlackAlgebraic(
        U=1.
                )])
    append!(node_list, [CSIMinimal(
        I_r=1.,Y_n=-.4
                )])
    append!(node_list, [SwingEqLVS(
                    H=0.1,P=1,D=0.02,Γ=0.5,V=1,Ω=2π*50,Y_n=-0.5*1im+0.5
                            )])
    line_list = []
    append!(line_list,[StaticLine(Y=-1/0.2*1im,from=1,to=2)])
    append!(line_list,[StaticLine(Y=-1/0.2*1im,from=2,to=3)])

    powergrid = PowerGrid(node_list,line_list)
end

operationpoint = find_operationpoint(powergrid)

# TODO: if I assign wrong node for perturbation error message is not understandable, sth with json error
solution = simulate(Perturbation(3, :ω, Inc(0.2)), powergrid, operationpoint, timespan = (0.0,1.))

using Plots
plot(solution,:,:p)

plot(solution,3,:ω)
