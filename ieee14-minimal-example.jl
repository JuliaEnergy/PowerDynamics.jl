using CSV
using DataFrames
using PowerDynBase
using LaTeXStrings
using Plots
using SparseArrays
using NetworkDynamics
using DifferentialEquations

powergrid = read_network_from_csv("IEEE14_busses.csv", "IEEE14_lines.csv")

# find the fixed point = normal operation point
system_size = systemsize(powergrid)
ic = find_operationpoint(powergrid.network_dynamics, ones(system_size))

begin
    # just ensure the correct admittance laplacian is used
    # in case the code was not executed in order
    #rhs = NetworkRHS(g)
    #rhs.LY[:] = LY_default
    # define the initial condition as a perturbation from the fixed point
    x0 = copy(ic)
    x0[3] += 0.2 # perturbation on the ω of the first node
    #x0[n, :int, i] : access to the i-th internal variables of the n-th node
    #x0[n, :u] : access to the complex voltage of the n-th node
    #x0[n, :v] : access to the magnitude of the voltage of the n-th node
    #x0[n, :φ] : access to the voltage angle of the n-th node

    timespan = (0.0,0.3)
    # solve it
    #solution = solve(powergrid, x0, timespan)
    problem = ODEProblem{true}(powergrid.network_dynamics,x0,timespan)
    solution = solve(problem, Rodas4(autodiff=false), force_dtmin=true)
    pl = plot(solution, vars = [s for s in powergrid.network_dynamics.syms if occursin("ω", string(s))])
    display(pl)
end

#swing_indices = findall(n -> isa(n, SwingEqLVS), powergrid.nodes)
#swing_indices = findall(busses_df[:type] .== "G")
#ω_colors = reshape(Plots.get_color_palette(:auto, plot_color(:white), 8)[swing_indices], (1,length(swing_indices)))
#ω_labels = reshape([latexstring(string(raw"\omega", "_{$i}")) for i=swing_indices], (1, length(swing_indices)))
#p_labels = reshape([latexstring(string(raw"p", "_{$i}","[$(busses_df[i,:type])]")) for i=1:length(node_list)], (1, length(node_list)))

################################################################
# plotting the network representing the power grid
# check the two below for plotting graphs
# using LightGraphs
# using GraphPlot
# g = Graph(Array(LY).!=0)
# gplot(g)
################################################################


begin
    #pl_v = plot(solution, :, :v, legend = (0.4, 1.), ylabel=L"V [p.u.]")
    #pl_p = plot(sol, :, :p, legend = (0.8, 0.95), ylabel=L"p [p.u.]", label=p_labels)
    #pl_ω = plot(solution, swing_indices, :ω, legend = (0.8, 0.7), ylabel=L"\omega \left[rad/s\right]", label=ω_labels, color=ω_colors)
    #pl = plot(
         #pl_v,pl_ω;
        #layout=(2,1),
        #size = (500, 500),
        #lw=3,
        #xlabel=L"t[s]"
    #)
    #savefig(pl, "ieee14-frequency-perturbation.pdf")
    #savefig(pl, "ieee14-frequency-perturbation.svg")
    #display(pl)
end


begin
    #fault_line = (1, 5)
    #@assert fault_line[1] < fault_line[2] "order important to kill the line in the data frame"
    #idxs = findall(.~((lines_df[:from] .== fault_line[1]) .& (lines_df[:to] .== fault_line[2])))
    #lines_df_fault = copy(lines_df)[idxs, :]
    #LY_fault = linedf2LY(lines_df_fault, length(node_list))
    #@show LY_default - LY_fault
    #rhs = NetworkRHS(g)
    #rhs.LY[:] = LY_fault
    # start from the fixed point of the original system
    #x0 = copy(fp)
    #timespan = (0.0,1.0)
    # # solve it
    #sol = solve(g, x0, timespan);
end



begin
    # pl_v = plot(sol, :, :v, legend = (0.4, 1.), ylabel=L"V [p.u.]")
    #pl_p = plot(sol, :, :p, legend = (0.8, 0.95), ylabel=L"p [p.u.]", label=p_labels)
    #pl_ω = plot(sol, swing_indices, :ω, legend = (0.8, 0.7), ylabel=L"\omega \left[rad/s\right]", label=ω_labels, color=ω_colors)
    #pl = plot(
        # pl_v,
        #pl_p, pl_ω;
        #layout=(2,1),
        #size = (500, 500),
        #lw=3,
        #xlabel=L"t[s]"
        #)
    #savefig(pl, "ieee14-line-tripping.pdf")
    #savefig(pl, "ieee14-line-tripping.svg")
    #display(pl)
end
