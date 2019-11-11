using PowerDynamics
using LaTeXStrings
using Plots
include("helpers.jl")

function plot_res(result, powergrid,disturbed_node,simulation_time)
    ω_indices = findall(n -> isa(n, SwingEqLVS), powergrid.nodes)
    append!(ω_indices,findall(n -> isa(n, VSIMinimal), powergrid.nodes))
    append!(ω_indices,findall(n -> isa(n, FourthOrderEq), powergrid.nodes))
    append!(ω_indices,findall(n -> isa(n, WindTurbineGenType4_RotorControl), powergrid.nodes))
    #append!(ω_indices,findall(n -> isa(n, CurtailedPowerPlantWithInertia), powergrid.nodes))


    ω_PLL_indices = findall(n -> isa(n, WindTurbineGenType4_RotorControl), powergrid.nodes)
    #append!(ω_PLL_indices,findall(n -> isa(n, CurtailedPowerPlantWithInertia), powergrid.nodes))

    ω_colors = reshape(Plots.get_color_palette(:auto, plot_color(:white), 8)[ω_indices], (1,length(ω_indices)))
    ω_labels = reshape([latexstring(string(raw"\omega", "_{$i}")) for i=ω_indices], (1, length(ω_indices)))
    ω_PLL_labels = reshape([latexstring(string(raw"\omega", "_{$i}")) for i=ω_PLL_indices], (1, length(ω_PLL_indices)))
    p_labels = reshape([latexstring(string(raw"p", "_{$i}")) for i=1:length(powergrid.nodes)], (1, length(powergrid.nodes)))
    pl_p = plot(result, :, :p, legend = (0.8, 0.95), ylabel=L"p [p.u.]", label=p_labels)
    pl_ω = plot(result, ω_indices, :ω, legend = (0.8, 0.7), ylabel=L"\omega \left[rad/s\right]", label=ω_labels, color=ω_colors)
    ω_PLL=[]
    for i in ω_PLL_indices
        append!(ω_PLL,[result.dqsol(result.dqsol.t,Val{1},idxs=variable_index(powergrid.nodes,i , :θ_PLL)).u])
    end
    pl_ω_PLL = plot(result.dqsol.t,ω_PLL, legend = (0.8, 0.7), ylabel=L"\omega_{PLL} \left[rad/s\right]", label=ω_PLL_labels)
    pl = plot(
        pl_p, pl_ω,pl_ω_PLL;
        layout=(3,1),
        size = (500, 500),
        lw=3,
        xlabel=L"t[s]"
    )
    display(pl)
end

function plot_res_compare(result1,result2,powergrid1,powergrid2,disturbed_node)
    ω_indices1 = findall(n -> isa(n, SwingEqLVS), powergrid1.nodes)
    append!(ω_indices1,findall(n -> isa(n, VSIMinimal), powergrid1.nodes))
    append!(ω_indices1,findall(n -> isa(n, FourthOrderEq), powergrid1.nodes))
    ω_colors1 = reshape(Plots.get_color_palette(:auto, plot_color(:white), 8)[ω_indices1], (1,length(ω_indices1)))
    ω_labels1 = reshape([latexstring(string(raw"\omega", "_{$i}")) for i=ω_indices1], (1, length(ω_indices1)))
    p_labels1 = reshape([latexstring(string(raw"p", "_{$i}")) for i=1:length(powergrid1.nodes)], (1, length(powergrid1.nodes)))
    pl_p1 = plot(result1, :, :p, legend = (0.8, 0.95), ylabel=L"p [p.u.]", label=p_labels1)
    pl_ω1 = plot(result1, ω_indices1, :ω, legend = (0.8, 0.7), ylabel=L"\omega \left[rad/s\right]", label=ω_labels1, color=ω_colors1)

    ω_indices2 = findall(n -> isa(n, SwingEqLVS), powergrid2.nodes)
    append!(ω_indices2,findall(n -> isa(n, VSIMinimal), powergrid2.nodes))
    append!(ω_indices2,findall(n -> isa(n, FourthOrderEq), powergrid2.nodes))
    ω_colors2 = reshape(Plots.get_color_palette(:auto, plot_color(:white), 8)[ω_indices2], (1,length(ω_indices2)))
    ω_labels2 = reshape([latexstring(string(raw"\omega", "_{$i}")) for i=ω_indices2], (1, length(ω_indices1)))
    p_labels2 = reshape([latexstring(string(raw"p", "_{$i}")) for i=1:length(powergrid2.nodes)], (1, length(powergrid2.nodes)))
    pl_p2 = plot(result2, :, :p, legend = (0.8, 0.95), ylabel=L"p [p.u.]", label=p_labels2)
    pl_ω2 = plot(result2, ω_indices2, :ω, legend = (0.8, 0.7), ylabel=L"\omega \left[rad/s\right]", label=ω_labels2, color=ω_colors2)

    pl = plot(
        pl_p1, pl_ω1,pl_p2,pl_ω2;
        layout=(2,2),
        size = (1000, 500),
        lw=3,
        xlabel=L"t[s]"
    )
    display(pl)
end
