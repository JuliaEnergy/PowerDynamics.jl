using PowerDynamics: symbolsof
using Plots: Measures,plot, get_color_palette, plot_color, savefig
using LaTeXStrings: latexstring, @L_str

function plot_res(result, powergrid,disturbed_node,scenario)
    ω_indices = findall(n -> isa(n, SwingEqLVS), powergrid.nodes)
    append!(ω_indices,findall(n -> isa(n, VSIMinimal), powergrid.nodes))
    #append!(ω_indices,findall(n -> isa(n, VSIMinimal_experimental), powergrid.nodes))
    append!(ω_indices,findall(n -> isa(n, PVInverterWithFrequencyControl), powergrid.nodes))
    append!(ω_indices,findall(n -> isa(n, WindTurbineGenType4_RotorControl), powergrid.nodes))
    append!(ω_indices,findall(n -> isa(n, FourthOrderEq), powergrid.nodes))
    #append!(ω_indices,findall(n -> isa(n, GridFormingTecnalia), powergrid.nodes))
    #append!(ω_indices,findall(n -> isa(n, GridFormingTecnalia_experimental), powergrid.nodes))

    #ω_colors = reshape(Plots.get_color_palette(:auto, plot_color(:white), 8)[ω_indices], (1,length(ω_indices)))
    ω_labels = reshape([latexstring(string(raw"\omega", "_{$i}")) for i=ω_indices], (1, length(ω_indices)))
    p_labels = reshape([latexstring(string(raw"p", "_{$i}")) for i=1:length(powergrid.nodes)], (1, length(powergrid.nodes)))
    pl_p = plot(result, :, :p, ylims=(-30.,30.), legend = (0.8, 0.95), ylabel=L"p [p.u.]", label=p_labels)
    pl_q = plot(result, :, :q, ylims=(-20.,20.), legend = (0.8, 0.95), ylabel=L"q [p.u.]", label=p_labels)
    pl_v = plot(result, :, :v, ylims=(0.95,1.05), legend = (0.8, 0.95), ylabel=L"v [p.u.]", label=p_labels)
    pl_ω = plot(result, ω_indices, :ω, legend = (0.8, 0.7), ylims = (-1.0,1.0), ylabel=L"\omega \left[rad/s\right]", label=ω_labels)#, color=ω_colors, left_margin = 4mm)
    pl = plot(
        pl_ω, pl_v, pl_p, pl_q;
        layout=(2,2),
        size = (1000, 700),
        lw=3,
        xlabel=L"t[s]"
    )
    savefig(pl,"Testcase_"*scenario*".pdf")
    display(pl)
end



function create_plot(sol,scenario)
    swing_indices = findall(n -> :ω ∈ symbolsof(n), sol.powergrid.nodes)
    ω_labels = reshape([latexstring(string(raw"\omega", "_{$i}")) for i=swing_indices], (1, length(swing_indices)))
    p_labels = reshape([latexstring(string(raw"p", "_{$i}")) for i=collect(keys(powergrid.nodes))], (1, length(collect(powergrid.nodes))))
    q_labels = reshape([latexstring(string(raw"q", "_{$i}")) for i=collect(keys(powergrid.nodes))], (1, length(collect(powergrid.nodes))))
    v_labels = reshape([latexstring(string(raw"v", "_{$i}")) for i=collect(keys(powergrid.nodes))], (1, length(collect(powergrid.nodes))))
    
    pl_p = plot(sol, :, :p, ylims=(-30.,30.), legend = (0.8, 0.95), ylabel=L"p [p.u.]", label=p_labels)
    pl_q = plot(sol, :, :q, ylims=(-20.,20.), legend = (0.8, 0.95), ylabel=L"q [p.u.]", label=p_labels)
    pl_v = plot(sol, :, :v, ylims=(0.95,1.05), legend = (0.8, 0.95), ylabel=L"v [p.u.]", label=p_labels)
    pl_ω = plot(sol, swing_indices, :ω, legend = (0.8, 0.7), ylabel=L"\omega \left[rad/s\right]", label=ω_labels)
    pl = plot(pl_ω, pl_v, pl_p, pl_q;
        layout=(2,2),
        size = (1000, 500),
        lw=3,
        xlabel=L"t[s]"
        )
    savefig(pl,"Testcase_"*scenario*".pdf")
    return pl
end