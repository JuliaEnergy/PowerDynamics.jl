using PowerDynamics
using LaTeXStrings
using Plots

function plot_res(result, powergrid,disturbed_node)
    ω_indices = findall(n -> isa(n, SwingEqLVS), powergrid.nodes)
    append!(ω_indices,findall(n -> isa(n, VSIMinimal), powergrid.nodes))
    #append!(ω_indices,findall(n -> isa(n, PVInverterWithFrequencyControl), powergrid.nodes))
    print(ω_indices)
    ω_colors = reshape(Plots.get_color_palette(:auto, plot_color(:white), 8)[ω_indices], (1,length(ω_indices)))
    ω_labels = reshape([latexstring(string(raw"\omega", "_{$i}")) for i=ω_indices], (1, length(ω_indices)))
    p_labels = reshape([latexstring(string(raw"p", "_{$i}")) for i=1:length(powergrid.nodes)], (1, length(powergrid.nodes)))
    pl_p = plot(result, [2], :p, legend = (0.8, 0.95), ylabel=L"p [p.u.]", label=p_labels)
    pl_ω = plot(result, ω_indices, :ω, legend = (0.8, 0.7), ylabel=L"\omega \left[rad/s\right]", label=ω_labels, color=ω_colors)
    pl_ϕ =
    pl = plot(
        pl_p,pl_ω;
        layout=(2,1),
        size = (500, 500),
        lw=3,
        xlabel=L"t[s]"
    )
    savefig(pl,"StudyCase_dist_ω"* string(disturbed_node) *"3.pdf")
    display(pl)
end
