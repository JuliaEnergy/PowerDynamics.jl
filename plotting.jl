using PowerDynBase
using LaTeXStrings
using Plots

function plot_res(result, powergrid)
    swing_indices = findall(n -> isa(n, SwingEqLVS), powergrid.nodes)
    ω_colors = reshape(Plots.get_color_palette(:auto, plot_color(:white), 8)[swing_indices], (1,length(swing_indices)))
    ω_labels = reshape([latexstring(string(raw"\omega", "_{$i}")) for i=swing_indices], (1, length(swing_indices)))
    p_labels = reshape([latexstring(string(raw"p", "_{$i}")) for i=1:length(powergrid.nodes)], (1, length(powergrid.nodes)))

    pl_p = plot(result, :, :p, legend = (0.8, 0.95), ylabel=L"p [p.u.]", label=p_labels)
    pl_ω = plot(result, swing_indices, :ω, legend = (0.8, 0.7), ylabel=L"\omega \left[rad/s\right]", label=ω_labels, color=ω_colors)
    pl = plot(
        pl_p, pl_ω;
        layout=(2,1),
        size = (500, 500),
        lw=3,
        xlabel=L"t[s]"
    )
    display(pl)
end
