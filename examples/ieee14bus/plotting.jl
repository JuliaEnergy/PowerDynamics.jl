using PowerDynamics: symbolsof
using Plots: plot, get_color_palette, plot_color
using LaTeXStrings: latexstring, @L_str

function create_plot(sol)
    swing_indices = findall(n -> :ω ∈ symbolsof(n), sol.powergrid.nodes);
    ω_labels = reshape([latexstring(string(raw"\omega", "_{$i}")) for i=swing_indices], (1, length(swing_indices)))
    p_labels = reshape([latexstring(string(raw"p", "_{$i}")) for i=swing_indices], (1, length(swing_indices)))
    q_labels = reshape([latexstring(string(raw"q", "_{$i}")) for i=swing_indices], (1, length(swing_indices)))
    v_labels = reshape([latexstring(string(raw"v", "_{$i}")) for i=swing_indices], (1, length(swing_indices)))

    pl_v = plot(sol, swing_indices, :v, legend = (0.8, 0.7), ylabel=L"V [p.u.]",label = v_labels)
    pl_p = plot(sol, swing_indices, :p, legend = (0.8, 0.7), ylabel=L"p [p.u.]", label=p_labels)
    pl_q = plot(sol, swing_indices, :q, legend = (0.8, 0.7), ylabel=L"q [p.u.]", label=q_labels)
    pl_ω = plot(sol, swing_indices, :ω, legend = (0.8, 0.7), ylabel=L"\omega \left[rad/s\right]", label=ω_labels)
    pl = plot(pl_ω, pl_v, pl_p, pl_q;
        layout=(2,2),
        size = (1000, 500),
        lw=3,
        xlabel=L"t[s]"
        )
    return pl
end
