using Plots
using MarinePowerDynamics: FourthOrderEq, LinearPTO

function create_plot(sol)
    generator_indices = findall(bus -> typeof(bus) == FourthOrderEq || typeof(bus) == LinearPTO, powergrid.nodes)
    labels = reshape(generator_indices,(1,length(generator_indices)))

    pl_v = plot(sol, generator_indices, :v, legend = (0.8, 0.7), ylabel="V [p.u.]",label = labels)
    pl_p = plot(sol, generator_indices, :p, legend = (0.8, 0.7), ylabel="p [p.u.]", label=labels)
    pl_q = plot(sol, generator_indices, :q, legend = (0.8, 0.7), ylabel="q [p.u.]", label=labels)
    pl_ω = plot(sol, generator_indices, :ω, legend = (0.8, 0.7), ylabel="ω [rad/s]", label=labels)
    pl_φ = plot(sol, generator_indices, :φ, legend = (0.8, 0.7), ylabel="φ [rad]", label=labels)

    pl = plot(pl_ω, pl_v, pl_p, pl_φ;
            layout=(2,2),
            size = (1000, 500),
            lw=3,
            xlabel="t[s]")

end

