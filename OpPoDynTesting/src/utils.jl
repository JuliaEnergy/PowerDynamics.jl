struct TrajectoriesOfInterest
    sol::SciMLBase.ODESolution
    syms::OrderedDict
end
function plottoi(tois...)
    @assert allequal(getproperty.(tois, :syms)) "All TOI musst have the same symbols of interest."
    toi = tois[1]
    row = 1
    fig = Figure(size=(1000,min(2000, 500*length(toi.syms))))
    for (title, series) in toi.syms
        ax = Axis(fig[row, 1], title=title)
        for (seriesi, (_label, sidx)) in enumerate(series)
            for (toii, toi) in enumerate(tois)
                ts = range(toi.sol.t[begin], toi.sol.t[end], length=1000)
                data = toi.sol(ts, idxs=sidx)
                label =  length(tois) > 1 ? _label*" (timeseries $toii)" : _label
                lines!(ax, data.t, data.u; label,
                    color=Cycled(seriesi),
                    linestyle=ax.scene.theme.palette.linestyle[][toii])
            end
        end
        axislegend(ax)
        row += 1
    end
    fig
end
function similartoi(toi1, toi2; tol=1e-5, verbose=false)
    @assert toi1.syms == toi2.syms "All TOI musst have the same symbols of interest."
    ts = Float64[]
    for toi in (toi1, toi2)
        append!(ts, toi.sol.t)
    end
    unique!(sort!(ts))

    absres = 0
    N = 0
    for (title, series) in toi1.syms
        verbose && println("Comparing $title:")
        for (label, sidx) in series
            data1 = toi1.sol(ts, idxs=sidx)
            data2 = toi2.sol(ts, idxs=sidx)
            res = sum(abs2.(data1.u .- data2.u)) / length(ts)
            verbose && println("  - $label: L²/N error = $res")
            N += 1
            absres += res
        end
    end
    absres = absres/N
    verbose && println("Total L²/N error = $absres")
    absres < tol
end
