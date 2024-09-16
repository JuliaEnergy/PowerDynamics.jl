struct TrajectoriesOfInterest{T<:SciMLBase.ODESolution,S<:OrderedDict}
    sol::T
    syms::S
end
Base.getindex(toi::TrajectoriesOfInterest, plot) = toi.syms[plot]
Base.setindex!(toi::TrajectoriesOfInterest, val, plot) = toi.syms[plot] = val

function plottoi(tois...; names=["timeseries $i" for i in 1:length(tois)])
    # @assert allequal(getproperty.(tois, :syms)) "All TOI musst have the same symbols of interest."
    allsyms = getproperty.(tois, :syms)
    syms = mergewith(merge, allsyms...)
    toi = tois[1]
    row = 1
    fig = Figure(size=(1000,min(2000, 500*length(syms))))
    for (title, series) in syms
        ax = Axis(fig[row, 1], title=title)
        for (seriesi, (_label, sidx)) in enumerate(series)
            for (toii, toi) in enumerate(tois)
                ts = range(toi.sol.t[begin], toi.sol.t[end], length=1000)
                data = try
                    toi.sol(ts, idxs=sidx)
                catch e
                    @warn "Could not extract data for $title: $e from series $(names[toii])"
                    continue
                end
                label =  length(tois) > 1 ? _label*" ("*names[toii]*")" : _label
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

function compare(toi1, toi2; verbose=false)
    if isnothing(toi1) || isnothing(toi2)
        return Inf
    end
    if toi1.syms != toi2.syms
        @warn "TOIs have different symbols of interest and are thus incomparable."
        return Inf
    end
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
    absres/N
end

function Base.show(io::IO, mime::MIME"text/plain", toi::TrajectoriesOfInterest)
    compact = get(io, :compact, false)::Bool
    print("ODESolution with trajectories of interest:")
    Nplots = length(toi.syms)
    for (i, plot) in enumerate(keys(toi.syms))
        s = i == Nplots ? "└" : "├"
        print(io, "\n ",s,"─ Plot: $plot")
        Nseries = length(toi.syms[plot])
        for (is, series) in enumerate(toi.syms[plot])
            pre = i == Nplots ? " " : "│"
            s = is == Nseries ? "└" : "├"
            print(io, "\n ", pre, "  ", s, "─ ", series.second, " (", series.first, ")")
        end
    end
end
