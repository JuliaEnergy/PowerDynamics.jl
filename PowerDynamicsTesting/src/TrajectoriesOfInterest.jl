"""
    TrajectoriesOfInterest(sol, syms; length=1000)

This is a containertype to wrap a solution together with a nested list of symbols of interest.
It can be used in `@reftest` or with `plottoi(toi...)` to plot and compare.
"""
struct TrajectoriesOfInterest{T<:SciMLBase.ODESolution,S<:OrderedDict,TS}
    sol::T
    syms::S
    ts::TS
end
TrajectoriesOfInterest(sol, syms; length=1000) = TrajectoriesOfInterest(sol, syms, range(sol.t[begin], sol.t[end]; length))

Base.getindex(toi::TrajectoriesOfInterest, plot) = toi.syms[plot]
Base.setindex!(toi::TrajectoriesOfInterest, val, plot) = toi.syms[plot] = val

struct FinalTrajectoriesOfInterest{TS,T}
    ts::TS
    syms::OrderedDict{String, OrderedDict{String, Vector{T}}}
end

Base.getindex(toi::FinalTrajectoriesOfInterest, plot) = toi.syms[plot]
Base.setindex!(toi::FinalTrajectoriesOfInterest, val, plot) = toi.syms[plot] = val

finalizetoi(toi::FinalTrajectoriesOfInterest) = toi
function finalizetoi(toi::TrajectoriesOfInterest)
    (; sol, syms, ts) = toi
    T = eltype(sol)
    _syms = OrderedDict{String, OrderedDict{String, Vector{T}}}()

    for (title, series) in syms
        _syms[title] = OrderedDict{String, Vector{T}}()
        for (_label, sidx) in series
            dat = sol(ts; idxs=sidx)
            _syms[title][_label] = dat.u
        end
    end
    FinalTrajectoriesOfInterest(ts, _syms)
end

function keysequal(toi1, toi2)
   if keys(toi1.syms) == keys(toi2.syms)
      for k in keys(toi1.syms)
          if keys(toi1[k]) != keys(toi2[k])
              return false
          end
      end
      return true
   end
   return false
end

function plottoi(tois...; names=["timeseries $i" for i in 1:length(tois)])
    tois = finalizetoi.(tois)
    allsyms = getproperty.(tois, :syms)
    syms = mergewith(merge, allsyms...)

    row = 1
    fig = Figure(size=(1000,min(2000, 500*length(syms))))
    for (title, series) in syms
        ax = Axis(fig[row, 1], title=title)
        for (seriesi, (_label, sidx)) in enumerate(series)
            for (toii, toi) in enumerate(tois)
                ts = toi.ts
                data = try
                    toi[title][_label]
                catch e
                    if e isa KeyError
                        @warn "Could not extract data for $title - $_label from $(names[toii])"
                        continue
                    else
                        rethrow(e)
                    end
                end
                label =  length(tois) > 1 ? _label*" ("*names[toii]*")" : _label
                lines!(ax, ts, data; label,
                    color=Cycled(seriesi),
                    linestyle=ax.scene.theme.palette.linestyle[][toii])
            end
        end
        axislegend(ax)
        row += 1
    end
    fig
end

function compare(_toi1, _toi2; verbose=false)
    toi1 = finalizetoi(_toi1)
    toi2 = finalizetoi(_toi2)
    if !keysequal(toi1, toi2)
        @warn "TOIs have different keys and are thus incomparable."
        return Inf
    end
    if toi1.ts != toi2.ts
        @warn "Timeseries are define und different timesteps and thus incomparable."
        return Inf
    end
    ts = toi1.ts
    absres = 0
    N = 0
    for (title, series) in toi1.syms
        verbose && println("Comparing $title:")
        for (label, sidx) in series
            data1 = toi1[title][label]
            data2 = toi2[title][label]
            res = sum(abs2.(data1 .- data2)) / length(ts)
            verbose && println("  - $label: L²/N error = $res")
            N += 1
            absres += res
        end
    end
    absres/N
end

function savetoi(path, toi)
    path = contains(path, r"\.jld2$") ? path : path*".jld2"
    ftoi = finalizetoi(toi)
    JLD2.save_object(path, ftoi)
end
function loadtoi(path)
    path = contains(path, r"\.jld2$") ? path : path*".jld2"
    JLD2.load_object(path)
end

function Base.show(io::IO, mime::MIME"text/plain", toi::TrajectoriesOfInterest)
    compact = get(io, :compact, false)::Bool
    print("TrajectoriesOfInterest: ODESolution container with plotspec")
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

function Base.show(io::IO, mime::MIME"text/plain", toi::FinalTrajectoriesOfInterest)
    compact = get(io, :compact, false)::Bool
    print("FinalTrajectoriesOfInterest: Evaluated timeseries for plotspec")
    Nplots = length(toi.syms)
    for (i, plot) in enumerate(keys(toi.syms))
        s = i == Nplots ? "└" : "├"
        print(io, "\n ",s,"─ Plot: $plot")
        Nseries = length(toi.syms[plot])
        for (is, series) in enumerate(toi.syms[plot])
            pre = i == Nplots ? " " : "│"
            s = is == Nseries ? "└" : "├"
            print(io, "\n ", pre, "  ", s, "─ ", series.first)
        end
    end
end
