const REFDIR = Ref(joinpath(pkgdir(PowerDynamics),"assets"))
const TMPDIR = Ref(tempdir())
const THREASHOLD = 200
const LAST = Ref("")

function set_reference_dir(path)
    if !isdir(path)
        @info "Create $path"
        mkdir(path)
    end
    REFDIR[] = path
end
function set_reference_dir(mod::Module, dir="assets")
    path = joinpath(pkgdir(mod), dir)
    set_reference_dir(path)
end

"""
    @reftest name TOI

Mark `TrajectoriesOfInterest` object as reference test. Those objects will be
saved using JLD2 and compared to the existing reference if available.
"""
macro reftest(name, toi)
    quote
        _toi = $(esc(toi))
        if VERSION < v"1.11-"
            @warn "Reference tests require julia 1.11 or higher"
            if _intestset()
                @test_broken false
            end
        else
            result = _save_and_compare($(esc(name)), _toi)
            if _intestset()
                @test result
            end
        end
        _toi
    end
end

_intestset() = Test.get_testset() !== Test.FallbackTestSet()

function _save_and_compare(name::String, toi; tol=1e-6)
    @assert !isempty(REFDIR[]) "Please set REFDIR[]!"
    @assert !isfile(REFDIR[]) "$(REFDIR[]) seems to be a file!"
    if !isdir(REFDIR[])
        @info "Create reference directory $(REFDIR[])"
        mkpath(REFDIR[])
    end

    refpath = joinpath(REFDIR[], name*".jld2")
    if !isfile(refpath)
        savetoi(refpath, toi)
        printstyled("Unknown reference! Stored new file as $(name).jld2\n"; color=:yellow)
        return false
    else # file allready exists
        tmppath = joinpath(TMPDIR[], name*".jld2")
        savetoi(tmppath, toi)
        reftoi = try
            loadtoi(refpath)
        catch e
           @warn "Could not deserialize reference file $refpath" e
           nothing
        end
        diff = isnothing(reftoi) ? Inf : compare(toi, reftoi)
        pluspath = replace(refpath, r"\.jld2$" => s"+.jld2")
        plotpath = replace(refpath, r"\.jld2$" => s"_comparison.png")
        if diff < tol
            _intestset() || printstyled("Test Passed"; color=:green, bold=true)
            # remove the tmp file
            rm(tmppath)
            # if the test succeds, delete the pluspath if it was there
            if isfile(pluspath)
                printstyled(": Remove (now obsolete) version $(name)+.jld2 and comparison plot."; color=:green, bold=true)
                rm(pluspath)
            end
            isfile(plotpath) && rm(plotpath)
            _intestset() || println()
            return true
        else
            LAST[] = name
            if isfile(pluspath)
                try
                    plustoi = loadtoi(pluspath)
                    if compare(plustoi, toi) < tol
                        printstyled("Test Failed: New $(name)+.jld2 matches the existing one! Call `refup()` to accept.\n"; color=:red, bold=true)
                    else
                        printstyled("Test Failed: New $(name)+.jld2 differs from existing one! Call `refup()` to accept.\n"; color=:red, bold=true)
                    end
                catch e
                    printstyled("Test Failed: stored new version as $(name)+.jld2, existing one could not be read. Call `refup()` to accept.\n"; color=:red, bold=true)
                end
            else
                printstyled("Test Failed: stored new version as $(name)+.jld2. Call `refup()` to accept.\n"; color=:red, bold=true)
            end
            mv(tmppath, pluspath, force=true)

            if !ismissing(Makie.current_backend())
                plotpath = replace(refpath, r"\.jld2$" => s"_comparison.png")
                try
                    fig = plottoi(toi, reftoi; names=["new", "ref"])
                    Makie.save(plotpath, fig)
                    println("Saved comparison plot $(name)_comparison.png.")
                catch e
                    @warn "Comparison plot failed" e
                end
            end

            return false
        end
    end
end

refup() = isempty(LAST[]) ? error("Don't now which reference to update.") : refup(LAST[])

function refup(name::AbstractString)
    @assert !isempty(REFDIR[]) "Please set REFDIR[]!"
    old = joinpath(REFDIR[], name*".jld2")
    new = joinpath(REFDIR[], name*"+.jld2")
    plot = joinpath(REFDIR[], name*"_comparison.png")
    @assert isfile(old) "There is no file for $name.jld2"
    @assert isfile(new) "There is no file for $name+.jld2"
    rm(old)
    isfile(plot) && rm(plot)
    mv(new, old)
    printstyled("Replaced $name.jld2 with $name+.jls!\n"; color=:green, bold=true)
end

function refup(s::Symbol)
    if s===:all
        map(readdir(REFDIR[])) do f
            m = match(r"^(.+)\+.jld2$", f)
            if !isnothing(m)
                refup(only(m.captures))
            end
        end
    else
        error("Invalid argument $s")
    end
    nothing
end
