const REFDIR = Ref("")
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
serialized and compared to the existing reference if available.
"""
macro reftest(name, toi)
    quote
        _toi = $(esc(toi))
        result = _save_and_compare($(esc(name)), _toi)
        if _intestset()
            @test result
        end
        _toi
    end
end

_intestset() = Test.get_testset() !== Test.FallbackTestSet()

function _save_and_compare(name::String, toi; tol=1e-6)
    @assert !isempty(REFDIR[]) "Please set REFDIR[]!"
    refpath = joinpath(REFDIR[], name*".ser")
    if !isfile(refpath)
        serialize(refpath, toi)
        printstyled("Unknown reference! Stored new file as $(name).ser\n"; color=:yellow)
        return false
    else # file allready exists
        tmppath = joinpath(TMPDIR[], name*".ser")
        serialize(tmppath, toi)
        reftoi = deserialize(refpath)
        diff = compare(toi, reftoi)
        pluspath = replace(refpath, r".ser$" => s"+.ser")
        plotpath = replace(refpath, r".ser$" => s"_comparison.png")
        if diff < tol
            _intestset() || printstyled("Test Passed"; color=:green, bold=true)
            # remove the tmp file
            rm(tmppath)
            # if the test succeds, delete the pluspath if it was there
            if isfile(pluspath)
                printstyled(": Remove conflicting version $(name)+.ser and comparison plot."; color=:green, bold=true)
                rm(pluspath)
            end
            isfile(plotpath) && rm(plotpath)
            _intestset() || println()
            return true
        else
            LAST[] = name
            if isfile(pluspath) && compare(deserialize(pluspath), toi) < tol
                printstyled("Test Failed: Created same $(name)+.ser again! Call `refup()` to accept.\n"; color=:yellow, bold=true)
                rm(tmppath)
            else
                mv(tmppath, pluspath, force=true)
                printstyled("Test Failed: stored new version as $(name)+.ser. Call `refup()` to accept.\n"; color=:red, bold=true)
            end

            if !ismissing(Makie.current_backend())
                plotpath = replace(refpath, r".ser$" => s"_comparison.png")
                try
                    fig = plottoi(toi, reftoi; names=["new", "ref"])
                    Makie.save(plotpath, fig)
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
    old = joinpath(REFDIR[], name*".ser")
    new = joinpath(REFDIR[], name*"+.ser")
    plot = joinpath(REFDIR[], name*"_comparison.png")
    @assert isfile(old) "There is no file for $name.ser"
    @assert isfile(new) "There is no file for $name+.ser"
    rm(old)
    isfile(plot) && rm(plot)
    mv(new, old)
    printstyled("Replaced $name.ser with $name+.ser!\n"; color=:green, bold=true)
end

function refup(s::Symbol)
    if s===:all
        map(readdir(REFDIR[])) do f
            m = match(r"^(.+)\+.ser$", f)
            if !isnothing(m)
                refup(only(m.captures))
            end
        end
    else
        error("Invalid argument $s")
    end
    nothing
end
