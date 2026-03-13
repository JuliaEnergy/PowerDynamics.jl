function freep(sys::System)
    return filter(p -> !haskey(ModelingToolkitBase.defaults(sys), p), ModelingToolkitBase.parameters(sys))
end
function freep(cf::NetworkDynamics.ComponentModel)
    return filter(p -> !NetworkDynamics.has_default(cf, p), NetworkDynamics.psym(cf))
end
function freep(nw::Network)
    setdiff(NetworkDynamics.SII.parameter_symbols(nw), keys(NetworkDynamics.SII.default_values(nw)))
end
function freeu(nw::Network)
    setdiff(NetworkDynamics.SII.variable_symbols(nw), keys(NetworkDynamics.SII.default_values(nw)))
end
function freeu(cf::NetworkDynamics.ComponentModel)
    return filter(u -> !NetworkDynamics.has_default(cf, u), NetworkDynamics.sym(cf))
end

set_voltage!(cf::VertexModel; mag, arg) = set_voltage!(cf, mag * exp(im * arg))
function set_voltage!(cf::VertexModel, c::Complex)
    set_default!(cf, :busbar₊u_r, real(c))
    set_default!(cf, :busbar₊u_i, imag(c))
    c
end

get_voltage(cf::VertexModel) = get_default(cf, :busbar₊u_r) + im * get_default(cf, :busbar₊u_i)

function set_current!(cf::VertexModel; P, Q)
    @assert has_default(cf, :busbar₊u_r) && has_default(cf, :busbar₊u_i)
    u = get_default(cf, :busbar₊u_r) + im * get_default(cf, :busbar₊u_i)
    i = conj((P + im * Q) / u)
    set_current!(cf, -i)
end
# set_current!(cf::VertexModel; mag, arg) = set_current!(cf, mag * exp(im * arg))
function set_current!(cf::VertexModel, c::Number)
    set_default!(cf, :busbar₊i_r, real(c))
    set_default!(cf, :busbar₊i_i, imag(c))
    c
end

get_current(cf::VertexModel) = -1 * (get_default(cf, :busbar₊i_r) + im * get_default(cf, :busbar₊i_i))
get_power(cf::VertexModel) = get_voltage(cf) * conj(get_current(cf))

function normalize_angle(θ)
    θ = mod2pi(θ)
    θ >= π ? θ - 2π : θ
end

"""
    load_pdtesting()

Internal function to dynamicially load the PowerDynamicsTesting code into Main.
If called in an interactive context and Revise.jl is available, it will use Revise.includet to load the code.
"""
function load_pdtesting(force=false)
    if !force && isdefined(Main, :PowerDynamicsTesting)
        println("PowerDynamicsTesting already loaded into Main")
        return
    end
    path = pdtesting_path()
    @eval Main begin
        if isinteractive() && isdefined(Main, :Revise)
            println("Loading PowerDynamicsTesting into Main with Revise")
            Revise.includet($path)
            dir = dirname($path)
            for file in readdir(joinpath(dir))
                if endswith(file, ".jl") && file != "PowerDynamicsTesting.jl"
                    path = PowerDynamics.pdtesting_path()
                    println(" -> track $file")
                    Revise.track(Main.PowerDynamicsTesting, joinpath(dir, file))
                end
            end
        else
            println("Loading PowerDynamicsTesting into Main")
            include($path)
        end
    end
end
pdtesting_path() = joinpath(pkgdir(@__MODULE__), "test", "PowerDynamicsTesting", "PowerDynamicsTesting.jl")

"""
    refine_timeseries(ts, factor=10)

Refine a time series by interpolating additional time points between existing ones.
Usefull for creating denser plots based on `refine_timeseries(sol.t)`.

# Arguments
- `ts`: Input time series vector
- `factor`: Number of subdivisions between each pair of consecutive time points (default: 10)
"""
function refine_timeseries(ts, factor=10)
    newts = Float64[]
    for i in 1:length(ts)-1
        t1 = ts[i]
        t2 = ts[i+1]
        push!(newts, t1)
        t1 == t2 && continue
        for j in 1:factor-1
            push!(newts, t1 + j*(t2 - t1)/factor)
        end
    end
    push!(newts, ts[end])
end


"""
    unwrap_deg(angles)

Unwrap phase angles in degrees. Detects and corrects discontinuities
greater than 180° to produce continuous phase trajectories.

Useful for Bode plots and other frequency response visualizations.
"""
unwrap_deg(angles) = _unwrap(angles, 360.0)

"""
    unwrap_rad(angles)

Unwrap phase angles in radians. Detects and corrects discontinuities
greater than π to produce continuous phase trajectories.
"""
unwrap_rad(angles) = _unwrap(angles, 2π)

function _unwrap(angles, range)
    unwrapped = similar(angles)
    unwrapped[1] = angles[1]

    for i in 2:length(angles)
        diff = angles[i] - unwrapped[i-1]
        # Round to nearest multiple of range and subtract to unwrap
        correction = round(diff / range) * range
        unwrapped[i] = angles[i] - correction
    end

    unwrapped
end

const GITHUB_REPO = "JuliaEnergy/PowerDynamics.jl"

# GITHUB_REF is evaluated at *precompile time* and baked into the compiled
# module.  ref_source_file uses it to embed source links in docstrings.
# Intended values per CI context:
#   - locally          → "v{version}" read from Project.toml
#   - PR build         → commit SHA (GITHUB_REF_NAME is "{number}/merge")
#   - push to main     → "main"
#   - tag push         → the tag name, e.g. "v4.4.0"
# IMPORTANT: the docs CI job must wipe the PowerDynamics precompile cache after
# restoring the depot cache and before building, so that each run re-precompiles
# with the correct GITHUB_REF_NAME rather than reusing a stale cached value.
const GITHUB_REF = let
    if haskey(ENV, "GITHUB_REF_NAME")
        ref_name = ENV["GITHUB_REF_NAME"]
        # PR refs look like "253/merge" — not valid for blob links; use commit SHA instead
        if occursin(r"^\d+/merge$", ref_name)
            get(ENV, "GITHUB_SHA", "main")
        else
            ref_name
        end
    else
        projecttoml = read(joinpath(pkgdir(PowerDynamics), "Project.toml"), String)
        versionstring = match(r"version\s?=\s?\"(.*)\"", projecttoml)[1]
        "v" * string(VersionNumber(versionstring))
    end
end
"""
    ref_source_file(f, line)

Generate a documentation string that points readers to the model's source code.
Intended to be called via string interpolation in a docstring, e.g.:

    \$(PowerDynamics.ref_source_file(@__FILE__, @__LINE__))

The GitHub ref embedded in the link is determined by [`GITHUB_REF`](@ref) which
is fixed at precompile time (see its definition for the per-context mapping):

- **locally**: `v{version}` from `Project.toml`
- **PR build**: commit SHA
- **push to `main`**: `"main"` (dev docs, always up to date)
- **tag push**: tag name (e.g. `"v4.4.0"`), so links stay valid in archived versions
"""
function ref_source_file(f, line)
    subf = replace(relpath(f, pkgdir(@__MODULE__)), '\\' => '/')

    link = "https://github.com/$GITHUB_REPO/blob/$GITHUB_REF/$subf#L$(line+2)"
    doctext = "For a concrete list of variables and parameters please check the model source"
    if haskey(ENV, "GITHUB_ACTIONS")
        # for online docs
        doctext *= " on [GitHub]($link)."
    else
        doctext *= "\n - online [$link]($link)"
        doctext *= "\n - local @ $f:$line"
    end
    doctext
end
