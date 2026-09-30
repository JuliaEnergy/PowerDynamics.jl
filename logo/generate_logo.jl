using Revise
Revise.includet("logo_utils.jl")

#=
Each image is a draw function which paints onto the current drawing and returns the pivots of the
moving parts. Keyword arguments are passed on to `powerdynamics_logo`, so the same function serves
the static exports, the GIF frames and the layers of the animated SVG.
=#
function draw_logo(; kwargs...)
    origin()
    powerdynamics_logo(; kwargs...)
end

function draw_preview(; foreground="black", layer=:all, kwargs...)
    origin(Point(200,200))
    pivots = powerdynamics_logo(; foreground, layer, kwargs...)
    if draws(layer, :static)
        sethue(foreground)
        fontsize = 120
        textpos = 230
        setfont("TamilMN", fontsize)
        settext("Power", Point(textpos,-fontsize/2); halign="left", valign="center")
        settext("Dynamics.jl", Point(textpos,fontsize/2); halign="left", valign="center")
    end
    pivots
end

function draw_banner(; foreground="black", layer=:all, kwargs...)
    origin(Point(200,200))
    pivots = powerdynamics_logo(; foreground, layer, kwargs...)
    if draws(layer, :static)
        sethue(foreground)
        setfont("TamilMN", 180)
        settext("PowerDynamics.jl", Point(230,0); halign="left", valign="center")
    end
    pivots
end

sizes = Dict(draw_logo => (400, 400), draw_preview => (1150, 400), draw_banner => (2050, 400))

function export_static(pathname, draw; background_color=nothing, kwargs...)
    Drawing(sizes[draw]..., pathname)
    isnothing(background_color) || background(background_color)
    draw(; kwargs...)
    finish()
    pathname
end

function compute_animation_duration(f_wt, f_sg, f_load)
    # Wind turbine has 3-fold symmetry, generator has 2-fold symmetry
    T_wt = rationalize(1 / (3 * f_wt))
    T_sg = rationalize(1 / (2 * f_sg))
    T_load = rationalize(1 / f_load)
    lcm(T_wt, T_sg, T_load)
end

function export_gif(pathname, draw;
        f_wt=-1/8, f_sg=1/4, f_load=1/16,
        load_amplitude=0.3, framerate=30,
        foreground="black", background_color="white")
    T = compute_animation_duration(f_wt, f_sg, f_load)
    nframes = round(Int, T * framerate)

    function frame(scene, framenumber)
        background(background_color)
        t = (framenumber - 1) / framerate
        wt_phase = -0.2 + 2π * f_wt * t
        gen_phase = 2π * f_sg * t
        load_scale = 1.0 + load_amplitude * sin(2π * f_load * t)
        draw(; wt_phase, gen_phase, load_scale, foreground)
    end

    movie = Movie(sizes[draw]..., "logo-animated", 1:nframes)
    animate(movie, [Scene(movie, frame, 1:nframes)];
        creategif=true, pathname, framerate)
end

# Render a single layer with Luxor and return the SVG markup inside the <svg> root.
function layer_markup(draw, layer; foreground)
    Drawing(sizes[draw]..., :svg)
    pivots = draw(; foreground, layer)
    finish()
    body = match(r"<svg[^>]*>(.*)</svg>"s, svgstring()).captures[1]
    strip(body), pivots
end

fmt(x) = string(round(x; digits=3))

# Sampled keyframes for a sine-driven animation. `f(s)` maps the load scale to a transform.
function sine_keyframes(name, f; amplitude, nsamples=32)
    frames = map(0:nsamples) do i
        s = 1.0 + amplitude * sin(2π * i / nsamples)
        "  $(fmt(100 * i / nsamples))% { transform: $(f(s)); }"
    end
    "@keyframes $name {\n" * join(frames, "\n") * "\n}"
end

"""
Write a self-animating SVG. The animation is plain CSS inside the file, so it runs when embedded
via <img> (GitHub READMEs, Documenter) and stops for `prefers-reduced-motion`.
"""
function export_animated_svg(pathname, draw;
        f_wt=-1/8, f_sg=1/4, f_load=1/16, load_amplitude=0.3, foreground="black")
    width, height = sizes[draw]
    static, pivots = layer_markup(draw, :static; foreground)
    rotor, _ = layer_markup(draw, :rotor; foreground)
    anchor, _ = layer_markup(draw, :anchor; foreground)
    shaft, _ = layer_markup(draw, :shaft; foreground)
    arrow, _ = layer_markup(draw, :arrow; foreground)

    # The shaft stretches from the bus downwards and the arrow head rides on its tip.
    # Constants mirror `load` in logo_utils.jl: shaft length 90*s, shortened by half the linewidth.
    shaft_scale(s) = "scaleY($(fmt((90s - 4) / 86)))"
    arrow_shift(s) = "translateY($(fmt(90 * (s - 1)))px)"

    origin_of(p) = "transform-origin: $(fmt(p.x))px $(fmt(p.y))px;"
    spin(name, f) = "animation: $name $(fmt(1 / abs(f)))s linear infinite;"

    css = """
    @keyframes rotor { to { transform: rotate($(sign(f_wt) * 360)deg); } }
    @keyframes anchor { to { transform: rotate($(sign(f_sg) * 360)deg); } }
    $(sine_keyframes("shaft", shaft_scale; amplitude=load_amplitude))
    $(sine_keyframes("arrow", arrow_shift; amplitude=load_amplitude))
    .rotor { $(origin_of(pivots.rotor)) $(spin("rotor", f_wt)) }
    .anchor { $(origin_of(pivots.anchor)) $(spin("anchor", f_sg)) }
    .shaft { $(origin_of(pivots.load_top)) animation: shaft $(fmt(1 / f_load))s linear infinite; }
    .arrow { animation: arrow $(fmt(1 / f_load))s linear infinite; }
    @media (prefers-reduced-motion: reduce) {
      .rotor, .anchor, .shaft, .arrow { animation: none; }
    }
    """

    svg = """
    <svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" width="$width" height="$height" viewBox="0 0 $width $height">
    <style>
    $css</style>
    <g>$static</g>
    <g class="rotor">$rotor</g>
    <g class="anchor">$anchor</g>
    <g class="shaft">$shaft</g>
    <g class="arrow">$arrow</g>
    </svg>
    """
    write(pathname, svg)
    pathname
end

asset(name) = joinpath(@__DIR__, name)

# logo and banner are transparent, preview and gif need an opaque background
for (suffix, foreground, bg) in (("", "black", "white"), ("-dark", "white", "black"))
    for ext in ("svg", "png")
        export_static(asset("logo$suffix.$ext"), draw_logo; foreground)
        export_static(asset("banner$suffix.$ext"), draw_banner; foreground)
    end
    export_static(asset("preview$suffix.png"), draw_preview; foreground, background_color=bg)
    export_animated_svg(asset("logo-animated$suffix.svg"), draw_logo; foreground)
    export_animated_svg(asset("banner-animated$suffix.svg"), draw_banner; foreground)
    export_gif(asset("logo-animated$suffix.gif"), draw_logo; foreground, background_color=bg)
end

# use imagick to convert svg to ico
# run(`convert -background none $(@__DIR__)/logo.svg -define icon:auto-resize=64,48,32,16 $(@__DIR__)/favicon.ico`)
# current logo is not suitable as favicon, far to thin
# maye later somethin like  turbine - gen - load in line?
