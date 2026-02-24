using Revise
Revise.includet("logo_utils.jl")

assetpath = joinpath(dirname(@__DIR__), "docs", "src", "assets")

function logo_svg(path, foreground="black")
    Drawing(400, 400, path)
    origin()
    # sethue(Luxor.julia_blue)
    # background("white")
    powerdynamics_logo(; foreground)
    finish()
    preview()
end
logo_svg(joinpath(assetpath, "logo.svg"), "black")
logo_svg(joinpath(assetpath, "logo-dark.svg"), "white")
# use imagick to convert svg to ico
# run(`convert -background none $(assetpath)/logo.svg -define icon:auto-resize=64,48,32,16 $(assetpath)/favicon.ico`)
# current logo is not suitable as favicon, far to thin
# maye later somethin like  turbine - gen - load in line?

function logo_text_svg(path, foreground="black", bg=nothing)
    Drawing(1150, 400, path)
    origin(Point(200,200))
    if bg !== nothing
        background(bg)
    end
    powerdynamics_logo(; s=1, foreground)

    sethue(foreground)
    fontsize = 120
    textpos = 230
    setfont("TamilMN", fontsize)
    # settext("PowerDynamics.jl", Point(200,0); halign="left", valign="center")
    settext("Power", Point(textpos,-fontsize/2); halign="left", valign="center")
    settext("Dynamics.jl", Point(textpos,fontsize/2); halign="left", valign="center")
    finish()
    preview()
end
logo_text_svg(joinpath(assetpath, "preview.png"), "black", "white")
# logo_text_svg(joinpath(assetpath, "preview-dark.png"), "white", "black")

function logo_text_wide_svg(path, foreground="black", bg=nothing)
    Drawing(2050, 400, path)
    origin(Point(200,200))
    if bg !== nothing
        background(bg)
    end
    powerdynamics_logo(; s=1, foreground)

    sethue(foreground)
    fontsize = 180
    textpos = 230
    setfont("TamilMN", fontsize)
    settext("PowerDynamics.jl", Point(textpos,0); halign="left", valign="center")
    finish()
    preview()
end
logo_text_wide_svg(joinpath(assetpath, "banner.png"))
logo_text_wide_svg(joinpath(assetpath, "banner-dark.png"), "white")

function compute_animation_duration(f_wt, f_sg, f_load)
    # Wind turbine has 3-fold symmetry, generator has 2-fold symmetry
    T_wt = rationalize(1 / (3 * f_wt))
    T_sg = rationalize(1 / (2 * f_sg))
    T_load = rationalize(1 / f_load)
    lcm(T_wt, T_sg, T_load)
end

function generate_animated_logo(pathname;
        f_wt=-1/8, f_sg=1/4, f_load=1/16,
        load_amplitude=0.3, framerate=30,
        foreground="black", background_color="white", imgsize=400)
    T = compute_animation_duration(f_wt, f_sg, f_load)
    nframes = round(Int, T * framerate)

    function frame(scene, framenumber)
        background(background_color)
        origin()
        t = (framenumber - 1) / framerate
        wt_phase = -0.2 + 2π * f_wt * t
        gen_phase = 2π * f_sg * t
        load_scale = 1.0 + load_amplitude * sin(2π * f_load * t)
        powerdynamics_logo(; wt_phase, gen_phase, load_scale, foreground)
    end

    movie = Movie(imgsize, imgsize, "logo-animated", 1:nframes)
    animate(movie, [Scene(movie, frame, 1:nframes)];
        creategif=true, pathname, framerate)
end

generate_animated_logo(joinpath(assetpath, "logo-animated.gif");
    foreground="black", background_color="white")
