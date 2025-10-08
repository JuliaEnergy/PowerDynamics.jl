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
