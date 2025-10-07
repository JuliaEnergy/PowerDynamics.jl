using Luxor

function blade(thicken)
    height = 100 + thicken
    thickness_bot = 9 + thicken
    thickness_top = 3 + thicken
    move(0, 0)
    line(Point(thickness_bot/2, 0))
    line(Point(thickness_top/2, -height))
    line(Point(-thickness_top/2, -height))
    # line(Point(-thickness_bot/2, 0))

    pstart = currentpoint()
    pstop = Point(-thickness_bot/2, 0)
    p1 = Point(-6-thicken, -0.1*height)
    p2 = Point(-6-thicken, -0.1*height)

    curve(p1, p2, pstop)
    closepath()
end
function rotor(;p=Point(0,0), s=1.0, α=0, thicken=0)
    for i in (0:2) .* 2*π./3 .+ α
        newsubpath()
        gsave()
            translate(p)
            # scale(s)
            rotate(i)
            blade(thicken)
        grestore()
    end
    path = storepath()
    newpath()
    path
end
function windturbine(;p=Point(0,0), s=1.0, α=0)
    gsave()
    translate(p)
    scale(s)
    thickness_bot = 18
    thickness_top = 7
    height = 170
    cpthickness = 10
    cpheight = 0.2*height

    move(0, 0)
    line(Point(thickness_bot/2, 0))
    curve(
        Point(cpthickness/2, -cpheight),
        Point(cpthickness/2, -cpheight),
        Point(thickness_top/2, -height)
    )
    line(Point(-thickness_top/2, -height))
    curve(
        Point(-cpthickness/2, -cpheight),
        Point(-cpthickness/2, -cpheight),
        Point(-thickness_bot/2, 0)
    )
    closepath()

    basepath = storepath()
    newpath()
    prot = Point(0, -height)
    rotorpath = rotor(; p=prot, s, α)
    # rotorpath_clip = rotor(; p=prot, s, α, thicken=3)

    # drawpath(rotorpath_clip, :clip)
    drawpath(basepath, :fill)
    drawpath(rotorpath, :fill)

    grestore()
end

function bus(;p, l, r)
    height = 10
    xmin = p.x-l
    ymin = p.y-height/2

    rect(xmin, ymin, l+r, height, :fill)
end

function load(; p)
    length = 90
    arrowwidth = 40
    linewidth = 6
    gsave()
    setline(linewidth)

    pm = p + Point(0, length)
    pm_red = pm - Point(0, linewidth/2)

    pl = pm + Point(-arrowwidth/2, -arrowwidth/2)
    pr = pm + Point(arrowwidth/2, -arrowwidth/2)
    poly([p, pm_red], action=:stroke, close=false)
    poly([pl, pm, pr], action=:stroke, close=false)
    grestore()

end

function gen(; p, α=0)
    offset = 70
    # linewidth = 7.5
    linewidth = 6
    gap = 3
    diameter = 75
    anchor_width = 30
    anchor_height = 25

    center = p + Point(0, offset)

    gsave()
    setline(linewidth)
    line(p, p + Point(0, offset-diameter/2), :stroke)
    circle(center, diameter/2, :stroke)

    secrad = diameter/2 - linewidth - gap

    translate(center)
    rotate(α)

    # middle part of anchor
    # p_tr = center + Point(anchor_width/2, -anchor_height/2)
    # p_br = center + Point(anchor_width/2, anchor_height/2)
    # p_bl = center + Point(-anchor_width/2, anchor_height/2)
    # p_tl = center + Point(-anchor_width/2, -anchor_height/2)
    p_tr = Point(anchor_width/2, -anchor_height/2)
    p_br = Point(anchor_width/2, anchor_height/2)
    p_bl = Point(-anchor_width/2, anchor_height/2)
    p_tl = Point(-anchor_width/2, -anchor_height/2)

    # anchor arms
    anchor_arm_yoffset = sqrt(secrad^2 - (anchor_width/2)^2)
    # p_ttr = center + Point(anchor_width/2, -anchor_arm_yoffset)
    # p_bbr = center + Point(anchor_width/2,  anchor_arm_yoffset)
    # p_bbl = center + Point(-anchor_width/2,  anchor_arm_yoffset)
    # p_ttl = center + Point(-anchor_width/2, -anchor_arm_yoffset)
    p_ttr = Point(anchor_width/2, -anchor_arm_yoffset)
    p_bbr = Point(anchor_width/2,  anchor_arm_yoffset)
    p_bbl = Point(-anchor_width/2,  anchor_arm_yoffset)
    p_ttl = Point(-anchor_width/2, -anchor_arm_yoffset)

    move(p_br)
    line(p_bbr)
    carc2r(Point(0,0), p_bbr, p_ttr)
    line(p_tr)
    line(p_tl)
    line(p_ttl)
    carc2r(Point(0,0), p_ttl, p_bbl)
    line(p_bl)
    closepath()
    strokepath()

    grestore()
end

function powerdynamics_logo(; p=Point(0,0), s=1.0)
    gsave()

    translate(p)
    scale(s)

    units = 50
    sethue(color)
    setline(5)
    move(-1.2units,0.2units)
    p_wt_con = currentpoint()
    rline(Point(0, -.75*units))
    rline(Point(1*units, -1*units))
    rline(Point(2*units, 0))
    rline(Point(0, .75*units))
    p_load_con = currentpoint()
    rline(Point(0, 0.75*units))
    rline(Point(-1*units, 1*units))
    rline(Point(0, .75*units))
    p_gen_conr = currentpoint()
    rline(Point(-2*units, 0))
    p_gen_conl = currentpoint()
    closepath()
    strokepath()

    sethue(Luxor.julia_green)
    bus(p=p_wt_con, l=2units, r=.75units)
    windturbine(; p=p_wt_con - Point(1.25*units,0), s=.75, α=-0.2)
    # windturbine(; p=p_wt_con - Point(1.25*units,0), s=.75, α=-π/6)

    sethue(Luxor.julia_purple)
    bus(p=p_load_con, l=.75units, r=2units)
    load(p=p_load_con + Point(1.25*units, 0))

    sethue(Luxor.julia_red)
    bus(p=p_gen_conl, l=.75units, r=2.75units)
    gen(p=(p_gen_conl + p_gen_conr)/2, α=0)

    grestore()
end

path = joinpath(dirname(@__DIR__), "docs", "src", "assets", "logo.svg")
path_dark = joinpath(dirname(@__DIR__), "docs", "src", "assets", "logo-dark.svg")

function logo_svg(path, color="black")
    Drawing(400, 400, path)
    origin()
    # sethue(Luxor.julia_blue)
    finish()
    preview()
end
logo_svg(path, "black")
logo_svg(path_dark, "white")
