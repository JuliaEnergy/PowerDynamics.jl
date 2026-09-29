using Luxor

# `layer` selects which part of the logo gets drawn, so the moving parts can be rendered on their
# own for the animated SVG. `:all` draws everything, `:static` everything that never moves.
draws(layer, part) = layer === :all || layer === part

function blade(thicken)
    height = 100 + thicken
    thickness_bot = 10 + thicken
    thickness_top = 4 + thicken
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
function windturbine(;p=Point(0,0), s=1.0, α=0, layer=:all)
    gsave()
    translate(p)
    scale(s)
    thickness_bot = 22
    thickness_top = 8
    height = 170
    cpthickness = 13
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
    draws(layer, :static) && drawpath(basepath, :fill)
    draws(layer, :rotor) && drawpath(rotorpath, :fill)

    pivot = getworldposition(prot; centered=false)
    grestore()
    pivot
end

function bus(;p, l, r)
    height = 10
    xmin = p.x-l
    ymin = p.y-height/2

    rect(xmin, ymin, l+r, height, :fill)
end

function load(; p, load_scale=1.0, layer=:all)
    length = 90 * load_scale
    arrowwidth = 40
    linewidth = 8
    gsave()
    setline(linewidth)

    pm = p + Point(0, length)
    pm_red = pm - Point(0, linewidth/2)

    pl = pm + Point(-arrowwidth/2, -arrowwidth/2)
    pr = pm + Point(arrowwidth/2, -arrowwidth/2)
    draws(layer, :shaft) && poly([p, pm_red], action=:stroke, close=false)
    draws(layer, :arrow) && poly([pl, pm, pr], action=:stroke, close=false)
    pivot = getworldposition(p; centered=false)
    grestore()
    pivot
end

function gen(; p, α=0, layer=:all)
    offset = 70
    # linewidth = 7.5
    linewidth = 8
    gap = 4
    diameter = 75
    # anchor_width = 30
    # anchor_height = 25
    anchor_width = 25
    anchor_height = 20

    center = p + Point(0, offset)

    gsave()
    setline(linewidth)
    if draws(layer, :static)
        line(p, p + Point(0, offset-diameter/2), :stroke)
        circle(center, diameter/2, :stroke)
    end

    secrad = diameter/2 - linewidth - gap

    translate(center)
    pivot = getworldposition(; centered=false)
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
    draws(layer, :anchor) ? strokepath() : newpath()

    grestore()
    pivot
end

"""
Draws the logo and returns the world positions of the pivots of the moving parts.
"""
function powerdynamics_logo(; p=Point(0,0), s=1.0, foreground="black",
                              wt_phase=-0.2, gen_phase=0.0, load_scale=1.0, layer=:all)
    gsave()

    translate(p)
    scale(s)

    units = 50
    sethue(foreground)
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
    draws(layer, :static) ? strokepath() : newpath()

    sethue(Luxor.julia_green)
    draws(layer, :static) && bus(p=p_wt_con, l=2units, r=.75units)
    rotor = windturbine(; p=p_wt_con - Point(1.25*units,0), s=.75, α=wt_phase, layer)
    # windturbine(; p=p_wt_con - Point(1.25*units,0), s=.75, α=-π/6)

    sethue(Luxor.julia_purple)
    draws(layer, :static) && bus(p=p_load_con, l=.75units, r=2units)
    load_top = load(; p=p_load_con + Point(1.25*units, 0), load_scale, layer)

    sethue(Luxor.julia_red)
    draws(layer, :static) && bus(p=p_gen_conl, l=.75units, r=2.75units)
    anchor = gen(; p=(p_gen_conl + p_gen_conr)/2, α=gen_phase, layer)

    grestore()
    (; rotor, anchor, load_top)
end
