using Gtk.ShortNames
import Graphics

# Global variables
points = []
curve = []
dilation = 5
iterations = 3
closed_curve = true

"Like `linspace` for arrays, but the last value is not included."
function Base.linspace{T<:AbstractFloat}(start :: Array{T}, stop :: Array{T}, n)
    a = Array(Array{T}, n)
    for i = 0:(n-1)
        alpha = i / n
        a[i+1] = start * (1 - alpha) + stop * alpha
    end
    a
end

"The return value is the discrete curvature at the central point."
function curvature(a, b, c)
    denom = norm(b - a) * norm(c - b) * norm(c - a)
    2 * det(hcat(b - a, c - b)) / denom
end

"""
Discrete clothoid curve generation, as in R. Schneider, L. Kobbelt,
_Discrete Fairing of Curves and Surfaces Based on Linear Curvature Distribution_,
Defence Technical Information Center, 2000.

Uses the indirect approach, with 0 curvature at the ends for open curves.
"""
function generate_curve()
    # Generate initial curve points
    global curve = []
    for i in 2:length(points)
        append!(curve, linspace(points[i-1], points[i], dilation))
    end
    if closed_curve
        append!(curve, linspace(points[end], points[1], dilation))
    end

    if length(points) < 3
        return
    end

    # Generate curvature values on the seed points
    point_curv = []
    if closed_curve
        push!(point_curv, curvature(points[end],points[1],points[2]))
    else
        push!(point_curv, 0.0)
    end
    for i in 2:length(points)-1
        push!(point_curv, curvature(points[i-1], points[i], points[i+1]))
    end
    if closed_curve
        push!(point_curv, curvature(points[end-1],points[end],points[1]))
    else
        push!(point_curv, 0.0)
    end

    # Generate target curvature values on the curve
    curve_curv = []
    for i in 2:length(point_curv)
        append!(curve_curv, linspace(point_curv[i-1], point_curv[i], dilation))
    end
    if closed_curve
        append!(curve_curv, linspace(point_curv[end], point_curv[1], dilation))
    end

    # Approximate clothoid curve
    for i in 1:iterations
        # TODO
    end
end


# Graphics

function draw_polygon(ctx, poly, closep = false)
    if isempty(poly)
        return
    end
    Graphics.new_path(ctx)
    Graphics.move_to(ctx, poly[1][1], poly[1][2])
    for p in poly[2:end]
        Graphics.line_to(ctx, p[1], p[2])
    end
    if closep && length(poly) > 2
        Graphics.line_to(ctx, poly[1][1], poly[1][2])
    end
    Graphics.stroke(ctx)
end

@guarded function draw_callback(canvas)
    ctx = Graphics.getgc(canvas)

    # White background
    Graphics.rectangle(ctx, 0, 0, Graphics.width(canvas), Graphics.height(canvas))
    Graphics.set_source_rgb(ctx, 1, 1, 1)
    Graphics.fill(ctx)

    # Input polygon
    Graphics.set_source_rgb(ctx, 1, 0, 0)
    Graphics.set_line_width(ctx, 1.0)
    draw_polygon(ctx, points, closed_curve)

    # Generated curve
    Graphics.set_source_rgb(ctx, 0, 0, 1)
    Graphics.set_line_width(ctx, 2.0)
    draw_polygon(ctx, curve)

    # Input points
    Graphics.set_source_rgb(ctx, 0, 0, 0)
    for p in points
        Graphics.arc(ctx, p[1], p[2], 5, 0, 2pi)
        Graphics.fill(ctx)
    end
end


# GUI

@guarded function mousedown_handler(canvas, event)
    p = [event.x, event.y]
    global clicked = findfirst(points) do q
        norm(p - q) < 10
    end
    if clicked == 0
        push!(points, p)
        clicked = length(points)
        generate_curve()
        draw(canvas)
    end
end

@guarded function mousemove_handler(canvas, event)
    global clicked
    points[clicked] = [event.x, event.y]
    generate_curve()
    draw(canvas)
end

function setup_gui()
    win = @Window("Discrete Clothoid")
    vbox = @Box(:v)

    # Canvas widget
    canvas = @Canvas(500, 500)
    canvas.mouse.button1press = mousedown_handler
    canvas.mouse.button1motion = mousemove_handler
    draw(draw_callback, canvas)
    push!(win, vbox)
    push!(vbox, canvas)

    # Reset button
    reset = @Button("Start Over")
    signal_connect(reset, "clicked") do _
        global points = []
        generate_curve()
        draw(canvas)
    end
    hbox = @Box(:h)
    setproperty!(hbox, :spacing, 10)
    push!(vbox, hbox)
    push!(hbox, reset)

    # Closed Checkbox
    closedp = @CheckButton("Closed curve")
    setproperty!(closedp, :active, closed_curve)
    signal_connect(closedp, "toggled") do cb
        global closed_curve = getproperty(cb, :active, Bool)
        generate_curve()
        draw(canvas)
    end
    push!(hbox, closedp)

    # Dilation Spinbutton
    dil = @SpinButton(5:5:30)
    setproperty!(dil, :value, dilation)
    signal_connect(dil, "value-changed") do sb
        global dilation = getproperty(sb, :value, Int)
        generate_curve()
        draw(canvas)
    end
    push!(hbox, @Label("Dilation:"))
    push!(hbox, dil)

    # Iterations Spinbutton
    its = @SpinButton(0:10)
    setproperty!(its, :value, iterations)
    signal_connect(its, "value-changed") do sb
        global iterations = getproperty(sb, :value, Int)
        generate_curve()
        draw(canvas)
    end
    push!(hbox, @Label("# of iterations:"))
    push!(hbox, its)

    showall(win)
end

setup_gui()
