using Gtk.ShortNames
import Graphics

# Global variables
points = []
curve = []
iterations = 0
closed_curve = true

distance(p, q) = norm(p - q)

function regenerate_curve()
    global curve
    curve = points
end


# Graphics

function draw_polygon(ctx, poly)
    if isempty(poly)
        return
    end
    Graphics.new_path(ctx)
    Graphics.move_to(ctx, poly[1][1], poly[1][2])
    for p in poly[2:end]
        Graphics.line_to(ctx, p[1], p[2])
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
    draw_polygon(ctx, points)

    # Generated curve
    Graphics.set_source_rgb(ctx, 0, 0, 1)
    Graphics.set_line_width(ctx, 1.6)
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
        distance(p, q) < 10
    end
    if clicked == 0
        push!(points, p)
        clicked = length(points)
        regenerate_curve()
        draw(canvas)
    end
end

@guarded function mousemove_handler(canvas, event)
    global clicked
    points[clicked] = [event.x, event.y]
    regenerate_curve()
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
        regenerate_curve()
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
        regenerate_curve()
        draw(canvas)
    end
    push!(hbox, closedp)

    # Iterations Spinbutton
    its = @SpinButton(0:10)
    setproperty!(its, :value, iterations)
    signal_connect(its, "value-changed") do sb
        global iterations = getproperty(sb, :value, Int)
        regenerate_curve()
        draw(canvas)
    end
    push!(hbox, @Label("# of iterations:"))
    push!(hbox, its)

    showall(win)
end

setup_gui()
