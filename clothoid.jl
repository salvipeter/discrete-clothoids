using Gtk.ShortNames
import Graphics

# Global variables
points = []
curve = []
subsampling = 15
iterations = 20
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
function curvature(prev, p, next)
    denom = norm(p - prev) * norm(next - p) * norm(next - prev)
    2 * det(hcat(p - prev, next - p)) / denom
end

"Subsamples the given data by `subsampling`, cycling through if the second argument is `true`."
function subsample(original, closedp)
    result = []
    for i in 2:length(original)
        append!(result, linspace(original[i-1], original[i], subsampling))
    end
    if closedp
        append!(result, linspace(original[end], original[1], subsampling))
    else
        push!(result, original[end])
    end
    result
end

"Updates one point, given its neighbors and the target curvature."
function update(prev, p, next, target)
    # target *= norm(p - prev) * norm(next - p) * norm(next - prev) / 2.0
    start = (next + prev) / 2
    dir = next - prev
    dir = [-dir[2], dir[1]]
    # the result is of the form start + x * dir
    # dumb search:
    result = start
    minerr = abs(curvature(prev, p, next) - target)
    for x in -1:0.1:1
        q = start + x * dir
        err = abs(curvature(prev, q, next) - target)
        if err < minerr
            err = minerr
            result = q
        end
    end
    result
end

"""
Discrete clothoid curve generation, as in R. Schneider, L. Kobbelt,
_Discrete Fairing of Curves and Surfaces Based on Linear Curvature Distribution_,
In Curve and Surface Design: Saint Malo, pp. 371-380, University Press, 2000.

Uses the indirect approach, with 0 curvature at the ends for open curves.
"""
function generate_curve()
    # Generate initial curve points
    global curve = subsample(points, closed_curve)

    if length(points) < 3
        return
    end

    # Approximate clothoid curve
    for it in 1:iterations
        # Generate curvature values at the seed points
        point_curv = []
        if closed_curve
            push!(point_curv, curvature(curve[end], curve[1], curve[2]))
        else
            push!(point_curv, 0.0)
        end
        for i in 2:length(points)-1
            j = 1 + (i-1) * subsampling
            push!(point_curv, curvature(curve[j-1], curve[j], curve[j+1]))
        end
        if closed_curve
            j = 1 + (length(points)-1) * subsampling
            push!(point_curv, curvature(curve[j-1], curve[j], curve[j+1]))
        else
            push!(point_curv, 0.0)
        end

        # Generate target curvature values on the curve
        curve_curv = subsample(point_curv, closed_curve)

        # Update curve points
        tmp = similar(curve)
        for i in eachindex(curve)
            if i % subsampling == 1
                tmp[i] = curve[i]
            elseif closed_curve && i == length(curve)
                tmp[i] = update(curve[end-1], curve[end], curve[1], curve_curv[end])
            else
                tmp[i] = update(curve[i-1], curve[i], curve[i+1], curve_curv[i])
            end
        end
        curve = tmp
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
    draw_polygon(ctx, curve, closed_curve)

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

    # Subsampling Spinbutton
    sss = @SpinButton(5:5:30)
    setproperty!(sss, :value, subsampling)
    signal_connect(sss, "value-changed") do sb
        global subsampling = getproperty(sb, :value, Int)
        generate_curve()
        draw(canvas)
    end
    push!(hbox, @Label("Subsampling:"))
    push!(hbox, sss)

    # Iterations Spinbutton
    its = @SpinButton(0:10:100)
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
