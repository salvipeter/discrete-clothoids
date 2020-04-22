module Clothoid

using Gtk
import Graphics
using LinearAlgebra

# Parameters
curvature_scaling = -2000
tangent_length = 80

# GUI parameters
closed_curve = false
subsampling = 25
iterations = 100
alpha = -1.0
precise = false
precision = 1.0e-5

# Global variables
points = []
curve = []
indices = []
curvature_comb = []

distance(p, q) = norm(p - q)

function perp2d(p)
    u = [-p[2], p[1]]
    unorm = norm(u)
    unorm > 0 ? u / unorm : u
end

"Like `linspace`, but the last value is not included."
linspace_nolast(start, stop, n) = collect(range(start, stop=stop, length=n+1))[1:n]

"The return value is the discrete curvature at the central point."
function curvature(prev, p, next)
    denom = distance(prev, p) * distance(p, next) * distance(prev, next)
    2 * det(hcat(p - prev, next - p)) / denom
end

"""
Subsamples the given data by `subsampling`, cycling through if the second argument is `true`.
"""
function subsample(original, closedp)
    if length(original) < 2
        return []
    end

    # Compute the average distance & sampling density
    distances = map(norm, original[2:end] - original[1:end-1])
    if closedp
        push!(distances, norm(original[1] - original[end]))
    end
    avg_distance = sum(distances) / length(distances)
    density = subsampling / avg_distance

    n = length(original)
    result = []
    global indices = [1]
    for i in 2:n
        resolution = Int(round(distances[i-1] * density))
        append!(result, linspace_nolast(original[i-1], original[i], resolution))
        push!(indices, length(result) + 1)
    end
    if closedp
        resolution = Int(round(distances[end] * density))
        append!(result, linspace_nolast(original[end], original[1], resolution))
        push!(indices, length(result) + 1)
    else
        push!(result, original[end])
    end

    result
end

function target_curvature(indices, curvatures, closedp)
    result = []
    curvatures = map(x -> sign(x)*abs(x)^(-alpha), curvatures)
    for i in 2:length(curvatures)
        resolution = indices[i] - indices[i-1]
        append!(result, linspace_nolast(curvatures[i-1], curvatures[i], resolution))
    end
    if closedp
        resolution = indices[end] - indices[end-1]
        append!(result, linspace_nolast(curvatures[end], curvatures[1], resolution))
    else
        push!(result, curvatures[end])
    end
    map(x -> sign(x)*abs(x)^(-1/alpha), result)
end

"Updates one point, given its neighbors and the target curvature."
function update(prev, p, next, target)
    # Rearrange the equation
    target *= distance(prev, p) * distance(p, next) * distance(prev, next) / 2.0
    target += prev[1] * next[2] - prev[2] * next[1]
    beta_x, beta_y = next[2] - prev[2], prev[1] - next[1]
    # Now we need: x beta_x + y beta_y = target
    start = (next + prev) / 2
    dir = perp2d(next - prev)
    # [x,y] is of the form start + alpha * dir
    target -= start[1] * beta_x + start[2] * beta_y
    alpha = target / (dir[1] * beta_x + dir[2] * beta_y)
    start + alpha * dir
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

    n = length(points)
    if n < 3
        return
    end

    # Approximate clothoid curve
    curve_curv = []
    iter = 0
    while precise || iter < iterations
        # Generate curvature values at the seed points
        point_curv = []
        if closed_curve
            push!(point_curv, curvature(curve[end], curve[1], curve[2]))
        else
            push!(point_curv, 0.0)
        end
        for i in 2:length(points)-1
            j = indices[i]
            push!(point_curv, curvature(curve[j-1], curve[j], curve[j+1]))
        end
        if closed_curve
            j = indices[n]
            push!(point_curv, curvature(curve[j-1], curve[j], curve[j+1]))
        else
            push!(point_curv, 0.0)
        end

        # Generate target curvature values on the curve
        curve_curv = target_curvature(indices, point_curv, closed_curve)

        # Update curve points
        tmp = similar(curve)
        tmp[1] = curve[1]
        for i in 2:n
            tmp[indices[i]] = curve[indices[i]]
            for j in indices[i-1]+1:indices[i]-1
                tmp[j] = update(curve[j-1], curve[j], curve[j+1], curve_curv[j])
            end
        end
        if closed_curve
            for j in indices[end-1]+1:indices[end]-2
                tmp[j] = update(curve[j-1], curve[j], curve[j+1], curve_curv[j])
            end
            tmp[end] = update(curve[end-1], curve[end], curve[1], curve_curv[end])
        end

        # Find maximal distance and exit when small
        change = maximum(norm, curve - tmp)
        curve = tmp
        if change < precision
            break
        end
        iter += 1
    end

    # Generate curvature comb
    normals = similar(curve)
    for i in eachindex(curve)
        if i == 1
            normals[1] = perp2d(curve[2] - curve[1])
        elseif i == length(curve)
            normals[end] = perp2d(curve[end] - curve[end-1])
        else
            normals[i] = perp2d(curve[i+1] - curve[i-1])
        end
    end
    if isempty(curve_curv)
        curve_curv = zeros(length(curve))
    end
    global curvature_comb = normals .* curve_curv * curvature_scaling
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

draw_callback = @guarded (canvas) -> begin
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

    # Approximate curvature comb
    if length(points) > 2
        Graphics.set_source_rgb(ctx, 0, 1, 0)
        Graphics.set_line_width(ctx, 1.0)
        Graphics.new_path(ctx)
        for i in eachindex(curve)
            Graphics.move_to(ctx, curve[i][1], curve[i][2])
            Graphics.rel_line_to(ctx, curvature_comb[i][1], curvature_comb[i][2])
        end
        Graphics.stroke(ctx)
    end

    if subsampling <= 10
        Graphics.set_source_rgb(ctx, 0.4, 0.4, 0.4)
        for p in curve
            Graphics.arc(ctx, p[1], p[2], 3, 0, 2pi)
            Graphics.fill(ctx)
        end
    end

    # Approximated tangents (does not account for non-equidistant sampling)
    Graphics.set_source_rgb(ctx, 1, 0, 1)
    Graphics.set_line_width(ctx, 1.0)
    for i in 2:length(points)-1
        j = indices[i]
        dir = curve[j-1] - curve[j+1]
        len = norm(dir)
        if len > 0
            dir /= len
        end
        p1 = points[i] - dir * tangent_length / 2
        p2 = points[i] + dir * tangent_length / 2
        Graphics.new_path(ctx)
        Graphics.move_to(ctx, p1[1], p1[2])
        Graphics.line_to(ctx, p2[1], p2[2])
        Graphics.stroke(ctx)
    end

    # Input points
    Graphics.set_source_rgb(ctx, 0, 0, 0)
    for p in points[2:end-1]
        Graphics.arc(ctx, p[1], p[2], 5, 0, 2pi)
        Graphics.fill(ctx)
    end
    if !isempty(points)
        Graphics.set_source_rgb(ctx, 0, 0.8, 0.8)
        for p in (points[1], points[end])
            Graphics.arc(ctx, p[1], p[2], 5, 0, 2pi)
            Graphics.fill(ctx)
        end
    end
end


# GUI

mousedown_handler = @guarded (canvas, event) -> begin
    p = [event.x, event.y]
    global clicked = findfirst(points) do q
        distance(p, q) < 10
    end
    if clicked === nothing
        push!(points, p)
        clicked = length(points)
        generate_curve()
        draw(canvas)
    end
end

mousemove_handler = @guarded (canvas, event) -> begin
    global clicked
    points[clicked] = [event.x, event.y]
    generate_curve()
    draw(canvas)
end

function setup_gui()
    win = GtkWindow("Discrete Clothoid")
    vbox = GtkBox(:v)

    # Canvas widget
    canvas = GtkCanvas(500, 500)
    canvas.mouse.button1press = mousedown_handler
    canvas.mouse.button1motion = mousemove_handler
    draw(draw_callback, canvas)
    push!(win, vbox)
    push!(vbox, canvas)

    # Reset button
    reset = GtkButton("Start Over")
    signal_connect(reset, "clicked") do _
        global points = []
        global curve = []
        draw(canvas)
    end
    hbox = GtkBox(:h)
    set_gtk_property!(hbox, :spacing, 10)
    push!(vbox, hbox)
    push!(hbox, reset)

    # Closed Checkbox
    closedp = GtkCheckButton("Closed curve")
    set_gtk_property!(closedp, :active, closed_curve)
    signal_connect(closedp, "toggled") do cb
        global closed_curve = get_gtk_property(cb, :active, Bool)
        generate_curve()
        draw(canvas)
    end
    push!(hbox, closedp)

    # Subsampling Spinbutton
    sss = GtkSpinButton(5:5:100)
    set_gtk_property!(sss, :value, subsampling)
    signal_connect(sss, "value-changed") do sb
        global subsampling = get_gtk_property(sb, :value, Int)
        generate_curve()
        draw(canvas)
    end
    push!(hbox, GtkLabel("Subsampling:"))
    push!(hbox, sss)

    # Iterations Spinbutton
    its = GtkSpinButton(0:20:1000)
    set_gtk_property!(its, :value, iterations)
    signal_connect(its, "value-changed") do sb
        global iterations = get_gtk_property(sb, :value, Int)
        generate_curve()
        draw(canvas)
    end
    push!(hbox, GtkLabel("# of iterations:"))
    push!(hbox, its)

    # Alpha Choices
    hbox = GtkBox(:h)
    push!(vbox, hbox)
    push!(hbox, GtkLabel("Alpha: "))
    choices = [-2.5, -1.0, -0.8, -0.5, -0.2, -0.1]
    radios = [GtkRadioButton(string(choice)) for choice in choices]
    set_gtk_property!(radios[1], :active, true)
    for r in radios
        set_gtk_property!(r, :group, radios[1])
        signal_connect(r, "toggled") do _
            global alpha = choices[findfirst(rb -> rb === r, radios)]
            generate_curve()
            draw(canvas)
        end
        push!(hbox, r)
    end

    # Precise Checkbox
    precisep = GtkCheckButton("Precision:")
    set_gtk_property!(precisep, :active, precise)
    signal_connect(precisep, "toggled") do cb
        global precise = get_gtk_property(cb, :active, Bool)
        generate_curve()
        draw(canvas)
    end
    push!(hbox, precisep)

    # Precision Spinbutton
    prs = GtkSpinButton(1:1:10)
    set_gtk_property!(prs, :value, Int(round(-log10(precision))))
    signal_connect(prs, "value-changed") do sb
        global precision = exp10(-get_gtk_property(sb, :value, Int))
        generate_curve()
        draw(canvas)
    end
    push!(hbox, prs)

    generate_curve()
    showall(win)
end

run() = begin setup_gui(); nothing end

end # module
