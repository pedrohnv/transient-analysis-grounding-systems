# -*- coding: utf-8 -*-
"""
Script to create a set of lines and break them at their intersection points.

Author: Pedro Henrique Nascimento Vieira
"""

using Plots
using LinearAlgebra

struct Point
    x::Real
    y::Real
end

struct Line
    start_point::Point
    end_point::Point
end

function crosses_point(line::Line, point::Point)
    """ Returns True if the point is upon the line. """
    if line.start_point.x <= line.end_point.x
        sx = line.start_point.x
        ex = line.end_point.x
    else
        ex = line.start_point.x
        sx = line.end_point.x
    end
    if line.start_point.y <= line.end_point.y
        sy = line.start_point.y
        ey = line.end_point.y
    else
        ey = line.start_point.y
        sy = line.end_point.y
    end
    xbound = (sx <= point.x <= ex)
    ybound = (sy <= point.y <= ey)
    return (xbound && ybound)
end

function intersection(line1::Line, line2::Line)
    """
    Returns the point at which the lines intersect, provided that they cross it
    Returns nothing otherwise or if there are infinite intersection points.
    """
    a1 = line1.end_point.y - line1.start_point.y
    b1 = line1.start_point.x - line1.end_point.x
    c1 = a1*line1.start_point.x + b1*line1.start_point.y
    a2 = line2.end_point.y - line2.start_point.y
    b2 = line2.start_point.x - line2.end_point.x
    c2 = a2*line2.start_point.x + b2*line2.start_point.y
    delta = a1*b2 - a2*b1
    x = (b2*c1 - b1*c2)
    y = (a1*c2 - a2*c1)
    point = Point(x/delta, y/delta)
    if crosses_point(line1, point) && crosses_point(line2, point)
        return point
    else
        return nothing
    end
end

function break_intersections(lines)
    """
    Breaks the lines in the set in every intersection point and return
    the new set.
    """
    function loop_set(lines)
        for li in lines  # TODO better search than linear?
            for lk in setdiff(lines, [li])
                intersect_point = intersection(li, lk)
                if intersect_point != nothing
                    old_length = new_length = length(lines)
                    for lm in [li, lk]
                        xp = intersect_point.x
                        yp = intersect_point.y
                        x0 = lm.start_point.x
                        y0 = lm.start_point.y
                        x1 = lm.end_point.x
                        y1 = lm.end_point.y
                        c1 = isapprox([xp, yp], [x0, y0])
                        c2 = isapprox([xp, yp], [x1, y1])
                        if !(c1 || c2)
                            x = lm.start_point.x
                            y = lm.start_point.y
                            l1 = Line(Point(x, y), intersect_point)
                            x = lm.end_point.x
                            y = lm.end_point.y
                            l2 = Line(intersect_point, Point(x, y))
                            setdiff!(lines, [lm])
                            union!(lines, [l1, l2])
                            new_length = length(lines)
                        end
                    end
                    if new_length > old_length
                        return lines
                    end
                end
            end
        end
        return lines
    end
    old_length = new_length = length(lines)
    first_run = true
    while new_length > old_length || first_run
        first_run = false
        old_length = new_length
        lines = loop_set(lines)
        new_length = length(lines)
    end
    return lines
end

function plot_lines(lines, nodes=false)
    fig = plot(leg=false, aspect_ratio=1)#, border=:none)
    for l in lines
        x0 = l.start_point.x
        y0 = l.start_point.y
        x1 = l.end_point.x
        y1 = l.end_point.y
        plot!([x0, x1], [y0, y1], line=(1, :black, :solid))
        nodes && scatter!([x0, x1], [y0, y1], markercolor=:black, markersize=2)
    end
    return fig
end

function segment_lines(lines, lmax)
    seg_lines = [Line(Point(0,0), Point(0,0))]
    for li in lines
        x0 = li.start_point.x
        y0 = li.start_point.y
        x1 = li.end_point.x
        y1 = li.end_point.y
        v = [x1 - x0, y1 - y0]
        L = norm(v)
        n = Int(cld(L, lmax))
        dl = v./n
        for k = 0:(n-1)
            p0 = [x0, y0] + k.*dl;
            p1 = p0 .+ dl;
            push!(seg_lines, Line(Point(p0[1], p0[2]), Point(p1[1], p1[2])))
        end
    end
    return seg_lines[2:end]
end

function main(nodes, lmax)
    lines = create_grid()
    seg_lines = segment_lines(lines, lmax)
    if nodes
        fig = plot_lines(seg_lines, nodes)
    else
        fig = plot_lines(lines, nodes)
    end
    savefig(fig, "show_case_grid.png")
    r = 10e-3  # radius
    z = 0.0
    println("number of lines = ", length(seg_lines))
    io = open("show_case_grid.txt", "w")
    for li in seg_lines
        x0 = li.start_point.x
        y0 = li.start_point.y
        x1 = li.end_point.x
        y1 = li.end_point.y
        write(io, join([string(x0), ", ", string(y0), ", ", string(z), ", ",
                        string(x1), ", ", string(y1), ", ", string(z), ", ",
                        string(r), "\n"]))
    end
    close(io)
    return lines, fig
end

function create_grid()
    lines = [
         Line(Point(0, 0), Point(170, 0)),
         Line(Point(0, 0), Point(0, 140)),
         Line(Point(0, 140), Point(90, 140)),
         Line(Point(90, 140), Point(90, 75)),
         Line(Point(85, 80), Point(170, 80)),
         Line(Point(170, 80), Point(170, 0)),
         Line(Point(0, 5), Point(170, 5)),
         Line(Point(5, 0), Point(5, 140)),
         Line(Point(0, 135), Point(90, 135)),
         Line(Point(85, 140), Point(85, 75)),
         Line(Point(85, 75), Point(170, 75)),
         Line(Point(165, 80), Point(165, 0)),
         Line(Point(15, 5), Point(15, 0)),
         Line(Point(35, 5), Point(35, 0)),
         Line(Point(55, 5), Point(55, 0)),
         Line(Point(75, 5), Point(75, 0)),
         Line(Point(95, 5), Point(95, 0)),
         Line(Point(115, 5), Point(115, 0)),
         Line(Point(135, 5), Point(135, 0)),
         Line(Point(155, 5), Point(155, 0)),
         Line(Point(15, 135), Point(15, 140)),
         Line(Point(35, 135), Point(35, 140)),
         Line(Point(55, 135), Point(55, 140)),
         Line(Point(75, 135), Point(75, 140)),
         Line(Point(105, 80), Point(105, 75)),
         Line(Point(120, 80), Point(120, 75)),
         Line(Point(140, 80), Point(140, 75)),
         Line(Point(155, 80), Point(155, 75)),
         Line(Point(0, 15), Point(5, 15)),
         Line(Point(0, 30), Point(5, 30)),
         Line(Point(0, 50), Point(5, 50)),
         Line(Point(0, 70), Point(5, 70)),
         Line(Point(0, 90), Point(5, 90)),
         Line(Point(0, 110), Point(5, 110)),
         Line(Point(0, 125), Point(5, 125)),
         Line(Point(85, 125), Point(90, 125)),
         Line(Point(85, 107.5), Point(90, 107.5)),
         Line(Point(85, 90), Point(90, 90)),
         Line(Point(165, 17.5), Point(170, 17.5)),
         Line(Point(165, 32.5), Point(170, 32.5)),
         Line(Point(165, 47.5), Point(170, 47.5)),
         Line(Point(165, 62.5), Point(170, 62.5)),
         Line(Point(5, 25), Point(165, 25)),
         Line(Point(5, 55), Point(165, 55)),
         Line(Point(25, 5), Point(25, 135)),
         Line(Point(65, 5), Point(65, 135)),
         Line(Point(5, 115), Point(85, 115)),
         Line(Point(5, 85), Point(85, 85)),
         Line(Point(145, 75), Point(145, 5)),
         Line(Point(105, 75), Point(105, 5)),
    ]
    return break_intersections(lines)
end

gr()
@time lines, fig = main(false, 1.0);
