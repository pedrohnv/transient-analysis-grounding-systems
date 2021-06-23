# Plot the results from example visacro57emc01
using Plots
using LinearAlgebra
using DataFrames
using CSV
import Printf.@sprintf

font1 = Plots.font("DejaVu Sans", 14)
font2 = Plots.font("DejaVu Sans", 14)
pyplot(guidefont=font1, xtickfont=font1, ytickfont=font1, legendfont=font2);

nomes = ["corner fast", "center fast", "corner slow", "center slow"];

function heidler(t, i0, tau1, tau2, n)
    eta = exp(-tau1 / tau2 * (n * tau2 / tau1)^(1/n))
    tt1n = (t / tau1)^n
    num = i0 * tt1n * exp(-t / tau2)
    den = eta * (1 + tt1n)
    return num / den
end

## Grid
mutable struct Electrode
    """ Defines a conductor segment. """
    start_point::Array{Float64,1}
    end_point::Array{Float64,1}
    middle_point::Array{Float64,1}
    length::Float64
    radius::Float64
end

import Base.==

function ==(sender::Electrode, receiver::Electrode)
    c1 = (sender.start_point == receiver.start_point)
    c2 = (sender.end_point == receiver.end_point)
    return c1 && c2
end

function new_electrode(start_point, end_point, radius)
    """ Creates a conductor segment. """
    return Electrode(start_point, end_point, (start_point + end_point)/2.0,
                     norm(start_point - end_point), radius)
end

function segment_electrode(electrode::Electrode, num_segments::Int)
    """ Segments a conductor. """
    nn = num_segments + 1;
    nodes = Array{Float64}(undef, nn, 3); # FIXME transpose all nodes
    startp = Array{Float64,1}(undef, 3);
    endp = Array{Float64,1}(undef, 3);
    for k = 1:3
        startp[k] = electrode.start_point[k];
        endp[k] = electrode.end_point[k];
    end
    increment = (endp - startp)/num_segments;
    for k = 0:num_segments
        nodes[k+1,:] = startp + k*increment;
    end
    segments = Array{Electrode,1}(undef, num_segments);
    for k = 1:num_segments
        segments[k] = new_electrode(nodes[k,:], nodes[k+1,:], electrode.radius);
    end
    return segments, nodes
end

function matchrow(a, B, atol=1e-9, rtol=0)
    """
    Returns the row in B that matches a. If there is no match, nothing is returned.
    taken and modified from https://stackoverflow.com/a/32740306/6152534
    """
    findfirst(i -> all(j -> isapprox(a[j], B[i,j], atol=atol, rtol=rtol),
                       1:size(B,2)), 1:size(B,1))
end

function matchcol(a, B, atol=1e-9, rtol=0)
    """
    Returns the col in B that matches a. If there is no match, nothing is returned.
    taken and modified from https://stackoverflow.com/a/32740306/6152534
    """
    findfirst(i -> all(j -> isapprox(a[j], B[j,i], atol=atol, rtol=rtol),
                       1:size(B,1)), 1:size(B,2))
end

function seg_electrode_list(electrodes, frac)
    """
    Segments a list of conductors such that they end up having at most 'L/frac'
    length.
    Return a list of the segmented conductors and their nodes.
    """
    num_elec = 0; #after segmentation
    for i=1:length(electrodes)
        #TODO store in array to avoid repeated calculations
        num_elec += Int(ceil(electrodes[i].length/frac));
    end
    elecs = Array{Electrode}(undef, num_elec);
    nodes = zeros(Float64, (2*num_elec, 3));
    e = 1;
    nodes = [];
    for i = 1:length(electrodes)
        ns = Int(ceil(electrodes[i].length/frac));
        new_elecs, new_nodes = segment_electrode(electrodes[i], ns);
        for k=1:ns
            elecs[e] = new_elecs[k];
            e += 1;
        end
        if (nodes == [])
            nodes = new_nodes;
        else
            for k = 1:size(new_nodes)[1]
                if (matchrow(new_nodes[k:k,:], nodes) == nothing)
                    nodes = cat(nodes, new_nodes[k:k,:], dims=1);
                end
            end
        end
    end
    return elecs, nodes
end

function nodes_from_elecs(electrodes)
    """ Given a list of electrodes, return their unique nodes. """
    ns = length(electrodes)
    nodes = rand(2ns, 3)
    n = 0
    for e in electrodes
        i = matchrow(e.start_point, nodes)
        if i == nothing
            n += 1
            nodes[n, :] = e.start_point
        end
        i = matchrow(e.end_point, nodes)
        if i == nothing
            n += 1
            nodes[n, :] = e.end_point
        end
    end
    return nodes[1:n, :]
end

struct Grid
    vertices_x::Int
    vertices_y::Int
    length_x::Float64
    length_y::Float64
    edge_segments_x::Int
    edge_segments_y::Int
    radius::Float64
    depth::Float64
end

function electrode_grid(grid)
    N = grid.edge_segments_x;
    Lx = grid.length_x;
    vx = grid.vertices_x;
    lx = Lx/(N*(vx - 1));
    M = grid.edge_segments_y;
    Ly = grid.length_y;
    vy = grid.vertices_y;
    ly = Ly/(M*(vy - 1));
    num_seg_horizontal = N*vy*(vx - 1);
    num_seg_vertical = M*vx*(vy - 1);
    num_seg = num_seg_horizontal + num_seg_vertical;

    num_elec = N*vy*(vx - 1) + M*vx*(vy - 1);
    num_nodes = vx*vy + vx*(vy - 1)*(M - 1) + vy*(vx - 1)*(N - 1);
    electrodes = Array{Electrode}(undef, num_elec);
    nodes = Array{Float64}(undef, 3, num_nodes);
    nd = 1;
    ed = 1;
    #Ns_hor = N*vy*(vx - 1);
    # Make horizontal electrodes
    for h = 1:vy
        for n = 1:(vx - 1)
            for k = 1:N
                #ed = (h*(vx - 1) + n)*N + k;
                x0 = lx*(N*(n - 1) + k - 1);
                y0 = ly*M*(h - 1);
                start_point = [x0, y0, grid.depth];
                end_point = [x0 + lx, y0, grid.depth];
                electrodes[ed] = new_electrode(start_point, end_point, grid.radius);
                ed += 1;
                if (n == 1 && k == 1)
                    nodes[1, nd] = start_point[1];
                    nodes[2, nd] = start_point[2];
                    nodes[3, nd] = start_point[3];
                    nd += 1;
                end
                nodes[1, nd] = end_point[1];
                nodes[2, nd] = end_point[2];
                nodes[3, nd] = end_point[3];
                nd += 1;
            end
        end
    end
    # Make vertical electrodes
    for g = 1:vx
        for m = 1:(vy - 1)
            for k = 1:M
                #ed = (g*(vy - 1) + m)*M + k + Ns_hor;
                x0 = lx*N*(g - 1);
                y0 = ly*(M*(m - 1) + k - 1);
                start_point = [x0, y0, grid.depth];
                end_point = [x0, y0 + ly, grid.depth];
                electrodes[ed] = new_electrode(start_point, end_point, grid.radius);
                ed += 1;
                if (k < M)
                    nodes[1, nd] = end_point[1];
                    nodes[2, nd] = end_point[2];
                    nodes[3, nd] = end_point[3];
                    nd += 1;
                end
            end
        end
    end
    return electrodes, transpose(nodes)
end

function plot_xy(electrodes, nodes=false)
	""" Plots the electrodes in the xy plane. """
    #fig = plot(leg=false, grid=false)#, aspect_ratio=1)# border=:none)
	#fig = plot(leg=false, grid=false, aspect_ratio=1)# border=:none)
    fig = plot(grid=false, aspect_ratio=1)# border=:none)
    for e in electrodes
        x0 = e.start_point[1]
		y0 = e.start_point[2]
        x1 = e.end_point[1]
		y1 = e.end_point[2]
        plot!([x0, x1], [y0, y1], line=(1, :black, :solid), label="")
        nodes && scatter!([x0, x1], [y0, y1], markercolor=:black, markersize=2, label="")
    end
    return fig
end

begin
    r = 5e-3
    h = -0.5
    lx = 20
    ly = 16
    grid = Grid(6, 5, lx, ly, 1, 1, r, h)
    electrodes, nodes = electrode_grid(grid)
    num_electrodes = ns = length(electrodes)
    num_nodes = nn = size(nodes)[1]
    plot_xy(electrodes, true)
    h = sqrt(4^2 + 5^2)
    cosθ = 5 / h
    sinθ = 4 / h
    x = -6:0.1:26
    y = x * 4/5
    plot!(x, y, label="ξ axis", xlabel="x [m]", ylabel="y [m]", color=:red,
          xlims=(-1, 29), ylims=(-1, 17))
    png("visacro57emc01_grid")
end

## GPR and injection
begin
    t = range(0, 30e-6, step=20e-9)
    i1 = @. heidler(t, 1.09610481e+00, 5.47446858e-07, 1.91552576e-06, 2.94082573e+00)
    i2 = @. heidler(t, 5.06486595e-01, 3.58100804e-06, 5.40426910e-05, 3.41650905e+00)
    inj_fast = i1 + i2
    i1 = @. heidler(t, 5.31736342e-01, 5.66828155e-06, 7.40316360e-07, 6.45026985e+00)
    i2 = @. heidler(t, 5.86866680e-01, 4.65351493e-06, 5.40167336e-05, 1.07878428e+01)
    inj_slow =  i1 + i2
    plot(t * 1e6, [inj_fast, inj_slow], label=["fast wave"  "slow wave"],
         xlabel="Time [μs]", ylabel="Current [A]")
    png("visacro57emc01_inj")
end

begin
    m = 1e6
    df = DataFrame(CSV.File("visacro57emc01_gpr.csv"))
    t = df.t * m
    nt = length(t)
    t2 = range(0, length=nt, step=20e-9) * m
    nend = nt - nt ÷ 20
    #1  ======================================================================
    p = plot(t[1:nend], df.i1[1:nend], label="corner, fast", linecolor=:blue,
             xlabel="time [μs]", ylabel="GPR [V]", ylims=(-1, 25), legend=:bottomright)
    plot!(t[1:nend], df.i2[1:nend], label="center, fast", linecolor=:orange)
    plot!(t[1:nend], df.i3[1:nend], label="corner, slow", linecolor=:green)
    plot!(t[1:nend], df.i4[1:nend], label="center, slow", linecolor=:black)
    png(p, "visacro57emc01_gpr")
end

## Electric Field
time_steps = [6, 26, 46, 76, 101, 126, 251];
begin
    df = DataFrame(CSV.File("visacro57emc01_efield.csv"))
    t = sort(unique(df.t))
    nt = length(t)
    ninj = (size(df)[2] - 3) ÷ 6  # number of injections
    for i = 1:ninj
        # do some Metaprogramming wizardry
        ecx = eval(Meta.parse("df.ecx$(i)"))
        ecy = eval(Meta.parse("df.ecy$(i)"))
        enoncx = eval(Meta.parse("df.encx$(i)"))
        enoncy = eval(Meta.parse("df.ency$(i)"))
        for k in time_steps
            c = isapprox.(df.t, t[k])
            x = df.x[c]
            y = df.y[c]
            Lx = x[end] - x[1]
            Ly = y[end] - y[1]
            Lr = sqrt(Lx^2 + Ly^2)
            cosθ = Lx / Lr
            sinθ = Ly / Lr
            ξ = x * cosθ + y * sinθ
            efc = ecx[c] * cosθ + ecy[c] * sinθ
            efnonc = enoncx[c] * cosθ + enoncy[c] * sinθ
            m = 1
            p = plot(ξ, efc * m, label="conservative", color=:blue,
                     ylabel="E [V/m]", xlabel="ξ [m]", legend=:bottomright)
                     #title="t = $(round(1.0e6 * t[k], sigdigits=2)) [μs]");
            plot!(ξ, efnonc * m, label="non-conservative", color=:red);
            plot!(ξ, (efc + efnonc) * m, label="total", color=:darkslategray)
            png(p, "visacro57emc01_efield$(i)_t$(k)")
        end
    end
end

## Step Voltage along a line in the +ξ direction
begin
    df = DataFrame(CSV.File("visacro57emc01_efield.csv"))
    t = sort(unique(df.t))
    nt = length(t)
    ninj = (size(df)[2] - 3) ÷ 6  # number of injections
    for i = 1:ninj
        ecx = eval(Meta.parse("df.ecx$(i)"))
        ecy = eval(Meta.parse("df.ecy$(i)"))
        enoncx = eval(Meta.parse("df.encx$(i)"))
        enoncy = eval(Meta.parse("df.ency$(i)"))
        for k in time_steps
            p = plot(ylabel="Step Voltage [V]", xlabel="ξ [m]", legend=:bottomright);
            m = 1
            c = isapprox.(df.t, t[k])
            x = df.x[c]
            y = df.y[c]
            Lx = x[end] - x[1]
            Ly = y[end] - y[1]
            Lr = sqrt(Lx^2 + Ly^2)
            cosθ = Lx / Lr
            sinθ = Ly / Lr
            ξ = x * cosθ + y * sinθ
            efc = ecx[c] * cosθ + ecy[c] * sinθ
            efnonc = enoncx[c] * cosθ + enoncy[c] * sinθ
            eftot = efc + efnonc
            dξ = ξ[2] - ξ[1]
            N = Int(1 ÷ dξ)
            nξ = length(ξ) - N
            dv1 = zeros(nξ)
            dv2 = zeros(nξ)
            for i = 1:nξ
                s = sum(efc[(i + 1):(i + N - 1)])
                dv1[i] = dξ * (s + (efc[i] + efc[i + N]) / 2)
                s = sum(eftot[(i + 1):(i + N - 1)])
                dv2[i] = dξ * (s + (eftot[i] + eftot[i + N]) / 2)
            end
            ξ = ξ[1:end-N]
            plot!(ξ, (dv1) * m,  label="pot. dif.", color=:blue)
            plot!(ξ, (dv2) * m, label="∫E⋅dℓ", color=:red)
            png(p, "visacro57emc01_stepv$(i)_t$(k)")
        end
    end
end

## RMSD(ξ)
@time begin
    df = DataFrame(CSV.File("visacro57emc01_efield.csv"))
    t = sort(unique(df.t))
    nt = length(t)
    nf = 401
    #nf = nt ÷ 2 + 1
    ninj = (size(df)[2] - 3) ÷ 6  # number of injections
    c = isapprox.(df.t, t[1])
    x = df.x[c]
    y = df.y[c]
    Lx = x[end] - x[1]
    Ly = y[end] - y[1]
    Lr = sqrt(Lx^2 + Ly^2)
    cosθ = Lx / Lr
    sinθ = Ly / Lr
    ξ = x * cosθ + y * sinθ
    dξ = ξ[2] - ξ[1]
    N = Int(1 ÷ dξ)
    nξ = length(ξ) - N
    ξ = ξ[1:end-N]
    dv1 = zeros(nξ)
    dv2 = zeros(nξ)
    rmsd = zeros(nf, ninj)
    for i = 1:ninj
        ecx = eval(Meta.parse("df.ecx$(i)"))
        ecy = eval(Meta.parse("df.ecy$(i)"))
        enoncx = eval(Meta.parse("df.encx$(i)"))
        enoncy = eval(Meta.parse("df.ency$(i)"))
        for k = 1:nf
            c = isapprox.(df.t, t[k])
            efc = ecx[c] * cosθ + ecy[c] * sinθ
            efnonc = enoncx[c] * cosθ + enoncy[c] * sinθ
            eftot = efc + efnonc
            for i = 1:nξ
                s = sum(efc[(i + 1):(i + N - 1)])
                dv1[i] = dξ * (s + (efc[i] + efc[i + N]) / 2)
                s = sum(eftot[(i + 1):(i + N - 1)])
                dv2[i] = dξ * (s + (eftot[i] + eftot[i + N]) / 2)
            end
            rmsd[k, i] = sqrt(sum((dv1 - dv2) .^ 2 / nξ))
        end
    end
    p = plot(t[1:nf] * 1e6, rmsd * 1,
             xlabel="time [μs]", ylabel="Step Voltage RMSD(ξ) [V]",
             label=["corner, fast"  "center, fast"  "corner, slow"  "center, slow"],
             color=[:blue  :red  :blue  :red],
             linestyle=[:solid  :solid  :dash  :dash])
    png(p, "visacro57emc01_rmsd.png")
end

## GPD profile
#=
import PlotlyJS
begin
    df = DataFrame(CSV.File("visacro57emc01_gpd.csv"))
    t = sort(unique(df.t))
    nt = length(t)
    ninj = (size(df)[2] - 3)  # number of injections
    for i = 1:ninj
        gpd = eval(Meta.parse("df.i$(i)"))
        for k in time_steps
            c = isapprox.(df.t, t[k])
            data = PlotlyJS.contour(; x=df.x[c], y=df.y[c], z=gpd[c],
                                    colorbar=PlotlyJS.attr(;title="Electric Potential [V]",
                                                  titleside="right"),
                                                  titlefont=PlotlyJS.attr(;size=16,
                                                  family="Arial, sans-serif"))
            p = PlotlyJS.plot(data)
            PlotlyJS.savefig(p, "visacro57emc01_gpd$(i)_t$(k).png")
        end
    end
end
=#
