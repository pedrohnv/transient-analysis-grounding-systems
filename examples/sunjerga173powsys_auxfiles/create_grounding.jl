#= creates and plots the wind grounding from:
SUNJERGA, Antonio et al. Isolated vs. interconnected wind turbine grounding
systems: Effect on the harmonic grounding impedance, ground potential rise
and step voltage. Electric Power Systems Research, v. 173, p. 230-239, 2019.
=#
using LinearAlgebra
using Plots

mutable struct Electrode
    """ Defines a conductor segment. """
    start_point::Array{Float64,1}
    end_point::Array{Float64,1}
    middle_point::Array{Float64,1}
    length::Float64
    radius::Float64
end

function new_electrode(start_point, end_point, radius)
    """ Creates a conductor segment. """
    return Electrode(start_point, end_point, (start_point + end_point)/2.0,
                     norm(start_point - end_point), radius)
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

function segment_electrode(electrode::Electrode, num_segments::Int)
    """ Segments a conductor. """
    nn = num_segments + 1;
    nodes = Array{Float64}(undef, 3, nn);
    startp = Array{Float64,1}(undef, 3);
    endp = Array{Float64,1}(undef, 3);
    for k = 1:3
        startp[k] = electrode.start_point[k];
        endp[k] = electrode.end_point[k];
    end
    increment = (endp - startp)/num_segments;
    for k = 0:num_segments
        nodes[:,k+1] = startp + k*increment;
    end
    segments = Array{Electrode,1}(undef, num_segments);
    for k = 1:num_segments
        segments[k] = new_electrode(nodes[:,k], nodes[:,k+1], electrode.radius);
    end
    return segments, nodes
end

function seg_electrode_list(electrodes, lmax)
    """
    Segments a list of conductors such that they end up having at most 'lmax'
    length.
    Return a list of the segmented conductors and their nodes.
    """
    num_elec = 0; #after segmentation
    for i=1:length(electrodes)
        #TODO store in array to avoid repeated calculations
        num_elec += Int(ceil(electrodes[i].length/lmax));
    end
    elecs = Array{Electrode}(undef, num_elec);
    nodes = zeros(Float64, (3, 2*num_elec));
    e = 1;
    nodes = [];
    for i = 1:length(electrodes)
        ns = Int(ceil(electrodes[i].length/lmax));
        new_elecs, new_nodes = segment_electrode(electrodes[i], ns);
        for k=1:ns
            elecs[e] = new_elecs[k];
            e += 1;
        end
        if (nodes == [])
            nodes = new_nodes;
        else
            for k = 1:size(new_nodes)[2]
                if (matchcol(new_nodes[:,k:k], nodes) == nothing)
                    nodes = cat(nodes, new_nodes[:,k:k], dims=2);
                end
            end
        end
    end
    return elecs, nodes
end

function incidence(electrodes, nodes; atol=0, rtol=1e-4)
    """
    Builds incidence matrices A and B for calculating the nodal admittance
    matrix: YN = AT*inv(ZL)*A + BT*inv(ZT)*B
    """
    ns = length(electrodes);
    nn = size(nodes)[2];
    a = zeros(ComplexF64, (ns,nn));
    b = zeros(ComplexF64, (ns,nn));
    for k = 1:nn
        for i = 1:ns
            if isapprox(collect(electrodes[i].start_point), nodes[:,k],
                        atol=atol, rtol=rtol)
                a[i,k] = 1.0;
                b[i,k] = 0.5;
            elseif isapprox(collect(electrodes[i].end_point), nodes[:,k],
                            atol=atol, rtol=rtol)
                a[i,k] = -1.0;
                b[i,k] = 0.5;
            end
        end
    end
    return a, b
end

function electrode_ring(r, segments::Int, z=0.0, radius=1e-3)
	""" Generates an Electrode ring circumscribed by a radius r circle. """
    dt = 2*pi/segments;
    angles = 0:dt:(2*pi - dt);
    nodes = [[r*cos(t), r*sin(t), z] for t in angles]
    s1 = segments - 1;
    #electrodes = [[nodes[i], nodes[i + 1], radius] for i in 1:1:s1];
	electrodes = [new_electrode(nodes[i], nodes[i + 1], radius) for i in 1:1:s1];
    push!(electrodes, new_electrode(nodes[segments], nodes[1], radius));
    return electrodes, nodes
end

function conect_rings(nodes1, nodes2, jump::Int=0, radius::Float64=1e-3)
	""" Conects nodes1[i] to nodes2[i] skipping every i = (jump+1). """
    n1 = length(nodes1);
    n2 = length(nodes2);
	n = (n1 < n2) ? n1 : n2;
	notjump(i) = !Bool((i - 1)%(jump + 1));
    #return [[nodes1[i], nodes2[i], radius] for i in 1:1:n if notjump(i)]
	return [new_electrode(nodes1[i], nodes2[i], radius) for i in 1:1:n if notjump(i)]
end

function plot_elecnodes(electrodes, nodes)
	"""
	Makes a 2D, xy-coordinates, plot of the electrodes and nodes.
	Call 'using Plots' and the backend of your choice before calling this function.
	"""
    num_electrodes = length(electrodes);
	points = nodes[1:2,:];
	p = scatter(points[1,:], points[2,:], legend=false, markercolor=:black, border=:none)
	for i=1:num_electrodes-1
	    e = electrodes[i];
	    plot!([e.start_point[1], e.end_point[1]], [e.start_point[2], e.end_point[2]],
	          line=(:black))
	end
	e = electrodes[end];
	plot!([e.start_point[1], e.end_point[1]], [e.start_point[2], e.end_point[2]],
	      line=(:black), legend=false, aspect_ratio=1)
	return p
end

"""
Makes a 3D plot of the electrodes and nodes.
Call 'using Plots' and the backend of your choice before calling this function.
"""
function plot_elecnodes3d(electrodes, nodes, camera=(10,45), msize=1)
    num_electrodes = length(electrodes);
	p = scatter(nodes[1,:], nodes[2,:], nodes[3,:], legend=false, markercolor=:black,
                camera=camera, msize=msize, fmt=:png)
    for i=1:num_electrodes
        e = electrodes[i];
        x = [e.start_point[1], e.end_point[1]];
        y = [e.start_point[2], e.end_point[2]];
        z = [e.start_point[3], e.end_point[3]];
        plot!(x, y, z, line=(:black))
    end
    e = electrodes[end];
    x = [e.start_point[1], e.end_point[1]];
    y = [e.start_point[2], e.end_point[2]];
    z = [e.start_point[3], e.end_point[3]];
    plot!(x, y, z, line=(:black), legend=false)#, aspect_ratio=1)
	return p
end

function sunjerga173powsys_seg(lmax1=1, lmax2=1)
	"""
	Creates wind grounding from:

	SUNJERGA, Antonio et al. Isolated vs. interconnected wind turbine grounding
	systems: Effect on the harmonic grounding impedance, ground potential rise
	and step voltage. Electric Power Systems Research, v. 173, p. 230-239, 2019.

	lmax1 : first ring's conductors max length
	lmax2 : other conductors max length
	"""
    radius = 1e-2;
    edges = 8;
    electrodes = Array{Electrode,1}();
    # Rings
	a = 1
	h = -0.05
	# first ring and cross

	elecs1, nodes1 = electrode_ring(2.6*a, edges, h, radius);
	elecs2, nodes2 = electrode_ring(2.6*a, edges, -1, radius);
    elecs3, nodes3 = electrode_ring(5.8*a, edges, -1.5, radius);
    elecs4, nodes4 = electrode_ring(9*a, edges, -2, radius);
    elecs5, nodes5 = electrode_ring(9*a, edges, -3, radius);
	# don't append elecs1
	append!(electrodes, elecs2);
    append!(electrodes, elecs3);
    append!(electrodes, elecs4);
    append!(electrodes, elecs5);
	# ring-connection conductors
	append!(electrodes, conect_rings(nodes1, nodes2, 1, radius));
	append!(electrodes, conect_rings(nodes2, nodes3, 1, radius));
	append!(electrodes, conect_rings(nodes3, nodes4, 1, radius));
	append!(electrodes, conect_rings(nodes4, nodes5, 1, radius));
	# vertical rods
    elecs6, nodes6 = electrode_ring(9*a, edges, -7, radius);
    rods = conect_rings(nodes5, nodes6, 1, radius);
	append!(electrodes, rods);
	# segment first ring and cross
	#append!(elecs1, [new_electrode([0.,0.,h], nodes1[i], radius) for i=1:2:length(nodes1)])
	#elecs1, nodes1 = seg_electrode_list(elecs1, lmax1)
	cross = [new_electrode([0.,0.,h], nodes1[i], radius) for i=1:2:length(nodes1)]
	cross, dummy = seg_electrode_list(cross, lmax1)
	append!(electrodes, elecs1);
    # segment other conductors
	electrodes, nodes = seg_electrode_list(electrodes, lmax2)
	#append!(electrodes, elecs1);
	append!(electrodes, cross);
	# "segment" to get list of unique nodes
	electrodes, nodes = seg_electrode_list(electrodes, 100)
	ne = num_electrodes = length(electrodes)
	nn = num_nodes = length(nodes)
    return electrodes, nodes
end

electrodes, nodes = sunjerga173powsys_seg(0.1, 0.5)
plot_elecnodes3d(electrodes, nodes, (10, 35))
nome = "sunjerga173powsys"
begin
	f = open("electrodes.csv", "w")
	for e in electrodes
		for i = 1:3
			write(f, "$(e.start_point[i]), ")
		end
		for i = 1:3
			write(f, "$(e.end_point[i]), ")
		end
		write(f, "$(e.radius)\n")
	end
	close(f)
end
begin
	f = open("nodes.csv", "w")
	for i = 1:length(nodes)
		write(f, "$(nodes[i])")
		if i % 3 == 0
			write(f, "\n")
		else
			write(f, ",")
		end
	end
	close(f)
end
