#=
Julia interface to the High Performance Hybrid Electromagnetic Model program.
=#
using LinearAlgebra;

struct Electrode
	#had to use tuple to successfully pass to C as a pointer
    start_point::NTuple{3,Cdouble}
    end_point::NTuple{3,Cdouble}
    middle_point::NTuple{3,Cdouble}
    length::Float64
    radius::Float64
    zi::Complex{Float64}
end;

@enum Integration_type begin
	INTG_NONE = 1
	INTG_DOUBLE = 2
	INTG_EXP_LOGNF = 3
	INTG_LOGNF = 4
end

function new_electrode(start_point, end_point, radius, internal_impedance)
    return Electrode(NTuple{3,Cdouble}(start_point), NTuple{3,Cdouble}(end_point),
                     NTuple{3,Cdouble}((start_point + end_point)/2.0),
                     norm(start_point - end_point), radius, internal_impedance)
end;

"""
Makes a 2D, xy-coordinates, plot of the electrodes and nodes.
Call 'using Plots' before calling this function.
"""
function plot_elecnodes(electrodes, nodes)
	pontos = nodes[:,1:2]
	scatter(pontos[:,1], pontos[:,2], legend=false)
	for i=1:num_electrodes-1
		e = electrodes[i];
		plot!([e.start_point[1], e.end_point[1]], [e.start_point[2], e.end_point[2]],
		      line=(:black))
	end
	e = electrodes[end];
	plot!([e.start_point[1], e.end_point[1]], [e.start_point[2], e.end_point[2]],
	      line=(:black), legend=false, aspect_ratio=1)
end;

function segment_electrode(electrode::Electrode, num_segments::Int)
    nn = num_segments + 1;
    nodes = Array{Float64}(undef, nn, 3);
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
        segments[k] = new_electrode(nodes[k,:], nodes[k+1,:], electrode.radius,
                                    electrode.zi);
    end
    return segments, nodes
end;

function matchrow(a, B, atol=1e-9, rtol=0)
	#= Returns the row in B that matches a. If there is no match, nothing is returned.
	taken and modified from from https://stackoverflow.com/a/32740306/6152534 =#
	findfirst(i -> all(j -> isapprox(a[j], B[i,j], atol=atol, rtol=rtol), 1:size(B,2)), 1:size(B,1))
end

function seg_electrode_list(electrodes, frac)
	num_elec = 0; #after segmentation
	for i=1:length(electrodes)
		#TODO store in array to avoid repeated calculations
		num_elec += Int(ceil(electrodes[i].length/frac));
	end
	elecs = Array{Electrode}(undef, num_elec);
	#nodes = zeros(Float64, (2*num_elec, 3));
	e = 1;
	nodes = [];
	for i=1:length(electrodes)
		ns = Int(ceil(electrodes[i].length/frac));
		new_elecs, new_nodes = segment_electrode(electrodes[i], ns);
		for k=1:ns
			elecs[e] = new_elecs[k];
			e += 1;
		end
		if (nodes == [])
			nodes = new_nodes;
		else
			for k=1:size(new_nodes)[1]
				if (matchrow(new_nodes[k:k,:], nodes) == nothing)
					nodes = cat(nodes, new_nodes[k:k,:], dims=1);
				end
			end
		end
	end
	return elecs, nodes
end;

function electrode_grid(a, n::Int, b, m::Int, h, r, zi=0.0im)
	#=
	Creates an electrode grid `h` coordinate below ground with each conductor
	having radius `r` and internal impedance `zi`.
	The grid has dimensions `a*b` with `n` and `m` divisions respectively.
	=#
	xx = 0:a/n:a;
	yy = 0:b/m:b;
	num_elec = n*(m + 1) + m*(n + 1);
	electrodes = Array{Electrode}(undef, num_elec);
	e = 1;
	for k=1:(m+1)
		for i=1:n
			electrodes[e] = new_electrode([xx[i], yy[k], h], [xx[i+1], yy[k], h], r, zi);
			e += 1;
		end
	end
	for k=1:(n+1)
		for i=1:m
			electrodes[e] = new_electrode([xx[k], yy[i], h], [xx[k], yy[i+1], h], r, zi);
			e += 1;
		end
	end
	# TODO return nodes as well?
	return electrodes
end;

function calculate_impedances(electrodes, gamma, s, mur, kappa, max_eval,
                              req_abs_error, req_rel_error, error_norm, intg_type)
    ns = length(electrodes);
    zl = zeros(Complex{Float64}, (ns,ns));
    zt = zeros(Complex{Float64}, (ns,ns));
    # path to the library must be a static symbol?
    # see https://stackoverflow.com/questions/35831775/issue-with-julia-ccall-interface-and-symbols
    ccall(("calculate_impedances", "C:\\Users\\pedro\\Documents\\codigos\\HP_HEM\\libhem_electrode.dll"), Int,
          (Ref{Electrode}, Int, Ref{Complex{Float64}}, Ref{Complex{Float64}},
          Complex{Float64}, Complex{Float64}, Float64, Complex{Float64},
          Int, Float64, Float64, Int, Int),
          electrodes, ns, zl, zt, gamma, s, mur, kappa,
          max_eval, req_abs_error, req_rel_error, error_norm, intg_type);
    return zl, zt
end;

function impedances_images(electrodes, images, zl, zt, gamma, s, mur, kappa,
						   ref_l, ref_t, max_eval, req_abs_error,
						   req_rel_error, error_norm, intg_type)
    ns = length(electrodes);
    # path to the library must be a static symbol?
    # see https://stackoverflow.com/questions/35831775/issue-with-julia-ccall-interface-and-symbols
    ccall(("impedances_images", "C:\\Users\\pedro\\Documents\\codigos\\HP_HEM\\libhem_electrode.dll"), Int,
          (Ref{Electrode}, Ref{Electrode}, Int, Ref{Complex{Float64}},
		  Ref{Complex{Float64}}, Complex{Float64}, Complex{Float64}, Float64,
		  Complex{Float64}, Complex{Float64}, Complex{Float64}, Int, Float64,
		  Float64, Int, Int),
          electrodes, images, ns, zl, zt, gamma, s, mur, kappa, ref_l, ref_t,
          max_eval, req_abs_error, req_rel_error, error_norm, intg_type);
    return zl, zt
end;

function incidence(electrodes::Vector{Electrode}, nodes::Matrix{Float64})
	# build incidence matrices for calculating 'YN = AT*inv(zt)*A + BT*inv(zl)*B'
    ns = length(electrodes);
    nn = size(nodes)[1];
    a = zeros(Float64, (ns,nn));
    b = zeros(Float64, (ns,nn));
    for i = 1:ns
        for k = 1:nn
            if isapprox(collect(electrodes[i].start_point), nodes[k,:])
                a[i,k] = 0.5;
                b[i,k] = 1.0;
            elseif isapprox(collect(electrodes[i].end_point), nodes[k,:])
                a[i,k] = 0.5;
                b[i,k] = -1.0;
            end
        end
    end
    return a, b
end;

function immittance(electrodes, nodes, zl, zt, ye);
	# build global immittance matrix
	m = length(electrodes);
    n = size(nodes)[1];
	N = 2*m + n;
    we = zeros(ComplexF64, (N,N));
    for k = 1:n
        for i = 1:m
            if isapprox(collect(electrodes[i].start_point), nodes[k,:])
                we[i+n, k] = -1.0;
                we[i+n+m, k] = -0.5;
				we[k, i+n] = 1.0;
            elseif isapprox(collect(electrodes[i].end_point), nodes[k,:])
				we[i+n, k] = 1.0;
                we[i+n+m, k] = -0.5;
                we[k, i+n+m] = 1.0;
            end
        end
    end
	for k = 1:n
		for i = 1:n
			we[i, k] = ye[i, k];
		end
	end
	for k = 1:m
		for i = 1:m
			zl2 = zl[i, k]/2;
			we[i+n, k+n] = zl2;
			we[i+n, k+n+m] = -zl2;
			we[i+n+m, k+n] = zt[i, k];
			we[i+n+m, k+n+m] = zt[i, k];
		end
	end
	return we
end;
