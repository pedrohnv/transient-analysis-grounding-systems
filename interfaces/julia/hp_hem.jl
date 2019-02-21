#=
Test: using the Julia interface
=#
using LinearAlgebra;

struct Electrode
    start_point::NTuple{3,Cdouble}
    end_point::NTuple{3,Cdouble}
    middle_point::NTuple{3,Cdouble}
    length::Float64
    radius::Float64
    zi::Complex{Float64}
end;

function new_electrode(start_point, end_point, radius, internal_impedance)
    return Electrode(NTuple{3,Cdouble}(start_point), NTuple{3,Cdouble}(end_point),
                     NTuple{3,Cdouble}((start_point + end_point)/2.0),
                     norm(start_point - end_point), radius, internal_impedance)
end;

function segment_electrode(electrode::Electrode, num_segments::Int)
    nn = num_segments + 1;
    nodes = Array{Float64}(undef, nn, 3)
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
        segments[k] = new_electrode(nodes[k,:], nodes[k+1,:], electrode.radius, electrode.zi);
    end
    return segments, nodes
end;

function calculate_impedances(electrodes, gamma, s, mur, kappa,
                              max_eval, req_abs_error, req_rel_error,
                              error_norm, intg_type)
    ns = length(electrodes);
    zl = zeros(Complex{Float64}, (ns,ns));
    zt = zeros(Complex{Float64}, (ns,ns));
    # path to the library must be a static symbol?
    # see https://stackoverflow.com/questions/35831775/issue-with-julia-ccall-interface-and-symbols
    ccall(("calculate_impedances", "/home/pedro/codigos/HP_HEM/libhem"), Int,
          (Ref{Electrode}, Int, Ref{Complex{Float64}}, Ref{Complex{Float64}},
          Complex{Float64}, Complex{Float64}, Float64, Complex{Float64},
          Int, Float64, Float64, Int, Int),
          electrodes, ns, zl, zt, gamma, s, mur, kappa,
          max_eval, req_abs_error, req_rel_error, error_norm, intg_type);
    return zl, zt
end;

r1 = 7e-3;
zi = 0.0;
electrodes = [
    new_electrode([0, 0, 0], [1, 0, 0], r1, zi),
    new_electrode([0, 0, 0], [0, 1, 0], r1, zi),
    new_electrode([0, 0, 0], [0, 0, 1], r1, zi)
];

zl = zeros(Complex{Float64}, (3,3));
zt = zeros(Complex{Float64}, (3,3));

ccall(("calculate_impedances", "/home/pedro/codigos/HP_HEM/libhem"), Int,
      (Ref{Electrode}, Int, Ref{Complex{Float64}}, Ref{Complex{Float64}},
      Complex{Float64}, Complex{Float64}, Float64, Complex{Float64},
      Int, Float64, Float64, Int, Int),
      electrodes, 3, zl, zt, 1.0, 1.0, 1.0, 1.0, 200, 1e-3, 1e-4, 2, 1);

println("zl =");
for i = 1:3
    for k = 1:3
        print(zl[i,k], "   ");
    end
    println();
end

println("zt =");
for i = 1:3
    for k = 1:3
        print(zt[i,k], "   ");
    end
    println();
end
#############
#zl =
#   9.3240e-07   0.0000e+00   0.0000e+00
#   0.0000e+00   9.3240e-07   0.0000e+00
#   0.0000e+00   0.0000e+00   9.3240e-07
#
#zt =
#   0.741977   0.085033   0.085033
#   0.085033   0.741977   0.085033
#   0.085033   0.085033   0.741977

#############
w = 2*pi*60;
e0 = 8.854187817620e-12;
mu0 = 1.256637061435917e-6;
mur = 1.0;
sigma = 1.0/600;
epsr = 15.0;
kappa = sigma + 1im*w*epsr*e0;
gamma = sqrt(1im*w*mur*mu0*kappa);
