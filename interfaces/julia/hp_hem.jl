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
    internal_impedance::Complex{Float64}
end;

function new_electrode(start_point, end_point, radius, internal_impedance)
    return Electrode(NTuple{3,Cfloat}(start_point), NTuple{3,Cfloat}(end_point),
                     NTuple{3,Cfloat}((start_point + end_point)/2.0),
                     norm(start_point - end_point), radius, internal_impedance)
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
