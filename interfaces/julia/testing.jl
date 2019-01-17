#=
Test: using the Julia interface
=#
struct Electrode
    start_point::Array{Float64, 1}
    end_point::Array{Float64, 1}
    middle_point::Array{Float64, 1}
    length::Float64
    radius::Float64
    internal_impedance::Complex128
end;


function new_electrode(start_point, end_point, radius, internal_impedance)
    return Electrode(start_point, end_point, (start_point + end_point)/2, norm(start_point - end_point), radius, internal_impedance)
end;

w = 2*pi*60;
e0 = 8.854187817620e-12;
sigma = 1.0/600;
epsr = 15.0;
gamma = sqrt(sigma + 1im*w*epsr*e0);
radius = 1e-2;

sender = new_electrode([0, 0, 0], [1, 0, 0], 0.1, 0.1);
receiver = new_electrode([0, 1, 0], [1, 1, 0], 0.1, 0.1);

maxEval = 0;
reqAbsError = 0.0;
reqRelError = 1e-4;

result = [0.1, 0.2];
error = [0.1, 0.2];

# path e nome devem setar direto no call
ccall(("integral_log_Nf", "/home/pedro/codigos/HEM_C/hem.so"), Int, (Ref{Electrode}, Ref{Electrode}, Complex128, UInt, Float64, Float64, Ref{Float64}, Ref{Float64}), sender, receiver, gamma, maxEval, reqAbsError, reqRelError, result, error);

#ccall(("integral_log_Nf", "/home/pedro/codigos/HEM_C/hem.so"), Int, (Electrode, Electrode, Complex128, UInt, Float64, Float64, Ref{Float64}, Ref{Float64}), sender, receiver, gamma, maxEval, reqAbsError, reqRelError, result, error);

val = [0.,0.];
err = [1.,1.];
ccall(("integrate", "/home/pedro/codigos/HEM_C/hem.so"), Int, (Ref{Float64}, Ref{Float64}), val, err);

#############
