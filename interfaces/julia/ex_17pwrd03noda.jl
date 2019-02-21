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

function incidence(electrodes::Vector{Electrode}, nodes::Matrix{Float64})
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
end

function laplace_transform(f::Vector{ComplexF64}, t::Vector{Float64}, s::Vector{ComplexF64})
    nt = length(t);
    nf = length(s);
    res = zeros(Complex{Float64}, nf);
    for k = 1:nf
        for i = 1:(nt-1)
            e0 = exp(s[k]*t[i]);
            e1 = exp(s[k]*t[i+1]);
            dt = t[i+1] - t[i];
            x = (f[i+1] - f[i])/s[k];
            res[k] += (e1*(f[i+1]*dt - x) - e0*(f[i]*dt - x))/dt;
        end
        res[k] = 1/s[k]*res[k];
    end
    return res
end;

# =======================================================================
## Parameters
mu0 = pi*4e-7;
mur = 1;
eps0 = 8.854e-12;
epsr = 1;
sigma1 = 0;
#rhoc = 1.9 10^-6;
rhoc = 1.68e-8;
sigmac = 1/rhoc;
rho_lead = 2.20e-7;
sigma_lead = 1/rho_lead;

rsource = 50.0;
gf = 1.0/rsource;
rh = 15e-3;
rv = 10e-3;
h = 0.5;
l = 4.0;

#= Frequencies
Due to numerical erros, to smooth the response, its necessary to use a
final time much greater than that up to which is desired.
=#
T = 0.7e-7*2;
dt = 2.0e-9;
n = T/dt;
t = collect(0.0:dt:(T-dt));
sc = log(n^2)/T;
kk = collect(0:1:n/2);
dw = 2.0*pi/(n*dt);
function sigma(j, alpha=0.53836)
    return alpha + (1 - alpha)*cos(2*pi*j/n)
end;
sk = -1im*sc*ones(length(kk)) + dw*kk;
nf = length(sk);
freq = real(sk)/(2*pi);
omega = 2*pi*freq[nf];
lambda = (2*pi/omega)*(1/sqrt( epsr*eps0*mu0/2*(1 + sqrt(1 + (sigma1/(omega*epsr*eps0))^2)) ));

# Integration Parameters
max_eval = 200;
req_abs_error = 1e-4;
req_rel_error = 1e-5;
error_norm = 1; #paired
intg_type = 1; #double integral

## Electrodes
x = 10;
nv = Int(ceil(h/(lambda/x)));
nh = Int(ceil(l/(lambda/x)));
vertical = new_electrode([0, 0, 0], [0, 0, 0.5], 10e-3, 0.0);
horizontal = new_electrode([0, 0, 0.5], [4.0, 0, 0.5], 15e-3, 0.0);
elecv, nodesv = segment_electrode(vertical, Int(nv));
elech, nodesh = segment_electrode(horizontal, Int(nh));

electrodes = elecv;
append!(electrodes, elech);
ns = length(electrodes);
nodes = cat(nodesv, nodesh, dims=1);
nn = size(nodes)[1];

mA, mB = incidence(electrodes, nodes);
mAT = transpose(mA);
mBT = transpose(mB);
yn = zeros(Complex{Float64}, (nn,nn));
zl = zeros(Complex{Float64}, (ns,ns));
zt = zeros(Complex{Float64}, (ns,ns));
exci = zeros(Complex{Float64}, (nn));
vout = zeros(Complex{Float64}, (nf,nn));

## Source input
using DataFrames, CSV
path = "/home/pedro/codigos/HP_HEM/interfaces/julia/";
input = string(path, "source.csv");
source = CSV.read(input, header=["t", "V"]);
source[:,1] = source[:,1]*1e-9;

input = string(path, "voltageArt.csv");
voltageArt = CSV.read(input, header=["t", "V"]);
input = string(path, "currentArt.csv");
currentArt = CSV.read(input, header=["t", "I"]);

ent_freq = laplace_transform(Vector{ComplexF64}(source.V), Vector{Float64}(source.t), -1.0im*sk);

## Freq. loop
println("starting loop")
for i = 1:nf
    kappa = 1im*sk[i]*eps0;
    k1 = sqrt(1im*sk[i]*mu0*kappa);
    zl, zt = calculate_impedances(electrodes, k1, sk[i], mur, kappa,
                                  max_eval, req_abs_error, req_rel_error,
                                  error_norm, intg_type);
    yn = mAT*inv(zt)*mA + mBT*inv(zl)*mB;
    exci[1] = ent_freq[i]*gf;
    vout[i,:] = yn\exci;
end

## Time response
tf = 40;
#vout = rand(nf,nn)+rand(nf,nn)*1.0im
outlow = map(i -> vout[:,1][Int(i+1)]*sigma(i), kk);
upperhalf = reverse(conj(outlow));
pop!(upperhalf);
lowerhalf = outlow;
pop!(lowerhalf);
append!(lowerhalf, upperhalf);
F = lowerhalf;

using FFTW
f = real(ifft(F));
out = map(i -> exp(sc*t[i])/dt*f[i], 1:length(t));
