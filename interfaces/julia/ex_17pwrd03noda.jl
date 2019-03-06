#=
Test case 17pwrd03noda

Reproducing the results in [1] of a time domain surge response.

[1] Noda, Taku, and Shigeru Yokoyama. "Thin wire representation in finite difference time domain
surge simulation." IEEE Transactions on Power Delivery 17.3 (2002): 840-847.
=#
using LinearAlgebra
using DataFrames, CSV
using FFTW
using Plots

include("hp_hem.jl");

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

#=
Put in a function as execution in global scope has poor performance.
The computation becomes really fast, but the pre-compiling Julia does
takes a lot of time...
=#
function simulate()
	## Parameters
	mu0 = pi*4e-7;
	mur = 1;
	eps0 = 8.854e-12;
	epsr = 1;
	sigma1 = 0;
	#rhoc = 1.9 10^-6;
	rhoc = 1.68e-8;
	sigma_cu = 1/rhoc;
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
	req_abs_error = 1e-3;
	req_rel_error = 1e-4;
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
	nodes = cat(nodesv[1:end-1,:], nodesh, dims=1);
	nn = size(nodes)[1];

	#create images
	images = Array{Electrode}(undef, ns);
	for i=1:ns
		start_point = [electrodes[i].start_point[1],
		               electrodes[i].start_point[2],
					   -electrodes[i].start_point[3]];
	    end_point = [electrodes[i].end_point[1],
		             electrodes[i].end_point[2],
	   		     	 -electrodes[i].end_point[3]];
	    r = electrodes[i].radius;
		zi = electrodes[i].zi;
		images[i] = new_electrode(start_point, end_point, r, zi);
	end

	mA, mB = incidence(electrodes, nodes);
	mAT = transpose(mA);
	mBT = transpose(mB);
	#yn = zeros(Complex{Float64}, (nn,nn));
	zl = zeros(Complex{Float64}, (ns,ns));
	zt = zeros(Complex{Float64}, (ns,ns));
	exci = zeros(Complex{Float64}, (nn));
	vout = zeros(Complex{Float64}, (nf,nn));

	## Source input
	#loading the files using CSV seems to be the slower part of the code
	#TODO search alternative
	path = "/home/pedro/codigos/HP_HEM/interfaces/julia/";
	input = string(path, "source.txt");
	source = CSV.read(input, header=["t", "V"]);
	source[:,1] = source[:,1]*1e-9;

	input = string(path, "voltageArt.txt");
	vout_art = CSV.read(input, header=["t", "V"]);
	input = string(path, "currentArt.txt");
	iout_art = CSV.read(input, header=["t", "I"]);

	ent_freq = laplace_transform(Vector{ComplexF64}(source.V), Vector{Float64}(source.t), -1.0im*sk);

	## Freq. loop
	for i = 1:nf
		jw = 1.0im*sk[i];
	    kappa = jw*eps0;
	    k1 = sqrt(jw*mu0*kappa);
	    zl, zt = calculate_impedances(electrodes, k1, jw, mur, kappa,
	                                  max_eval, req_abs_error, req_rel_error,
	                                  error_norm, intg_type);
		kappa_cu = sigma_cu + jw*epsr*eps0;
		ref_t = (kappa - kappa_cu)/(kappa + kappa_cu);
		ref_l = ref_t;
		zl, zt = impedances_images(electrodes, images, zl, zt, k1, jw, mur, kappa,
								   ref_l, ref_t, max_eval, req_abs_error,
								   req_rel_error, error_norm, intg_type)
	    yn = mAT*inv(zt)*mA + mBT*inv(zl)*mB;
		yn[1,1] += gf;
	    exci[1] = ent_freq[i]*gf;
	    vout[i,:] = yn\exci;
	end;

	## Time response
	outlow = map(i -> vout[:,1][Int(i+1)]*sigma(i), kk);
	upperhalf = reverse(conj(outlow));
	pop!(upperhalf);
	lowerhalf = outlow;
	pop!(lowerhalf);
	append!(lowerhalf, upperhalf);
	F = lowerhalf;
	f = real(ifft(F));
	outv = map(i -> exp(sc*t[i])/dt*f[i], 1:length(t));
	# ======
	iout = -(vout[:,1] - ent_freq)*gf;
	outlow = map(i -> iout[:,1][Int(i+1)]*sigma(i), kk);
	upperhalf = reverse(conj(outlow));
	pop!(upperhalf);
	lowerhalf = outlow;
	pop!(lowerhalf);
	append!(lowerhalf, upperhalf);
	F = lowerhalf;
	f = real(ifft(F));
	outi = map(i -> exp(sc*t[i])/dt*f[i], 1:length(t));
	return outv, outi, source, vout_art, iout_art, t
end;

outv, outi, source, vout_art, iout_art, t = @time simulate();
plotly()
#pyplot()
display(plot([t*1e9, source.t*1e9, vout_art.t], [outv, source.V, vout_art.V],
     	     xlims = (0, 50), ylims = (0, 80), xlabel="t (ns)", ylabel="V (V)",
			 label=["calculated" "source" "article"],
			 color=["red" "green" "blue"], marker=true, title="Vout"))

display(plot([t*1e9, iout_art.t], [outi, iout_art.I],
		     xlims = (0, 50), ylims = (-0.2, 0.5), xlabel="t (ns)", ylabel="I (A)",
			 label=["calculated" "article"],
			 color=["red" "blue"], marker=true, title="Iout"))
