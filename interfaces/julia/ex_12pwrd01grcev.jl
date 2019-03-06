#=
Test case 12pwrd01grcev-a

Reproducing the results in [1] for a grounding grid.
Execution time is expected to be up to 20 min.

[1] L. D. Grcev and M. Heimbach, "Frequency dependent and transient
characteristics of substation grounding systems," in IEEE Transactions on
Power Delivery, vol. 12, no. 1, pp. 172-178, Jan. 1997.
doi: 10.1109/61.568238
=#
using Plots

include("hp_hem.jl");

function simulate(gs::Int, nf::Int)
	## Parameters
	# Soil
	mu0 = pi*4e-7;
	mur= 1.0;
	eps0 = 8.854e-12;
	epsr = 10;
	sigma1 = 1.0/1000.0;

	# Integration Parameters
	max_eval = 200;
	req_abs_error = 1e-3;
	req_rel_error = 1e-4;
	error_norm = 1; #paired
	intg_type = 1; #double integral

	# Frequencies
	freq = exp10.(range(2, stop=6.4, length=nf)); #logspace
	omega = 2*pi*freq[nf];
	lambda = (2*pi/omega)*(1/sqrt( epsr*eps0*mu0/2*(1 + sqrt(1 + (sigma1/(omega*epsr*eps0))^2)) ));
	frac = lambda/10; #for segmentation

	# Grid
	r = 7e-3;
	h = -0.5;
	l = gs;
	n = Int(gs/10);
	elecs = electrode_grid(l, n, l, n, h, r);
	electrodes, nodes =  seg_electrode_list(elecs, frac);
	num_electrodes= ns = length(electrodes)
	num_nodes = size(nodes)[1]
	inj_node = matchrow([0.,0.,h], nodes)

	#create images
	images = Array{Electrode}(undef, ns);
	for i=1:ns
		images;
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
	exci = zeros(ComplexF64, num_nodes);
	exci[inj_node] = 1.0;
	zh = zeros(ComplexF64, nf);

	## Frequency loop
	for i = 1:nf
		#println("i = ", i)
		jw = 1.0im*2*pi*freq[i];
	    kappa = sigma1 + jw*epsr*eps0;
	    k1 = sqrt(jw*mu0*kappa);
	    zl, zt = calculate_impedances(electrodes, k1, jw, mur, kappa,
	                                  max_eval, req_abs_error, req_rel_error,
	                                  error_norm, intg_type);
		kappa_air = jw*eps0;
		ref_t = (kappa - kappa_air)/(kappa + kappa_air);
		ref_l = ref_t;
		zl, zt = impedances_images(electrodes, images, zl, zt, k1, jw, mur, kappa,
								   ref_l, ref_t, max_eval, req_abs_error,
								   req_rel_error, error_norm, intg_type)
	    yn = mAT*inv(zt)*mA + mBT*inv(zl)*mB;
		vn = yn\exci;
	    zh[i] = vn[inj_node];
	end;
	return zh
end;

nf = 150;
freq = exp10.(range(2, stop=6.4, length=nf)); #logspace
zh = zeros(ComplexF64, (nf, 5)); #Julia is Column Major.
i = 1;
tot_time = @timed for gs in [10, 20, 30, 60, 120]
	println("gs = ", gs)
	global zh[:,i] = @time simulate(gs, nf);
	global i += 1;
end;
println("total elapsed time: ", tot_time[2]/60, " min")

plotly()
#pyplot()
display(plot(freq, [map(i -> abs(i), zh[:,k]) for k=1:5],
             xaxis=:log, xlabel="f (Hz)", ylabel="|Zin| (Ohms)",
			 title="Harmonic Impedance",
			 label=["GS 10" "GS 20" "GS 30" "GS 60" "GS 120"]))
display(plot(freq, [map(i -> angle(i)*180/pi, zh[:,k]) for k=1:5],
             xaxis=:log, xlabel="f (Hz)", ylabel="phase Zin (deg)",
 			 title="Harmonic Impedance",
 			 label=["GS 10" "GS 20" "GS 30" "GS 60" "GS 120"]))
