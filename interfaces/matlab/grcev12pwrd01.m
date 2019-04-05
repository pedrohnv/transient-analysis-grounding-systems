% Test case 12pwrd01grcev-a
% 
% Reproducing the results in [1] for a grounding grid.
% Execution time is expected to be up to 20 min.
% 
% [1] L. D. Grcev and M. Heimbach, "Frequency dependent and transient
% characteristics of substation grounding systems," in IEEE Transactions on
% Power Delivery, vol. 12, no. 1, pp. 172-178, Jan. 1997.
% doi: 10.1109/61.568238

% setup() % To compile MEX files

%% Parameters
gs = 30 % grid size in meters
% gs = 60 -> 30 seconds
% gs = 120 -> 685 seconds
% Soil
mu0 = pi*4e-7;
mur = 1.0;
eps0 = 8.854e-12;
epsr = 10;
sigma1 = 1.0/1000.0;

% Integration Parameters
max_eval = 200;
req_abs_error = 1e-3;
req_rel_error = 1e-4;
error_norm = Error_norm.PAIRED;
intg_type = Integration_type.DOUBLE;

% Frequencies
nf = 150;
freq = logspace(2, 6.4, nf);
omega = 2*pi*freq(nf);
lambda = (2*pi/omega)*(1/sqrt( epsr*eps0*mu0/2*(1 + sqrt(1 + (sigma1/(omega*epsr*eps0))^2)) ));
frac = lambda/6; %for segmentation

% Grid
r = 7e-3;
h = -0.5;
l = gs;
n = gs/10;
elecs = electrode_grid(l, n, l, n, h, r);
[electrodes, nodes] =  seg_electrode_list(elecs, frac);
num_electrodes = length(electrodes)
ns = num_electrodes;
dims = size(nodes);
num_nodes = dims(1)
for i=1:num_nodes
    if ismember([0.,0.,h], nodes(i,:), 'row')
       inj_node = i;
       break
    end
end
images = electrodes;
for i=1:ns
    images(i).start_point(3) = -images(i).start_point(3);
    images(i).middle_point(3) = -images(i).middle_point(3);
    images(i).end_point(3) = -images(i).end_point(3);
end

[mA, mB] = incidence(electrodes, nodes);
mAT = transpose(mA);
mBT = transpose(mB);
exci = zeros(num_nodes,1);
exci(inj_node) = 1.0;
zh = zeros(nf, 1, 'like', 1+1j);
zl = zeros(ns, ns, 'like', 1+1j);
zt = zl;

%% Frequency loop
tic
for i = 1:nf
    jw = 1.0j*2*pi*freq(i);
    kappa = sigma1 + jw*epsr*eps0;
    k1 = sqrt(jw*mu0*kappa);
    kappa_air = jw*eps0;
    ref_t = (kappa - kappa_air)/(kappa + kappa_air);
    ref_l = ref_t;
    [zl, zt] = calculate_impedances(electrodes, k1, jw, mur, kappa, ...
                                    max_eval, req_abs_error, req_rel_error, ...
                                    error_norm, intg_type);
    [zl, zt] = impedances_images(electrodes, images, zl, zt, k1, jw, mur, kappa, ...
                                 ref_l, ref_t, max_eval, req_abs_error, ...
                                 req_rel_error, error_norm, intg_type);
    yn = mAT*inv(zt)*mA + mBT*inv(zl)*mB;
    vn = yn\exci;
    zh(i) = vn(inj_node);
end
toc
%% PLOT
figure()
semilogx(freq, abs(zh))
title('Harmonic Impedance')
ylabel('|Z| (\Omega)')
xlabel('f (Hz)')

figure()
semilogx(freq, angle(zh)*180/pi)
title('Harmonic Impedance')
ylabel('phase(Z) (deg)')
xlabel('f (Hz)')