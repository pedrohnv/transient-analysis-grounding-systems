% Test case 17pwrd03noda
%
% Reproducing the results in [1] of a time domain surge response.
%
% [1] Noda, Taku, and Shigeru Yokoyama. "Thin wire representation in finite difference time domain
% surge simulation." IEEE Transactions on Power Delivery 17.3 (2002): 840-847.

% setup() % To compile MEX files

%% Parameters
mu0 = pi*4e-7;
mur = 1.0;
eps0 = 8.854e-12;
epsr = 1.0;
sigma1 = 0.0 + 0.0j;
%rhoc = 1.9 10^-6;
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

% 	Due to numerical erros, to smooth the response, its necessary to use a
% 	final time much greater than that up to which is desired.
T = 0.7e-7*2;
dt = 2.0e-9;
n = T/dt;
t = 0.0:dt:(T-dt);
sc = log(n^2)/T;
kk = 0:1:n/2;
dw = 2.0*pi/(n*dt);
alpha = 0.53836;
sigma = @(j) alpha + (1 - alpha)*cos(2*pi*j/n);
sk = transpose(-1j*sc*ones(length(kk),1)) + dw*kk;
nf = length(sk);
freq = real(sk)/(2*pi);
omega = 2*pi*freq(nf);
lambda = (2*pi/omega)*(1/sqrt( epsr*eps0*mu0/2*(1 + sqrt(1 + (sigma1/(omega*epsr*eps0))^2)) ));

% Integration Parameters
max_eval = 200;
req_abs_error = 1e-3;
req_rel_error = 1e-4;
error_norm = Error_norm.PAIRED;
intg_type = Integration_type.DOUBLE;

%% Electrodes
x = 10;
nv = ceil(h/(lambda/x));
nh = ceil(l/(lambda/x));
vertical = new_electrode([0, 0, 0], [0, 0, 0.5], 10e-3, 0.0);
horizontal = new_electrode([0, 0, 0.5], [4.0, 0, 0.5], 15e-3, 0.0);
[elecv, nodesv] = segment_electrode(vertical, nv);
[elech, nodesh] = segment_electrode(horizontal, nh);
electrodes = cat(1, elecv, elech);
ns = length(electrodes);
nodes = cat(1, nodesv(1:end-1,:), nodesh);
nn = length(nodes);
num_nodes = nn;

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
inj_node = 1;
exci(inj_node) = 1.0;
zh = zeros(nf, 1, 'like', 1+1j);
zl = zeros(ns, ns, 'like', 1+1j);
zt = zl;

vout = zeros(nf, nn, 'like', zh);

%% Source input
if ispc % windows?
    path = '..\\..\\examples\\noda17pwrd03_auxfiles\\';
else
    path = '../../examples/noda17pwrd03_auxfiles/';
end

fileID = fopen([path, 'source.txt'],'r');
A = fscanf(fileID,'%f,%f');
source = transpose(reshape(A, 2,12));
source(:,1) = source(:,1)*1e-9;

fileID = fopen([path, 'voltage.txt'],'r');
A = fscanf(fileID,'%f,%f');
vout_art = transpose(reshape(A, 2,37));

fileID = fopen([path, 'current.txt'],'r');
A = fscanf(fileID,'%f,%f');
iout_art = transpose(reshape(A, 2,39));

ent_freq = laplace_transform(source(:,2), source(:,1), -1.0j*sk);

%% Freq. loop
for i = 1:nf
    jw = 1.0j*sk(i);
    kappa = jw*eps0;
    k1 = sqrt(jw*mu0*kappa);
    kappa_cu = sigma_cu + jw*epsr*eps0;
    ref_t = (kappa - kappa_cu)/(kappa + kappa_cu) + 0.0j;
    ref_l = ref_t + 0.0j;
    [zl, zt] = calculate_impedances(electrodes, k1, jw, mur, kappa, ...
                                    max_eval, req_abs_error, req_rel_error, ...
                                    error_norm, intg_type);
    [zl, zt] = impedances_images(electrodes, images, zl, zt, k1, jw, mur, kappa, ...
                                 ref_l, ref_t, max_eval, req_abs_error, ...
                                 req_rel_error, error_norm, intg_type);
    yn = mAT*inv(zt)*mA + mBT*inv(zl)*mB;
    yn(inj_node,inj_node) = yn(inj_node,inj_node) + gf;
    exci(inj_node) = ent_freq(i)*gf;
    vout(i,:) = yn\exci;
end

%% Time response
for m = length(kk):-1:1
    i = kk(m);
    outlow(m) = vout((i+1), 1)*sigma(i);
end
upperhalf = fliplr(conj(outlow));
upperhalf = upperhalf(1:end-1);
lowerhalf = outlow;
lowerhalf = lowerhalf(1:end-1);
lowerhalf = [lowerhalf, upperhalf];
F = lowerhalf;
f = real(ifft(F));
for i = length(t):-1:1
    outv(i) = exp(sc*t(i))/dt*f(i);
end
% ======
iout = -(vout(:,1) - ent_freq)*gf;
for m = length(kk):-1:1
    i = kk(m);
    outlow(m) = iout((i+1), 1)*sigma(i);
end
upperhalf = fliplr(conj(outlow));
upperhalf = upperhalf(1:end-1);
lowerhalf = outlow;
lowerhalf = lowerhalf(1:end-1);
lowerhalf = [lowerhalf, upperhalf];
F = lowerhalf;
f = real(ifft(F));
for i = length(t):-1:1
    outi(i) = exp(sc*t(i))/dt*f(i);
end
%% PLOT
figure()
plot(t*1e9, outv)
hold on
plot(vout_art(:,1), vout_art(:,2), 'r')
hold on
plot(source(:,1)*1e9, source(:,2), 'g')
ylabel('V (V)')
xlabel('t (ns)')
xlim([0, 50])
ylim([0, 80])
legend('calculated', 'source', 'article')

figure()
plot(t*1e9, outi)
hold on
plot(iout_art(:,1), iout_art(:,2), 'r')
ylabel('I (A)')
xlabel('t (ns)')
xlim([0, 50])
ylim([-0.2, 0.5])
legend('calculated', 'article')
