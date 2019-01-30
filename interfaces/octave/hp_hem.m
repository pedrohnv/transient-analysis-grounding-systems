hem_octave;

function electrodes_p = electrodes_pointer(electrodes)
% electrodes_pointer: (EVIL) Creates a pointer to an array of Electrodes.
%
% The pointer must be %deleted later by the user.
% This funtion is EVIL. Don't dabble with it unless you do not fear
% losing your soul.
    num_electrodes = length(electrodes);
    electrodes_p = new_electrodep(num_electrodes);
    start_point = new_doublep(3);
    end_point = new_doublep(3);
    for el=1:(num_electrodes)
        for k=1:3
            doublep_assign(start_point, k-1, electrodes(el).start_point(k));
            doublep_assign(end_point, k-1, electrodes(el).end_point(k));
        endfor
        electrodep_assign(electrodes_p, el-1, start_point, end_point,
                          electrodes(el).radius, electrodes(el).internal_impedance);
    endfor
    %delete_doublep(start_point);
    %delete_doublep(end_point);
endfunction  % electrodes_pointer


function complex_arrayp = complexp_array_assign(array)
% complexp_array_assign: (EVIL) Creates a pointer to a complex array.
%
% The pointer must be %deleted later by the user.
% This funtion is EVIL. Don't dabble with it unless you do not fear
% losing your soul.
    n = length(array)
    complex_arrayp = new_complexp(n);
    for k=1:n
        complexp_assign(complex_arrayp, k-1, array(k));
        complexp_assign(complex_arrayp, k-1, array(k));
    endfor
endfunction  % complexp_array_assign


function electrode = Electrode(start_point, end_point, radius, zi)
% Electrode: object to represent an Electrode
% attributes:
%    start_point: real vector (x,y,z)
%    end_point: real vector (x,y,z)
%    middle_point: real vector (x,y,z)
%    radius: real scalar
%    internal_impedance: complex scalar
    if (nargin == 0)
        electrode.start_point = [0.0, 0.0, 0.0];
        electrode.end_point = [0.0, 0.0, 0.0];
        electrode.middle_point = [0.0, 0.0, 0.0];
        electrode.radius = 0.0;
        electrode.internal_impedance = 0.0*1i;
    else
        if !(isvector (start_point) && isreal(start_point))
            error ("@Electrode: start_point must be a real vector");
        endif
        if !(isvector (end_point) && isreal(end_point))
            error ("@Electrode: end_point must be a real vector");
        endif
        if !(isscalar (radius) && isreal(radius))
            error ("@Electrode: radius must be a real scalar");
        endif
        if !(isscalar (zi))
            error ("@Electrode: zi must be a complex scalar");
        endif
        electrode.start_point = start_point;
        electrode.end_point = end_point;
        electrode.middle_point = (start_point + end_point)/2.0;
        electrode.radius = radius;
        electrode.internal_impedance = zi;
    endif
endfunction % Electrode


function [zl, zt] = calculate_impedances(electrodes, gamma1, s, mur, kappa,
                                         max_eval, req_abs_error, req_rel_error)
% calculate_impedances: Calculates the impedance matrices ZL and ZT.
%
% parameters
%   electrodes: array of structure Electrode
%   gamma1: complex scalar, medium propagation constant
%   s: complex scalar, complex angular frequency `c + I*w` (rad/s)
%   mur: real scalar, relative magnetic permeability of the medium
%   kappa: complex scalar, medium complex conductivity `(sigma + I*w*eps)` in S/m
%   max_eval: int scalar, specifies a maximum number of function evaluations (0 for no limit)
%   req_abs_error: real scalar, the absolute error requested (0 to ignore)
%   req_rel_error: real scalar, the relative error requested (0 to ignore)
    error_norm = 2; % cubature: PAIRED
    integration_type = 1; % DOUBLE
    num_electrodes = length(electrodes);
    electrodes_p = electrodes_pointer(electrodes);
    zl_p = new_complexp(num_electrodes);
    zt_p = new_complexp(num_electrodes);
    failure = Ocalculate_impedances(electrodes_p, num_electrodes, zl_p, zt_p,
                                    gamma1, s, mur, kappa, max_eval, req_abs_error,
                                    req_rel_error, error_norm, integration_type);
    zl = zeros(num_electrodes);
    zt = zeros(num_electrodes);
    for m=0:(num_electrodes - 1)
        for k=0:(num_electrodes - 1)
            zl(m+1, k+1) = complexp_value(zl_p, m*num_electrodes + k);
            zt(m+1, k+1) = complexp_value(zt_p, m*num_electrodes + k);
        endfor
    endfor
    %delete_electrodep(electrodes_p);
    %delete_complexp(zl_p);
    %delete_complexp(zt_p);
endfunction % calculate_impedances


function [zl, zt] = impedances_images(zl, zt, ref_l, ref_t, electrodes, images,
                                      gamma1, s, mur, kappa, max_eval,
                                      req_abs_error, req_rel_error)
% impedances_images: Updates the impedance matrices ZL and ZT with the images' influence.
%
% parameters
%   zl: matrix of longitudinal impedances
%   zt: matrix of transversal impedances
%   ref_l: longitudinal reflection coefficient
%   ref_l: transversal reflection coefficient
%   electrodes: array of structure Electrode
%   images: array of structure Electrode
%   gamma1: complex scalar, medium propagation constant
%   s: complex scalar, complex angular frequency `c + I*w` (rad/s)
%   mur: real scalar, relative magnetic permeability of the medium
%   kappa: complex scalar, medium complex conductivity `(sigma + I*w*eps)` in S/m
%   max_eval: int scalar, specifies a maximum number of function evaluations (0 for no limit)
%   req_abs_error: real scalar, the absolute error requested (0 to ignore)
%   req_rel_error: real scalar, the relative error requested (0 to ignore)
    error_norm = 2; % cubature: PAIRED
    integration_type = 1; % DOUBLE
    num_electrodes = length(electrodes);
    electrodes_p = electrodes_pointer(electrodes);
    images_p = electrodes_pointer(images);
    zl_p = complexp_array_assign(zl);
    zt_p = complexp_array_assign(zt);
    failure = Oimpedances_images(electrodes_p, images_p, num_electrodes,
                                 zl_p, zt_p, gamma1, s, mur, kappa,
                                 ref_l, ref_t, max_eval, req_abs_error,
                                 req_rel_error, error_norm, integration_type);
    num_nodes = length(zl);
    for m=0:(num_nodes - 1)
        for k=0:(num_nodes - 1)
            zl(m+1, k+1) = complexp_value(zl_p, m*num_electrodes + k);
            zt(m+1, k+1) = complexp_value(zt_p, m*num_electrodes + k);
        endfor
    endfor
    %delete_electrodep(electrodes_p);
    %delete_electrodep(images_p);
    %delete_complexp(zl_p);
    %delete_complexp(zt_p);
endfunction % impedances_images


function zh = harmonic_impedance1(s, kappa1, kappa2, gamma1, electrodes, images,
                                  nodes, max_eval, req_abs_error, req_rel_error,
                                  rsource)
% harmonic_impedance1:
%       Calculates the harmonic impedance of a copper electrode system in a
%       two layer medium. Electrodes are considered to be in medium 1.
%       No segmentation of the electrodes is done.
%       Injection node is considered the first.
%
% parameters
%   s: complex array, complex angular frequency `c + I*w` (rad/s)
%   kappa1: complex array, medium 1 complex conductivity `(sigma + I*w*eps)` in S/m
%   kappa2: complex array, medium 2 complex conductivity `(sigma + I*w*eps)` in S/m
%   gamma1: complex array, medium 1 propagation constant
%   electrodes: array of structure Electrode
%   images: array of structure Electrode
%   nodes: matrix (num_nodes, 3) of nodes
%   max_eval: int scalar, specifies a maximum number of function evaluations (0 for no limit)
%   req_abs_error: real scalar, the absolute error requested (0 to ignore)
%   req_rel_error: real scalar, the relative error requested (0 to ignore)
%   rsource: source internal resistence to consider.
    error_norm = 2; % cubature: PAIRED
    integration_type = 1; % DOUBLE
    num_nodes = length(nodes);
    num_electrodes = length(electrodes);
    ns = length(s);
    electrodes_p = electrodes_pointer(electrodes);
    images_p = electrodes_pointer(images);
    s_p = complexp_array_assign(s);
    kappa1_p = complexp_array_assign(kappa1);
    kappa2_p = complexp_array_assign(kappa2);
    gamma1_p = complexp_array_assign(gamma1);
    zh_p = new_complexp(ns);
    failure = Oharmonic_impedance1(ns, s_p, kappa1_p, kappa2_p, gamma1_p,
                                   electrodes_p, images_p, num_electrodes,
                                   nodes_p, num_nodes, max_eval, req_abs_error,
                                   req_rel_error, error_norm, rsource, zh_p);
    zh = zeros(ns);
    for k=0:(ns - 1)
        zh(k+1) = complexp_value(zh_p, k);
    endfor
    %delete_electrodep(electrodes_p);
    %delete_electrodep(images_p);
    %delete_complexp(s_p);
    %delete_complexp(kappa1_p);
    %delete_complexp(kappa2_p);
    %delete_complexp(gamma1_p);
    %delete_complexp(zh_p);
endfunction  % harmonic_impedance1


% =================
n = 10;
z = new_complexp(n);
for k=0:(n-1)
    complexp_assign(z, k, 0.2*k + (k - 3.0)*1i);
end
sum = sumall(z, n)
%delete_complexp(z);
complexp_value(z, 2) % pointer not deallocated OR undefined behavior
% =================
elarr = [
    Electrode([0,0,0], [1,0,0], 7e-3, 0),
    Electrode([0,0,0], [0,1,0], 7e-3, 0)
];
[zl, zt] = calculate_impedances(elarr, 1.0, 1.0, 1.0, 1.0, 200, 1e-3, 1e-4);
zl
zt
length(zl)
disp("END: success")
% FIXME call to pcubature leads to segfault on call
% FIXME call to hcubature leads to segfault on exit
% FIXME: deleting the pointers inside functions lead to segfault
