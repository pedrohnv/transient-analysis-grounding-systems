% Add the images effect to the impedance matrices.
%
% The integration is done using the C Cubature package
% see https://github.com/stevengj/cubature
% 
% Calling syntax:
% [zl, zt] = impedances_images(electrodes, images, zl, zt, gamma, s, mur, ...
%                              kappa, ref_l, ref_t, max_eval, req_abs_error, ...
%                              req_rel_error, error_norm, intg_type);
% Parameters
% ----------
%   electrodes: array of Electrode struct
%       see new_electrode
%   images: array of Electrode struct
%       see new_electrode
%   zl: complex matrix
%       longitudinal impedance
%   zt: complex matrix
%       transversal impedance
%   gamma: complex
%       medium propagation constant
%   s: complex
%       angular frequency `c + I*w` (rad/s)
%   mur: real
%       relative magnetic permeability of the medium
%   kappa: complex
%       medium conductivity `(sigma + I*w*eps)` in S/m
%   ref_l: complex
%       longitudinal current reflection coefficient
%   ref_t: complex
%       transversal current reflection coefficient
%   max_eval: integer
%       specifies a maximum number of function evaluations (0 for no limit)
%   req_abs_error: real
%       the absolute error requested (0 to ignore)
%   req_rel_error: real
%       the relative error requested (0 to ignore)
%   error_norm: Enumeration
%       error checking scheme (see Error_norm)
%   integration_type: Enumeration
%       type of integration to be done (see Integration_type)
%   
% Returns
% -------
%   zl: complex matrix
%       longitudinal impedance
%   zt: complex matrix
%       transversal impedance
