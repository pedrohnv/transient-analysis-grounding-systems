% Add the images effect to the impedance matrices.
%
% The integration is done using the C Cubature package
% see https://github.com/stevengj/cubature
%
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
%   max_eval: positive integer, optional, default = 200
%       specifies a maximum number of function evaluations (0 for no limit)
%   req_abs_error: positive real, optional, default = 1e-3
%       the absolute error requested (0 to ignore)
%   req_rel_error: positive real, optional, default = 1e-4
%       the relative error requested (0 to ignore)
%   error_norm: Enumeration, optional, default = Error_norm.PAIRED
%       error checking scheme (see Error_norm)
%   intg_type: Enumeration, optional, default = Integration_type.DOUBLE
%       type of integration to be done (see Integration_type)
%   
% Returns
% -------
%   zl: complex matrix
%       longitudinal impedance
%   zt: complex matrix
%       transversal impedance
function [zl, zt] = impedances_images(electrodes, images, zl, zt, gamma, ...
                                      s, mur, kappa, ref_l, ref_t, ...
                                      max_eval, req_abs_error, ...
                                      req_rel_error, error_norm, intg_type);
    if nargin < 15
        intg_type = Integration_type.DOUBLE;
    elseif ~isa(intg_type, 'Integration_type')
        error('intg_type must be of class Integration_type, an Enumeration')
    end
    if nargin < 14
        error_norm = Error_norm.PAIRED;
    elseif ~isa(error_norm, 'Error_norm')
        error('error_norm must be of class Error_norm, an Enumeration')
    end
    if nargin < 13
        req_rel_error = 1e-4;
    end
    if nargin < 12
        req_abs_error = 1e-3;
    end
    if nargin < 11
        max_eval = 200;
    end
    [zl, zt] = Mimpedances_images(electrodes, images, zl, zt, gamma, s, ...
                                  mur, kappa, ref_l, ref_t, max_eval, ...
                                  req_abs_error, req_rel_error, ...
                                  int16(error_norm), int16(intg_type));
end