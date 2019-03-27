% Computes the numerical integration of integral2(exp(-gamma*r)/r, dls, dlr)
% using pure MATLAB routines.
% Parameters
% ----------
%   sender: Electrode struct
%       see new_electrode
%   receiver: Electrode struct
%       see new_electrode
%   gamma: complex
%       medium propagation constant
%   req_abs_error: real
%       the absolute error requested (0 to ignore)
%   req_rel_error: real
%       the relative error requested (0 to ignore)
%   intg_type: Enumeration
%       type of integration to be done (see Integration_type)
%   
% Returns
% -------
%   intg: complex
%       integral result
function intg = integral_hem(sender, receiver, gamma, req_abs_error, req_rel_error, intg_type)
    if nargin < 6
        intg_type = Integration_type.INTG_LOGNF;
    end
    if nargin < 5
        req_rel_error = 1e-4;
    end
    if nargin < 4
        req_abs_error = 1e-3;
    end
    ls = sender.length;
    lr = receiver.length;
    point_s = @(ts) transpose(ts).*(sender.end_point - sender.start_point) + sender.start_point;
    point_r = @(tr) transpose(tr).*(receiver.end_point - receiver.start_point) + receiver.start_point;
    r = @(ts,tr) transpose(vecnorm(point_s(ts) - point_r(tr), 2, 2));
    if intg_type == Integration_type.INTG_NONE
        rbar = vecnorm(sender.middle_point - receiver.middle_point, 2, 2);
        intg = ls*lr*exp(-gamma*rbar)/rbar;
    elseif intg_type == Integration_type.INTG_DOUBLE
        % FIXME error
        fun = @(ts,tr) exp(-gamma*r(ts,tr))./r(ts,tr);
        intg = integral2(fun, 0.0, 1.0, 0.0, 1.0, ...
                         'AbsTol', req_abs_error, 'RelTol', req_rel_error);
        intg = intg*ls*lr;
    elseif intg_type == Integration_type.INTG_EXP_LOGNF
        r1 = @(tr) transpose(vecnorm(sender.start_point - point_r(tr), 2, 2));
        r2 = @(tr) transpose(vecnorm(sender.end_point - point_r(tr), 2, 2));
        rbar = @(tr) transpose(vecnorm(sender.middle_point - point_r(tr), 2, 2));
        fun = @(tr) exp(-gamma*rbar(tr)).*log( (r1(tr) + r2(tr) + ls)./(r1(tr) + r2(tr) - ls) );
        intg = integral(fun, 0.0, 1.0, 'AbsTol', req_abs_error, 'RelTol', req_rel_error);
        intg = intg*ls;
    elseif intg_type == Integration_type.INTG_LOGNF
        r1 = @(tr) transpose(vecnorm(sender.start_point - point_r(tr), 2, 2));
        r2 = @(tr) transpose(vecnorm(sender.end_point - point_r(tr), 2, 2));
        fun = @(tr) log( (r1(tr) + r2(tr) + ls)./(r1(tr) + r2(tr) - ls) );
        intg = integral(fun, 0.0, 1.0, 'AbsTol', req_abs_error, 'RelTol', req_rel_error);
        rbar = vecnorm(sender.middle_point - receiver.middle_point, 2, 2);
        intg = exp(-gamma*rbar)*intg*ls;
    end
end