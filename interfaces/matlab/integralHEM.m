% exp(-gamma*r)/r
function intg = integralHEM(sender, receiver, gamma, req_abs_error, req_rel_error, intg_type)
    if nargin < 6
        intg_type = 2;
    end
    if nargin < 5
        req_rel_error = 1e-4;
    end
    if nargin < 4
        req_abs_error = 1e-3;
    end
    ls = sender.length;
    lr = receiver.length;
    point_s = @(ts) transpose(ts/ls).*(sender.end_point - sender.start_point) + sender.start_point;
    point_r = @(tr) transpose(tr/lr).*(receiver.end_point - receiver.start_point) + receiver.start_point;
    r = @(ts,tr) transpose(vecnorm(point_s(ts) - point_r(tr), 2, 2));
    if intg_type == 0
        rbar = vecnorm(sender.middle_point - receiver.middle_point, 2, 2);
        intg = ls*lr*exp(-gamma*rbar)/rbar;
    elseif intg_type == 1
        % FIXME error
        fun = @(ts,tr) exp(-gamma*r(ts,tr))./r(ts,tr);
        intg = integral2(fun, 0.0, ls, 0.0, lr, ...
                         'AbsTol', req_abs_error, 'RelTol', req_rel_error);
    elseif intg_type == 2
        r1 = @(tr) transpose(vecnorm(sender.start_point - point_r(tr), 2, 2));
        r2 = @(tr) transpose(vecnorm(sender.end_point - point_r(tr), 2, 2));
        fun = @(tr) log( (r1(tr) + r2(tr) + ls)./(r1(tr) + r2(tr) - ls) );
        intg = integral(fun, 0.0, lr, 'AbsTol', req_abs_error, 'RelTol', req_rel_error);
        rbar = vecnorm(sender.middle_point - receiver.middle_point, 2, 2);
        intg = exp(-gamma*rbar)*intg;
    end
end