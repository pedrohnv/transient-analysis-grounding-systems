% Numerical Laplace Transform of an arbitrary function f(t).
% 
% Parameters
% ----------
%   f : vector of the sampled function in time domain
%   t : vector of the time of the samples
%   s : vector of "frequency" points of interest
% 
% Returns
% -------
%   F(s) : vector of the sampled function in "frequency" domain
function res = laplace_transform(f, t, s)
    nt = length(t);
    nf = length(s);
    res = zeros(nf, 1, 'like', 0.0+1.0j);
    for k = 1:nf
        for i = 1:(nt-1)
            e0 = exp(s(k)*t(i));
            e1 = exp(s(k)*t(i+1));
            dt = t(i+1) - t(i);
            x = (f(i+1) - f(i))/s(k);
            res(k) = res(k) + (e1*(f(i+1)*dt - x) - e0*(f(i)*dt - x))/dt;
        end
        res(k) = 1/s(k)*res(k);
    end
end