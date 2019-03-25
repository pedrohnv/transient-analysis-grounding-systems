function [a, b] = incidence(electrodes, nodes)
	% build incidence matrices for calculating 'YN = AT*inv(zt)*A + BT*inv(zl)*B'
    ns = length(electrodes);
    dim = size(nodes);
    nn = dim(1);
    a = zeros(ns,nn);
    b = zeros(ns,nn);
    for k = 1:nn
        noincidence = true;
        for i = 1:ns
            if vecnorm(electrodes(i).start_point - nodes(k,:)) < eps
                a(i,k) = 0.5;
                b(i,k) = 1.0;
                noincidence = false;
            elseif vecnorm(electrodes(i).end_point - nodes(k,:)) < eps
                a(i,k) = 0.5;
                b(i,k) = -1.0;
                noincidence = false;
            end
        end
        if noincidence
            error('no incidence on nodes(%i,:)', k)
        end
    end
end