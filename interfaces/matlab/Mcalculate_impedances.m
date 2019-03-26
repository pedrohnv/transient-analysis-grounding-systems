function [zl, zt] = Mcalculate_impedances(electrodes, gamma, s, mur, kappa, ...
                                          req_abs_error, req_rel_error, ...
                                          intg_type)
    ns = length(electrodes);
    zl = zeros(ns, ns, 'like', (0.0 + 0.0j));
    zt = zl;
    FOUR_PI = 12.56637061435917;
    MU0 = 1.256637061435917e-6; %permeability vac.
    iwu_4pi = s*mur*MU0/(FOUR_PI);
    one_4pik = 1.0/(FOUR_PI*kappa);
    k1 = [0., 0., 0.];
    for i=1:ns
        ls = electrodes(i).length;
        k1(1) = electrodes(i).radius/ls;
        k2 = sqrt(1.0 + k1(1)*k1(1));
        cost = 2.0*(log( (k2 + 1.)/k1(1) ) - k2 + k1(1));
        zl(i,i) = iwu_4pi*ls*cost + electrodes(i).zi;
        zt(i,i) = one_4pik/ls*cost;
        for m=1:3
            k1(m) = (electrodes(i).end_point(m) - electrodes(i).start_point(m));
        end
        for k=(i+1):ns
            lr = electrodes(k).length;
            cost = 0.0;
            for m=1:3
                k2 = (electrodes(k).end_point(m) - electrodes(k).start_point(m));
                cost = cost + k1(m)*k2;
            end
            cost = abs(cost/(ls*lr));
            intg = integral_hem(electrodes(i), electrodes(k), gamma, req_abs_error, req_rel_error, intg_type);
            zl(k,i) = iwu_4pi*intg*cost;
            zt(k,i) = one_4pik/(ls*lr)*intg;
            zl(i,k) = zl(k,i);
            zt(i,k) = zt(k,i);
        end
    end
end