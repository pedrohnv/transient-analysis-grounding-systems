function [zl, zt] = Mimpedances_images(electrodes, images, zl, zt, gamma, ...
                                       s, mur, kappa, ref_l, ref_t, ...
                                       req_abs_error, req_rel_error, ...
                                       intg_type)
    FOUR_PI = 12.56637061435917;
    MU0 = 1.256637061435917e-6; %permeability vac.
    iwu_4pi = s*mur*MU0/(FOUR_PI);
    one_4pik = 1.0/(FOUR_PI*kappa);
    ns = length(electrodes);
    k1 = [0., 0., 0.];
    for i=1:ns
        ls = electrodes(i).length;
        for m=1:3
            k1(m) = (electrodes(i).end_point(m) - electrodes(i).start_point(m));
        end
        for k=i:ns
            lr = electrodes(k).length;
            cost = 0.0;
            for m=1:3
                k2 = (electrodes(k).end_point(m) - electrodes(k).start_point(m));
                cost = cost + k1(m)*k2;
            end
            cost = abs(cost/(ls*lr));
            intg = integral_hem(electrodes(i), images(k), gamma, req_abs_error, req_rel_error, intg_type);
            zl(k,i) = zl(k,i) + ref_l*iwu_4pi*intg*cost;
            zt(k,i) = zt(k,i) + ref_t*one_4pik/(ls*lr)*intg;
            zl(i,k) = zl(k,i);
            zt(i,k) = zt(k,i);
        end
    end
end