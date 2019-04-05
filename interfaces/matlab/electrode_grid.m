% Creates an electrode grid `h` coordinate below ground with each conductor
% having radius `r` and internal impedance `zi=0.0`.
% The grid has dimensions `a*b` with `n` and `m` divisions respectively.
function electrodes = electrode_grid(a, n, b, m, h, r)
    xx = 0:a/n:a;
    yy = 0:b/m:b;
    num_elec = n*(m + 1) + m*(n + 1);
    electrodes = repmat(new_electrode(), num_elec, 1); %pre-allocate
    e = 1;
    for k=1:(m+1)
		for i=1:n
			electrodes(e) = new_electrode([xx(i), yy(k), h], [xx(i+1), yy(k), h], r);
			e = e + 1;
		end
	end
	for k=1:(n+1)
		for i=1:m
			electrodes(e) = new_electrode([xx(k), yy(i), h], [xx(k), yy(i+1), h], r);
			e = e + 1;
		end
    end
end