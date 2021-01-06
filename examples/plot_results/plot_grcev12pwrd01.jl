# Plot the results from example grcev12pwrd01
using Plots
using LinearAlgebra

function load_vector(fname)
	stringvec = split(read(fname, String), "\n")
	try
		return map(x -> parse(ComplexF64, x), stringvec)
	catch
		return map(x -> parse(ComplexF64, x), stringvec[1:end-1])
	end
end

begin
	p1 = plot(xlabel="Frequency [Hz]", ylabel="Harmonic Impedance [Ω]", legend=:topleft)
	p2 = plot(xlabel="Frequency [Hz]", ylabel="Phase [deg]", legend=:topleft)
	for gs in [10, 20, 30, 60, 120]
		zh = load_vector("gs$(gs).csv")
		nf = length(zh)
		freq = exp10.(range(2, 6.4, length=nf))
		plot!(p1, freq, abs.(zh), xaxis=:log, label="GS $(gs)")
		plot!(p2, freq, rad2deg.(angle.(zh)), xaxis=:log, label="GS $(gs)")
	end
	display(p1)
	display(p2)
end
