using LinearAlgebra
using DataFrames
using CSV
import Printf.@sprintf
using Plots; pyplot()

## Electric Field
for rho in [100, 1000]
    df = DataFrame(CSV.File("alipio83powsys_efield_$(rho).csv"))
    f = sort(unique(df.f))
    nf = length(f)
    x = sort(unique(df.x))

    p = plot(ylabel="|Ex| [V/m]", xlabel="x [m]", title="σ = $(1/rho * 1e3) mS", legend=:best)
    labels = ["100 Hz", "500 kHz", "1 MHz", "2 MHz"]
    colors = [:blue, :orange, :green, :black]
    for k = 1:nf
        c = isapprox.(df.f, f[k])
        efc = (parse.(ComplexF64, df.ecx1[c]))
        efnonc = (parse.(ComplexF64, df.encx1[c]))
        #plot!(x, abs.(efc), label="EC")
        #plot!(x, abs.(efnonc), label="ENC")
        plot!(x, abs.(efc + efnonc), label=labels[k], linecolor=colors[k])
        #plot!(x, abs.(efnonc ./ efc), label=labels[k])
    end
    display(p)
    png(p, "alipio83powsys_efield_$(rho)")
end



## Step Voltage along the line y = 0, +x direction
for rho in [100, 1000]
    freq = ["100 Hz", "500 kHz", "1 MHz", "2 MHz"]
    for k = 1:4
        title = @sprintf("f = %s, σ = %.1f mS", freq[k], 1/rho * 1e3)
        # GPD difference =====================================================
        df = DataFrame(CSV.File("alipio83powsys_gpd_$(rho).csv"))
        f = sort(unique(df.f))
        nf = length(f)
        x = sort(unique(df.x))
        dx = x[2] - x[1]
        c = isapprox.(f[k], df.f)
        N = Int(1 / dx)
        nx = length(x) - N
        v = parse.(ComplexF64, df.i1)[c]
        dv = [(v[i] - v[i + N]) for i = 1:nx]
        x = x[1:end-N]
        p = plot(x, abs.(dv), title=title, ylabel="Step Voltage [V]",
                 xlabel="x [m]", label="pot. dif.")

        # integral of the E_tot field ========================================
#=
        df = DataFrame(CSV.File("alipio83powsys_efield_$(rho).csv"))
        f = sort(unique(df.f))
        nf = length(f)
        x = sort(unique(df.x))
        dx = x[2] - x[1]
        c = isapprox.(df.f, f[k])
        ef = parse.(ComplexF64, df.ecx1[c]) + parse.(ComplexF64, df.encx1[c])
        N = Int(1 / dx)
        nx = length(x) - N
        dv = zeros(ComplexF64, nx)
        for i = 1:nx
            s = sum(ef[(i + 1):(i + N - 1)])
            dv[i] = dx * (s + (ef[i] + ef[i + N]) / 2)
        end
        x = x[1:end-N]
        plot!(x, abs.(dv), label="voltage")
=#
        # Voltage ============================================================
        df = DataFrame(CSV.File("alipio83powsys_voltage_$(rho).csv"))
        c = isapprox.(df.f, f[k])
        dv = (parse.(ComplexF64, df.v1[c]))
        plot!(df.x[c], abs.(dv), label="voltage")

        display(p)
        png(p, "alipio83powsys_voltage$(k)_$(rho)")
    end
end
