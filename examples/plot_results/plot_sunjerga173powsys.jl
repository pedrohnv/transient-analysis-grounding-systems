using LinearAlgebra
using DataFrames
using CSV
import Printf.@sprintf
using Plots; pyplot()

function heidler(t, imax, tau1, tau2, n)
    xi = exp( -(tau1 / tau2) * ((n * tau2 / tau1)^(1.0 / n)) )
    tt1n = (t / tau1)^n
    return imax / xi * tt1n / (1 + tt1n) * exp(-t / tau2)
end

## GPR
begin
    df = DataFrame(CSV.File("sunjerga173powsys_gpr_100.csv"))
    t = df.t
    nt = length(t)
    tk_max = nt - div(nt, 8)
    inj1 = heidler.(t, 28e3, 1.8e-6, 95e-6, 2)
    inj2 = heidler.(t, 10.7e3, 0.25e-6, 2.5e-6, 2) + heidler.(t, 6.5e3, 2e-6, 230e-6, 2)
    p = plot(t * 1e6, inj1 * 1e-3, label="First stroke",
             xlabel="time [μs]", ylabel="Current [kA]", ylims=(0, 36));
    plot!(t * 1e6, inj2 * 1e-3, label="Subsequent stroke")
    png(p, "sunjerga173powsys_inj")
    for rho in [100, 1000]
        if rho != 100
            df = DataFrame(CSV.File("sunjerga173powsys_gpr_$(rho).csv"))
            t = df.t
        end
        title = @sprintf("σ = %.1f [mS]", 1/rho * 1e3)
        p = plot(t[1:tk_max] * 1e6, df.i1[1:tk_max] * 1e-3, label="First stroke",
                 xlabel="time [μs]", ylabel="GPR [kV]", title=title);
        plot!(t[1:tk_max] * 1e6, df.i2[1:tk_max] * 1e-3, label="Subsequent stroke")
        png(p, "sunjerga173powsys_gpr_$(rho)")
    end
end

## Electric Field
time_steps = [8, 81, 176]
for rho in [100, 1000]
    df = DataFrame(CSV.File("sunjerga173powsys_efield_$(rho).csv"))
    t = sort(unique(df.t))
    nt = length(t)
    ninj = (size(df)[2] - 3) ÷ 2  # number of injections
    strike = ["First strike", "Subsequent strike"]
    for i = 1:ninj
        for k in time_steps
            title = @sprintf("%s, t = %.2f [μs], σ = %.1f [mS]",
                             strike[i], t[k] * 1e6, 1/rho * 1e3)
            c = isapprox.(df.t, t[k])
            # do some Metaprogramming wizardry
            efc = eval(Meta.parse("df.ecx$(i)"))
            efnonc = eval(Meta.parse("df.encx$(i)"))
            m = 1e-3
            p = plot(df.x[c], efc[c] * m, title=title, ylabel="Ex  [kV/m]",
                     xlabel="x [m]", label="E conservative");
            plot!(df.x[c], efnonc[c] * m, label="E induced");
            plot!(df.x[c], (efc + efnonc)[c] * m, label="E total")
            png(p, "sunjerga173powsys_efield$(i)_$(rho)_t$(k)")
        end
    end
end

## Step Voltage along the line y = 0, +x direction
for rho in [100, 1000]
    df = DataFrame(CSV.File("sunjerga173powsys_gpd_$(rho).csv"))
    t = sort(unique(df.t))
    nt = length(t)
    m = 1e-3
    ninj = (size(df)[2] - 3)  # number of injections
    strike = ["First strike", "Subsequent strike"]
    for i = 1:ninj
        for k in time_steps
            title = @sprintf("%s, t = %g [μs], σ = %.1f [mS]",
                             strike[i], t[k] * 1e6, 1/rho * 1e3)
            # GPD difference
            df = DataFrame(CSV.File("sunjerga173powsys_gpd_$(rho).csv"))
            x = sort(unique(df.x))
            dx = x[2] - x[1]
            c = isapprox.(t[k], df.t)
            N = Int(1 ÷ dx)
            nx = length(x) - N
            v = df.i1[c]
            dv = [(v[i] - v[i + N]) for i = 1:nx]
            x = df.x[c][1:end-N]
            p = plot(x, (dv) * m, title=title, ylabel="Step Voltage [kV]",
                     xlabel="x [m]", label="pot. dif.");

            # integral of the E_tot field ========================================
            df = DataFrame(CSV.File("sunjerga173powsys_efield_$(rho).csv"))
            x = sort(unique(df.x))
            dx = x[2] - x[1]
            c = isapprox.(df.t, t[k])
            # Metaprogramming sorcery
            ecx = eval(Meta.parse("df.ecx$(i)"))
            encx = eval(Meta.parse("df.encx$(i)"))
            ef = ecx[c] + encx[c]
            N = Int(1 ÷ dx)
            nx = length(x) - N
            dv = zeros(nx)
            for i = 1:nx
                s = sum(ef[(i + 1):(i + N - 1)])
                dv[i] = dx * (s + (ef[i] + ef[i + N]) / 2)
            end
            x = x[1:end-N]
            plot!(x, (dv) * m, label="voltage")
            png(p, "sunjerga173powsys_stepv$(i)_$(rho)_t$(k)")
        end
    end
end
