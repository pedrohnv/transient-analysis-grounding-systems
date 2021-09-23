using LinearAlgebra
using DataFrames
using CSV
using FFTW
import Printf.@sprintf
using Plots

font1 = Plots.font("DejaVu Sans", 16)
font2 = Plots.font("DejaVu Sans", 16)
pyplot(guidefont=font1, xtickfont=font1, ytickfont=font1, legendfont=font2);
colors = [:blue, :darkorange, :blue, :darkorange]
linestyles = [:solid, :solid, :dash, :dash]

connect = true ? "_connected" : "_isolated"

# Harmonic Analysis
## Zh
begin
    p = plot(xlabel="Frequency [Hz]", ylabel="Harmonic Impedance [Ω]",
             legend=:topleft, background_color_legend=nothing)
    line = 1
    for rho in [100, 1000]
        df = DataFrame(CSV.File("three_grids_freq_gpr_$(rho)$(connect).csv"))
        f = sort( unique(df.f) )
        nf = length(f)
        zh = parse.(ComplexF64, df.i1)
        c = @. isapprox(df.x, 0.0) & isapprox(df.y, 0.0)
        plot!(f, abs.(zh[c]), xaxis=:log, label="(0,0), σ = $(1/rho * 1e3) mS")
        c = @. isapprox(df.x, 18.0) & isapprox(df.y, 0.0)
        plot!(f, abs.(zh[c]), xaxis=:log, label="(18,0), σ = $(1/rho * 1e3) mS")
        c = @. isapprox(df.x, 36.0) & isapprox(df.y, 0.0)
        plot!(f, abs.(zh[c]), xaxis=:log, label="(36,0), σ = $(1/rho * 1e3) mS")
        line += 1
    end
    png(p, "three_grids_freq_zh_abs$(connect)")
    # ========================================
    p = plot(xlabel="Frequency [Hz]", ylabel="Harmonic Impedance Phase [deg]",
             legend=:topleft, background_color_legend=nothing)
    line = 1
    for rho in [100, 1000]
        df = DataFrame(CSV.File("three_grids_freq_gpr_$(rho)$(connect).csv"))
        f = sort(unique(df.f))
        nf = length(f)
        zh = parse.(ComplexF64, df.i1)
        c = @. isapprox(df.x, 0.0) & isapprox(df.y, 0.0)
        plot!(f, rad2deg.(angle.(zh[c])), xaxis=:log, label="(0,0), σ = $(1/rho * 1e3) mS")
        c = @. isapprox(df.x, 18.0) & isapprox(df.y, 0.0)
        plot!(f, rad2deg.(angle.(zh[c])), xaxis=:log, label="(18,0), σ = $(1/rho * 1e3) mS")
        c = @. isapprox(df.x, 36.0) & isapprox(df.y, 0.0)
        plot!(f, rad2deg.(angle.(zh[c])), xaxis=:log, label="(36,0), σ = $(1/rho * 1e3) mS")
        line += 1
    end
    png(p, "three_grids_freq_zh_arg$(connect)")
end

## Step Voltage along the line y = 0, +x direction
colors = [:blue, :blue, :darkorange, :darkorange]
linestyles = [:solid, :dash, :solid, :dash]
for k in [1, 70, 80, 90]
    #title = @sprintf("f = %s", freq[k])
    title = ""
    p = plot(title=title, xlabel="x [m]", ylabel="Step Voltage [V]",
             legend=:topright, background_color_legend=nothing);
    line = 1
    for rho in [100, 1000]
        # integral of the E_C field ========================================
        df = DataFrame(CSV.File("three_grids_freq_efield_$(rho)$(connect).csv"))
        f = sort(unique(df.f))
        nf = length(f)
        x = sort(unique(df.x))
        dx = x[2] - x[1]
        c = isapprox.(df.f, f[k])
        ef = parse.(ComplexF64, df.ecx1[c])
        N = Int(1 ÷ dx)
        nx = length(x) - N
        dv1 = zeros(ComplexF64, nx)
        for i = 1:nx
            s = sum(ef[(i + 1):(i + N - 1)])
            dv1[i] = dx * (s + (ef[i] + ef[i + N]) / 2)
        end
        x1 = x[1:end-N]
        plot!(p,x1, abs.(dv1), label="pot. dif., σ = $(1/rho * 1e3) mS",
              color=colors[line], linestyles=linestyles[line])
        line += 1
        # integral of the E_tot field ========================================
        ef = parse.(ComplexF64, df.ecx1[c]) + parse.(ComplexF64, df.encx1[c])
        dv2 = zeros(ComplexF64, nx)
        for i = 1:nx
            s = sum(ef[(i + 1):(i + N - 1)])
            dv2[i] = dx * (s + (ef[i] + ef[i + N]) / 2)
        end
        x2 = x[1:end-N]
        plot!(p, x2, abs.(dv2), label="voltage, σ = $(1/rho * 1e3) mS",
              color=colors[line], linestyles=linestyles[line])
        line += 1
    end
    png(p, "three_grids_freq_stepv_$(k)$(connect)")
end

# Time domain
## GPR
begin
    p = plot(xlabel="Time [us]", ylabel="GPR [V]",
             legend=:topleft, background_color_legend=nothing)
    line = 1
    for rho in [100, 1000]
        df = DataFrame(CSV.File("three_grids_time_gpr_$(rho)$(connect).csv"))
        t = sort(unique(df.t)) * 1e6
        nt = length(t)
        zh = df.i1
        c = @. isapprox(df.x, 0.0) & isapprox(df.y, 0.0)
        plot!(t, (zh[c]), label="(0,0), σ = $(1/rho * 1e3) mS")
        c = @. isapprox(df.x, 18.0) & isapprox(df.y, 0.0)
        plot!(t, (zh[c]), label="(18,0), σ = $(1/rho * 1e3) mS")
        c = @. isapprox(df.x, 36.0) & isapprox(df.y, 0.0)
        plot!(t, (zh[c]), label="(36,0), σ = $(1/rho * 1e3) mS")
    end
    png(p, "three_grids_time_gpr")
end

## Step Voltage along the line y = 0, +x direction
colors = [:blue, :blue, :darkorange, :darkorange]
linestyles = [:solid, :dash, :solid, :dash]
for k in [1, 70, 80, 90]
    #title = @sprintf("f = %s", freq[k])
    title = ""
    p = plot(title=title, xlabel="x [m]", ylabel="Step Voltage [V]",
             legend=:topright, background_color_legend=nothing);
    line = 1
    for rho in [100, 1000]
        # integral of the E_C field ========================================
        df = DataFrame(CSV.File("three_grids_time_efield_$(rho)$(connect).csv"))
        t = sort(unique(df.t))
        nt = length(t)
        x = sort(unique(df.x))
        dx = x[2] - x[1]
        c = isapprox.(df.t, t[k])
        ef = df.ecx1[c]
        N = Int(1 ÷ dx)
        nx = length(x) - N
        dv1 = zeros(Float64, nx)
        for i = 1:nx
            s = sum(ef[(i + 1):(i + N - 1)])
            dv1[i] = dx * (s + (ef[i] + ef[i + N]) / 2)
        end
        x1 = x[1:end-N]
        plot!(p, x1, dv1, label="pot. dif., σ = $(1/rho * 1e3) mS",
              color=colors[line], linestyles=linestyles[line])
        line += 1
        # integral of the E_tot field ========================================
        ef = df.ecx1[c] + df.encx1[c]
        dv2 = zeros(Float64, nx)
        for i = 1:nx
            s = sum(ef[(i + 1):(i + N - 1)])
            dv2[i] = dx * (s + (ef[i] + ef[i + N]) / 2)
        end
        x2 = x[1:end-N]
        plot!(p, x2, dv2, label="voltage, σ = $(1/rho * 1e3) mS",
              color=colors[line], linestyles=linestyles[line])
        line += 1
    end
    png(p, "three_grids_time_stepv_$k")
end
