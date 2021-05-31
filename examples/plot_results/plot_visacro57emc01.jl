# Plot the results from example visacro57emc01
using LinearAlgebra
using DataFrames
using CSV
import Printf.@sprintf
using Plots
pyplot()

nomes = ["corner fast", "centro fast", "corner slow", "centro slow"];
time_steps = [8, 81, 176];

## GPR and injection
begin
    #1  ======================================================================
    df = DataFrame(CSV.File("visacro57emc01_gpr.csv"))
    t = df.t
    nt = length(t)
    p = plot(t * 1e6, df.i1, label="simulated", linestyle=:solid, linecolor=:blue,
             xlabel="time [μs]", ylabel="GPR [V]", ylims=(-1, 25), xlims=(0,25))
    df2 = DataFrame(CSV.File("examples/visacro57emc01_auxfiles/vt_corner_fast.csv", header=false))
    plot!(t * 1e6, df2[:,1][1:nt], label="measured")
    png(p, "visacro57emc01_gpr_corner_fast")
    #2  ======================================================================
    p = plot(t * 1e6, df.i2, label="simulated", linestyle=:solid, linecolor=:blue,
             xlabel="time [μs]", ylabel="GPR [V]", ylims=(-1, 25), xlims=(0,25))
    df2 = DataFrame(CSV.File("examples/visacro57emc01_auxfiles/vt_center_fast.csv", header=false))
    plot!(t * 1e6, df2[:,1][1:nt], label="measured")
    png(p, "visacro57emc01_gpr_center_fast")
    #3  ======================================================================
    p = plot(t * 1e6, df.i3, label="simulated", linestyle=:solid, linecolor=:blue,
             xlabel="time [μs]", ylabel="GPR [V]", ylims=(-1, 25), xlims=(0,25))
    df2 = DataFrame(CSV.File("examples/visacro57emc01_auxfiles/vt_corner_slow.csv", header=false))
    plot!(t * 1e6, df2[:,1][1:nt], label="measured")
    png(p, "visacro57emc01_gpr_corner_slow")
    #4  ======================================================================
    p = plot(t * 1e6, df.i4, label="simulated", linestyle=:solid, linecolor=:blue,
             xlabel="time [μs]", ylabel="GPR [V]", ylims=(-1, 25), xlims=(0,25))
    df2 = DataFrame(CSV.File("examples/visacro57emc01_auxfiles/vt_center_slow.csv", header=false))
    plot!(t * 1e6, df2[:,1][1:nt], label="measured")
    png(p, "visacro57emc01_gpr_center_slow")
end

## Electric Field
begin
    df = DataFrame(CSV.File("visacro57emc01_efield.csv"))
    t = sort(unique(df.t))
    nt = length(t)
    ninj = (size(df)[2] - 3) ÷ 2  # number of injections
    for i = 1:ninj
        # do some Metaprogramming wizardry
        ecx = eval(Meta.parse("df.ecx$(i)"))
        encx = eval(Meta.parse("df.encx$(i)"))
        for yi in [0, 8]
            for k in time_steps
                c = isapprox.(df.t, t[k]) .& isapprox.(df.y, yi)
                efc = ecx[c]
                efnonc = encx[c]
                x = df.x[c]
                p = plot(x, efc, ylabel="Ex  [V/m]", xlabel="x [m]",
                         label="E conservative", title="t = $(round(1.0e6 * t[k], sigdigits=2)) [μs]");
                plot!(x, efnonc, label="E induced");
                plot!(x, (efc + efnonc), label="E total")
                png(p, "visacro57emc01_efield$(i)_y$(yi)_t$(k)")
            end
        end
    end
end


## Step Voltage along a line in the +x direction
begin
    df_efield = DataFrame(CSV.File("visacro57emc01_efield.csv"))
    t = sort(unique(df_efield.t))
    nt = length(t)
    ninj = (size(df)[2] - 3) ÷ 2  # number of injections
    for i = 1:ninj
        ecx = eval(Meta.parse("df_efield.ecx$(i)"))
        encx = eval(Meta.parse("df_efield.ecx$(i)"))
        for yi in [0, 8]
            for k in time_steps
                # integral of the E conservative field =======================
                x = sort(unique(df_efield.x))
                dx = x[2] - x[1]
                c = isapprox.(df_efield.t, t[k]) .& isapprox.(df_efield.y, yi)
                ef = ecx[c] #+ encx[c]
                N = Int(1 ÷ dx)
                nx = length(x) - N
                dv = zeros(nx)
                for i = 1:nx
                    s = sum(ef[(i + 1):(i + N - 1)])
                    dv[i] = dx * (s + (ef[i] + ef[i + N]) / 2)
                end
                x = x[1:end-N]
                p = plot(x, (dv), ylabel="Step Voltage [V]",
                         xlabel="x [m]", label="pot. dif.");

                # integral of the E total field ==============================
                x = sort(unique(df_efield.x))
                dx = x[2] - x[1]
                c = isapprox.(df_efield.t, t[k]) .& isapprox.(df_efield.y, yi)
                ef = ecx[c] + encx[c]
                N = Int(1 ÷ dx)
                nx = length(x) - N
                dv = zeros(nx)
                for i = 1:nx
                    s = sum(ef[(i + 1):(i + N - 1)])
                    dv[i] = dx * (s + (ef[i] + ef[i + N]) / 2)
                end
                x = x[1:end-N]
                plot!(x, (dv), label="Voltage")
                png(p, "visacro57emc01_stepv$(i)_y$(yi)_t$(k)")
            end
        end
    end
end

## GPD profile
import PlotlyJS
begin
    df = DataFrame(CSV.File("visacro57emc01_gpd.csv"))
    t = sort(unique(df.t))
    nt = length(t)
    ninj = (size(df)[2] - 3)  # number of injections
    for i = 1:ninj
        gpd = eval(Meta.parse("df.i$(i)"))
        for k in time_steps
            c = isapprox.(df.t, t[k])
            data = PlotlyJS.contour(; x=df.x[c], y=df.y[c], z=gpd[c],
                                    colorbar=PlotlyJS.attr(;title="Electric Potential [V]",
                                                  titleside="right"),
                                                  titlefont=PlotlyJS.attr(;size=16,
                                                  family="Arial, sans-serif"))
            p = PlotlyJS.plot(data)
            PlotlyJS.savefig(p, "visacro57emc01_gpd$(i)_t$(k).png")
        end
    end
end
