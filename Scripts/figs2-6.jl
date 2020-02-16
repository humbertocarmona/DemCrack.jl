#---
using Plots;
using Plots.PlotMeasures
pyplot();
using CSV, DataFrames
using HDF5, JLD
using LsqFit
using Statistics
using LaTeXStrings
PyPlot.rc("text", usetex = "true")
PyPlot.rc("font", family = "CMU Serif")

#---
function τ(h; σ = 3, δ = 1.4, L = 500)
    tau = (1.0 .+ (σ^2.0) * (4 * δ / L)^(2.0 * h))^0.5 # tortuosity
end

#---
function g(filename; hurst = collect(0.3:0.1:0.9))
    df = CSV.read(filename)
    shapes = [
        :xcross,
        :circle,
        :utriangle,
        :square,
        :dtriangle,
        :diamond,
        :cross,
    ]
    ms = [12, 11, 11, 11, 11, 11, 11] .* 0.6
    cl = [1, :orange, :green, :red, :purple, :brown, :lightgray]

    p1 = plot(
        size = (400, 300),
        legendfontsize = 9,
        legend = :inside,
        top_margin = 3mm,
        bottom_margin = 3mm,
        left_margin=3mm,
        right_margin=3mm,
        fg_legend = :white,
        dpi = 300,
        framestyle = :box,
        #xaxis = (L"Re", [0.4, 700],  [1, 10, 100, 1000], font(15), :log),
        #yaxis = (L"G", [11.35, 1000], [10, 100, 1000], font(15), :log),
        xaxis = (L"\log(\mathrm{Re})", [-0.37, 2.8],  [0, 1, 2, 3], font(15)),
        yaxis = (L"\log(G)", [1.05, 2.95], [ 1, 2, 3], font(15)),
        grid = false,
        # thickness_scaling=1.1
    )

    @. model(x, p) = p[1] + p[2] * x + p[3] * x^2 # model to fit G

    coefs = DataFrame(
            h = Float64[],
            α = Float64[],
            β = Float64[],
            γ = Float64[],
        )
    for (i, h) in enumerate(hurst)
        x = df.Re
        gcol = "G0$(Int64(10*h))"
        gcol = df[!, Symbol(gcol)]
        y = gcol

        # ν = 10
        # w = 40.0
        # L = 550.0
        # Le = L*τ(h; η=10, δ=3)
        # dpcol = "DP0$(Int64(10*h))"
        # dpcol = df[!,Symbol(dpcol)]

        # vcol  = "V0$(Int64(10*h))"
        # vcol  = df[!,Symbol(vcol)]
        # y = (w^2).*dpcol./((L*ν).*vcol)


        y0 = mean(y[1:5])
        lb = [y0 - 0.01, -Inf, -Inf]
        ub = [y0 + 0.01, Inf, Inf]
        p0 = [y0, 1.2, 1.2]
        p = curve_fit(model, x, y, p0, lower = lb, upper = ub)
        push!(coefs, [h, p.param[1], p.param[2], p.param[3]])
        xfit = range(minimum(x), stop = maximum(x), length = 200)
        yfit = model(xfit, p.param)
        # p1 = plot!(log10.(xfit), log10.(yfit), linewidth=1, color=cl[i], label="")

        # p1 = scatter!(log10.(x), log10.(y),
        #           markershape=shapes[i],
        #           markersize=ms[i],
        #           markerstrokewidth =0.7,
        #           markerstrokealpha = 1,
        #           markerstrokecolor = :black,
        #           markeralpha=0.9,
        #           color=cl[i],
        #           label="H=$h")

        p1 = plot!(log10.(xfit), log10.(yfit), linewidth = 1, color = cl[i], label = "")

        p1 = scatter!(
            log10.(x),
            log10.(y),
            markershape = shapes[i],
            markersize = ms[i],
            markerstrokewidth = 0.7,
            markerstrokealpha = 1,
            markerstrokecolor = :black,
            markeralpha = 0.9,
            color = cl[i],
            label = L"H="*"$h",
            )
    end

    # fig = p1.o
    # ax = fig[:get_axes]()[1]
    # ax[:spines]["left"][:set_linewidth](20)
    # # PyPlot.draw()
    # ax[:tick_params](axis="x", pad=800)

    return p1, coefs
end

#---
# p1, coefs = g("data/G-all_1.csv")
# savefig(p1, "../Paper/fig2.pdf")
# plot(p1)

#---
function collapse(filename;
                  hurst = collect(0.3:0.1:0.9),
                  plt = true,
                  δs = 0.15:0.05:0.25,
                  )
    """
        this uses the analytical expression for the tortuosity
        fixing the standard deviation η, and using the spacing δ
        as a fitting parameter, using least square...
    """

    shapes = [
        :xcross,
        :circle,
        :utriangle,
        :square,
        :dtriangle,
        :diamond,
        :cross,
    ]
    ms = [11, 10, 10, 10, 10, 10, 10] .* 0.7
    cl = [1, :orange, :green, :red, :purple, :brown, :lightgray]
    p1 = plot(
        size = (400, 300),
        legendfontsize = 10,
        legend = :inside,
        top_margin = 3mm,
        bottom_margin = 3mm,
        left_margin=3mm,
        right_margin=3mm,
        fg_legend = :white,
        dpi = 150,
        frame = :box,
        xaxis = (L"\log(\mathrm{Re}/H)",[0, 3.4],font(15)),
        yaxis = (L"\log(G/\alpha w^2)", [0, 1], [0, 0.5, 1], font(15)),
        grid = false,
    )


    df = CSV.read(filename)
    @. model(x, p) = p[1] + p[2] * x + p[3] * x^2 # model to fit G
    residuals = []
    δs = collect(δs)

    for δ in δs
        Lp = τ.(hurst; σ = 3, δ = 1)
        xmin = 100000000
        xmax = -1000
        for (i, h) in enumerate(hurst)
            x = df.Re .* Lp[i]
            xmin = min(xmin, minimum(x))
            xmax = max(xmax, maximum(x))
        end
        xfit = exp10.(range(log10(xmin), stop = log10(xmax), length = 50))

        # compute Gfit = collapsed G
        Gfit = [[] for _ in hurst]
        for (i, h) in enumerate(hurst)
            x = df.Re / h
            y = df[:, i+1]

            y0 = mean(y[1:5])
            lb = [y0 - 0.01, -Inf, -Inf]
            ub = [y0 + 0.01, Inf, Inf]
            p0 = [y0, 1.2, 1.2]
            p = curve_fit(model, x, y, p0, lower = lb, upper = ub)
            α = p.param[1]
            if plt
                scatter!(
                    log10.(x),
                    log10.(y ./ α),
                    markershape = shapes[i],
                    markersize = ms[i],
                    markerstrokewidth = 0.7,
                    markerstrokealpha = 1,
                    markerstrokecolor = :black,
                    markeralpha = 0.9,
                    color = cl[i],
                    label = "H=$h",
                )
            end
            Gfit[i] = model(xfit, p.param) ./ α
            # p1 = plot!(log10.(xfit), log10.(Gfit[i]), linewidth=3, color=:auto, label="$h")
        end

        # compute resituals
        n = length(Gfit)
        ymean = Gfit[1]
        for i = 2:n
            ymean = ymean + Gfit[i]
        end
        ymean = ymean ./ length(Gfit)

        meanres = 0.0
        for i = 1:n
            res = Gfit[i] - ymean
            meanres += sum(res .^ 2)
        end
        push!(residuals, meanres / n)
    end
    display(δs)
    display(residuals)
    if plt
        return p1
    else
        return δs, residuals
    end
end
#---
#p2 = collapse("data/G-all_1.csv", plt=true, δs=0.025)
#savefig("../Paper/fig4.pdf")
#plot(p2)

#---
function plot_tau(coefs)
    p1 = plot(
        size = (400, 300),
        legendfontsize = 8,
        legend = :inside,
        top_margin = 3mm,
        bottom_margin = 3mm,
        left_margin=3mm,
        right_margin=3mm,
        fg_legend = :lightgray,
        dpi = 150,
        frame = :box,
        xaxis = (L"H", font(11)),
        yaxis = (L"\tau", font(11)),
        grid = false,
        background_color_outside = :transparent,
        background_color_inside = :transparent,
        background_color_legend = :transparent,
    )

    coefs
    df = CSV.read("./data/surf.csv")
    xfit = 0.29:0.01:0.92
    Lp = τ.(xfit, σ = 12, δ = 1.4, L = 510)
    scatter!(coefs.h, df.Lp, label = "num", markersize = 7)
    plot!(xfit, Lp, color = :black, label = "Fit", linestyle = :dash)
    return p1
end
# plot_tau(coefs)

#---
function plotalphaH(coefs)
    """
    simply \alpha x H
    """
    df = CSV.read("./data/surf.csv")
    x = coefs.h
    y = coefs.α
    p4 = plot(
        size = (400, 300),
        legendfontsize = 9,
        legend = :inside,
        top_margin = 3mm,
        bottom_margin = 3mm,
        left_margin=3mm,
        right_margin=3mm,
        fg_legend = :white,
        dpi = 150,
        frame = :box,
        grid = false,
        xaxis = (L"H", font(15), [0.3, 0.9]),
        yaxis = (L"\alpha w^2", [11.65, 68.35], font(15)),
        # background_color_outside = :transparent,
        # background_color_inside = :transparent,
        # background_color_legend = :transparent,
    )
    p4 = plot!(
        x,
        y,
        markershape = :circle,
        linewidth = 0.0,
        markersize = 9,
        label = "",
    )


    plot!(
        [12],
        seriestype = :hline,
        linewidth = 1,
        color = :black,
        linestyle = :dash,
        label = L"\alpha w^2=12",
    )
    xfit = range(0.29, stop = 0.92, length = 200)
    yfit = 58.9 * τ.(xfit, σ = 5, δ = 2, L = 500) .- 42.4
    plot!(xfit, yfit, label = "", color = :black, linewidth = 1)
    return p4
end
#pα = plotalphaH(coefs)
#pα
#savefig("../Paper/fig3_inset.pdf")

#---
function plotalpha(coefs)
    """
    this used "experimental"
    """
    df = CSV.read("./data/surf.csv")
    #x = df.area
    Lp = τ.(coefs.h, σ = 5, δ = 2, L = 500)
    x = Lp
    #     x = df.Lp
    y = coefs.α
    p4 = plot(
        size = (400, 300),
        legendfontsize = 9,
        legend = :inside,
        top_margin = 3mm,
        bottom_margin = 3mm,
        left_margin=3mm,
        right_margin=3mm,
        fg_legend = :white,
        dpi = 300,
        frame = :box,
        grid = false,
        xaxis = (L"\tau", font(15)),
        yaxis = (L"\alpha w^2", font(15)),
        # background_color_outside = :transparent,
        # background_color_inside = :transparent,
        # background_color_legend = :transparent,
    )
    p4 = plot!(
        x,
        y,
        markershape = :circle,
        linewidth = 0.0,
        markersize = 8,
        label = "",
    )
    @. model(x, p) = p[1] .+ p[2] .* x
    fit = curve_fit(model, x[1:5], y[1:5], [1.0, 1.0])

    yav = mean(y)
    ss_tot = sum((y .- yav) .^ 2)
    fi = model(x, fit.param)
    ss_res = sum((fi .- yav) .^ 2)

    R2 = ss_res / ss_tot
    println(ss_tot, " ", ss_res, " ", R2)

    xfit = 0.97:0.1:1.9
    yfit = model(xfit, fit.param)
    plot!(xfit, yfit, color = :black, label = "")
    return p4, fit.param
end
#---
#pα, params = plotalpha(coefs)
#pα
#savefig("../Paper/fig3.pdf")
# savefig("../../paper18122019/fig3.png")

#---
function plotbeta(coefs)
    y = coefs.β
    x = coefs.α ./ coefs.h

    @. model(x, p) = p[1] * x + p[2]
    p4 = plot(
        size = (400, 300),
        legendfontsize = 9,
        legend = :inside,
        top_margin = 3mm,
        bottom_margin = 3mm,
        left_margin=3mm,
        right_margin=3mm,
        fg_legend = :white,
        dpi = 150,
        frame = :box,
        grid = false,
        xaxis = (L"\alpha w^2/H", font(15)),
        yaxis = (L"\beta w", font(15)),
        # background_color_outside = :transparent,
        # background_color_inside = :transparent,
        # background_color_legend = :transparent,
    )
    p4 = scatter!(
        x,
        y,
        markershape = :circle,
        markersize = 8,
        linewidth = 0.5,
        label = "",
    )
    p0 = [0.02, 0.0]
    lb = [0.0, -0.1]
    ub = [0.1, 0.0]

    #     p0 = [0.02]
    #     lb = [0.0]
    #     ub = [0.1]

    fit = curve_fit(model, x, y, p0, lower = lb, upper = ub)
    xfit = 0:210
    yfit = model(xfit, fit.param)

    yav = mean(y)
    ss_tot = sum((y .- yav) .^ 2)
    fi = model(x, fit.param)
    ss_res = sum((fi .- yav) .^ 2)

    R2 = ss_res / ss_tot
    println(ss_tot, " ", ss_res, " ", R2)

    println(coef(fit))
    println(stderror(fit))
    plot!(xfit, yfit, color = :black, label = "")
    return p4, fit.param
end

#---
# pβ, params = plotbeta(coefs)
# pβ
# savefig("../Paper/fig5.pdf")

#---
function plotgamma(coefs)

    x = coefs.α ./ coefs.h .^ 2
    y = coefs.γ
    p4 = plot(
        size = (400, 300),
        legendfontsize = 9,
        legend = :inside,
        top_margin = 6mm,
        bottom_margin = 1mm,
        left_margin=3mm,
        right_margin=3mm,
        fg_legend = :white,
        dpi = 150,
        frame = :box,
        grid = false,
        xaxis = (L"\alpha w^2/H^2", [20,700],font(15)),
        yaxis = (L"\gamma", [0, 2.3e-4], font(15)),
        yticks = ([0,1e-4,2e-4, 2.5e-4], [L"0", L"1", L"2"])
                  #background_color_outside=:transparent,
                  #background_color_inside=:transparent,
    )

    annotate!([(25, 2.57e-4, text(L"\times 10^{-4}", font(9)))])

    p4 = plot!(
        x,
        y,
        markershape = :utriangle,
        markersize = 9,
        color = 2,
        linewidth = 0.8,
        label = "",
    )
    @. model(x, p) = p[1] * x + p[2]
    fit = curve_fit(
        model,
        x[end-4:end],
        y[end-4:end],
        [1.0, 0.0],
        lower = [0.0, -0.1],
        upper = [Inf, 0.0],
    )
    xfit = 0:240
    yfit = model(xfit, fit.param)
    plot!(xfit, yfit, lw = 1.2, color = :black, label = "")


    yav = mean(y[end-4:end])
    ss_tot = sum((y[end-4:end] .- yav) .^ 2)
    fi = model(x[end-4:end], fit.param)
    ss_res = sum((fi .- yav) .^ 2)

    R2 = ss_res / ss_tot
    println(ss_tot, " ", ss_res, " ", R2)

    println(stderror(fit))
    println(coef(fit))
    println(x)
    return p4, fit.param
end
#---
#pγ, params = plotgamma(coefs)
#pγ
#savefig("../Paper/fig5b.pdf")
# print(params)

#---

function fig6(fname="data/pi_all.csv")
    df = CSV.read(fname)

    # filter Re

    df3d = filter(row->(0<row[:v3d]<22),df)
    df3dm = aggregate(df3d, :H, mean)


    p = plot(size = (400, 300),
            legendfontsize = 11,
            legend = :bottomright,
            top_margin = 3mm,
            bottom_margin = 3mm,
            left_margin=3mm,
            right_margin=3mm,
            fg_legend = :white,
            background_color_legend = :transparent,
            dpi = 150,
            framestyle = :box,grid = false,
            xaxis=(L"H", [0.2,0.9], font(15)),
            yaxis=(L"\pi", [0.4,0.99], font(15)))

    p = scatter!(df3dm[:H], df3dm[:pi3d_mean], marker=:circle, markersize=10, color=2, label="Bulk")
    p = scatter!(df3dm[:H], df3dm[:pi2d_mean], marker=:utriangle, markersize=10,
    color=1, label=L"w/2"*" level")
    plot!(
        [0.7],
        seriestype = :hline,
        linewidth = 0.7,
        color = :black,
        linestyle = :dash,
        label = "")
    return p
    end
p=fig6()
savefig("../Paper/fig6.pdf")
