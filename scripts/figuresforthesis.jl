
# Runs code from other scripts and produces figures for manuscript 

using DrWatson
@quickactivate :ImmuneBoostingODEs
using CairoMakie, MakieTeX, PlotFormatting
CairoMakie.activate!()  # allows figures to be saved as vector files

include("immuneduration.jl") # provides immunedurationmodel
include("equilibriumvalues.jl") # provides equilSs, equilIs, equilRs, equilparms, 
    # critvector, bifurcationI_1_5, bifurcationI_5, bifurcationI_10, bifurcationI_15, 
    # basic5, basic1_5, basic5_, stablesim, unstablein12, unstableout12, unstablein5, 
    # unstableout5, limitcycle12, limitcycle5
include("fouriertransforms.jl") # provides simparms, unforced6mfreqs, unforced400dfreqs, 
    # unforced25freqs, unforced6mdensities, unforced400ddensities, unforced25densities, 
    # forced6mfreqs, forced400dfreqs, forced25freqs
include("npisimulation.jl") # provides crgtdata, reduceday, increaseday, npisim_phi0, 
    # npisim_phi5, npisim_phi13_2, npiparms
include("displayrsvanalysis.jl") # provides data, rsvsim_psi0, rsvsim_psi5, rsvsim_psi13_2,
    # agedata, rsvparms

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Model flowchart 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

introplot = let 
    fig = Figure(; size=( 500, 100 ))
    lt1 = LTeX(
        fig[1, 1], 
        TeXDocument(read(scriptsdir("tikz_sirrrs.tex"), String)); 
        tellwidth=false
    )

    fig
end

# save the figure
safesave(plotsdir("introplot.pdf"), introplot)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Proportions in each compartment at equilibrium 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

equilibriumproportionsplot = with_theme(theme_latexfonts()) do
    @unpack γ, μ, psis, ω, R0s = equilparms
    betas = R0s .* (γ + μ)

    values = [ equilSs, equilIs, equilRs ]
    labels = [ "Susceptible", "Infectious", "Resistant (immune)" ]

    fig = Figure(; size = ( 500, 200 ))
    axs = [ Axis(fig[2, i]) for i ∈ 1:3 ]
    for (i, val) ∈ enumerate(values) 
        if i == 2 
            cp = contourf!(axs[i], R0s, psis, val; levels=0.0:0.0016:0.016)
            cb = Colorbar(fig[1, i], cp; ticks=[ 0, 0.01 ], vertical=false)
        else 
            cp = contourf!(axs[i], R0s, psis, val; levels=0.0:0.1:1.0)
            cb = Colorbar(fig[1, i], cp; vertical=false)
        end
        formataxis!(cb; horizontal=true)
        #cb.topspinevisible=false
        #cb.topspinecolor=( :white, 0 )
        cb.flipaxis = false
    end

    for (i, lbl) ∈ enumerate(labels)
        Label(fig[0, i], lbl; halign=:left, fontsize=11.84, tellwidth=false)
        formataxis!(
            axs[i]; 
            hidey=(i != 1), hidespines=( :r, :t ), 
            setorigin=true, setpoint=( 15, 0 ),
            trimspines=true,
        )
    end

    Label(
        fig[3, 1:3], L"Basic reproduction ratio, $\mathcal{R}_0$"; 
        fontsize=11.84, tellwidth=false
    )
    Label(
        fig[2, 0], L"Boosting parameter, $\psi$"; 
        fontsize=11.84, rotation=π/2, tellheight=false
    )
    
    colgap!(fig.layout, 1, 5)
    for r ∈ [ 1, 3 ] rowgap!(fig.layout, r, 5) end
    rowgap!(fig.layout, 2, 8)

    fig
end

safesave(plotsdir("equilibriumproportionsplot.pdf"), equilibriumproportionsplot)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Critical levels of psi 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

critpsiplot = with_theme(theme_latexfonts()) do
    @unpack γ, μ, durations = equilparms
    R0s = 0:0.1:17.5 

    md = [ 
        [ rand(truncated(Normal(65, 3), 0, 80)), rand(truncated(Normal(15.25, 0.5), 1, 17.5)) ] 
        for _ ∈ 1:100_000 
    ]
    mp = [ 
        [ rand(truncated(Normal(30, 3), 0, 80)), rand(truncated(Normal(14.65, 0.5), 1, 17.5)) ] 
        for _ ∈ 1:100_000 
    ]
    mv = [ 
        [ rand(truncated(Normal(50, 3), 0, 80)), rand(truncated(Normal(8.75, 0.5), 1, 17.5)) ] 
        for _ ∈ 1:100_000 
    ]
    mr = [ 
        [ rand(truncated(Normal(2.6, 2), 0, 80)), rand(truncated(Normal(3.5, 1), 1, 17.5)) ] 
        for _ ∈ 1:100_000 
    ]

    fig = Figure(; size = ( 500, 350 ))
    ax = Axis(fig[1, 1]; yticks=[ 0, 1, 5, 10, 15 ])
    cp = contourf!(ax, durations, R0s, critpsi'; levels=0:01:10, extendhigh=:auto)
    cb = Colorbar(fig[1, 2], cp)
    scatter!(ax, Point2f.(md); markersize=1, color=( :red, 0.2 ))
    scatter!(ax, Point2f.(mp); markersize=1, color=( :red, 0.2 ))
    scatter!(ax, Point2f.(mv); markersize=1, color=( :red, 0.2 ))
    scatter!(ax, Point2f.(mr); markersize=1, color=( :red, 0.2 ))
    text!(ax, 65, 15.25; text="Measles", align=( :center, :center ), fontsize=10)
    text!(ax, 30, 14.65; text=L"$$\textit{B. pertussis}", align=( :center, :center ), fontsize=10)
    text!(ax, 50, 8.75; text="Varicella\nzoster", align=( :center, :center ), fontsize=10)
    text!(ax, 2.6, 3.5; text="Respiratory\nviruses", align=( :center, :center ), fontsize=10, rotation=3π/2)
    hlines!(ax, 1; color=RGBAf(0, 0, 0, 0.12), linestyle=( :dot, :dense ), linewidth =1,)

    formataxis!(ax; setpoint=( 0, 1 ), trimspines=true, hidespines=( :r, :t ))
    formataxis!(cb)

    Label(
        fig[2, 1], "Mean duration of immunity without boosting, years"; 
        fontsize=11.84, tellwidth=false
    )
    Label(
        fig[1, 0], "Basic reproduction ratio"; 
        fontsize=11.84, rotation=π/2, tellheight=false
    )
    Label(
        fig[1, 3], L"Critical boosting coefficient, $\psi^*$"; 
        fontsize=11.84, rotation=3π/2, tellheight=false
    )
    
    for c ∈ [1 , 3 ] colgap!(fig.layout, c, 5) end
    colgap!(fig.layout, 2, 8)
    rowgap!(fig.layout, 1, 5)

    fig
end

safesave(plotsdir("critpsiplot.pdf"), critpsiplot)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Effect of boosting on equilibria 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

equilibriumplot = with_theme(theme_latexfonts()) do 
    fig = Figure(; size=( 500, 300 ))
    ga = GridLayout(fig[1, 1])
    plotequilibriumc!(ga, basic5, basic1_5, basic5_)
    plotsi!(
        ga, stablesim, 
        unstablein12, unstableout12, unstablein5, unstableout5, 
        limitcycle12, limitcycle5;
        row=5
    )
    rowsize!(ga, 5, Auto(2.25))
    labelplots!([ "A", "B" ], ga; rows=[ 0, 5 ])
    fig
end

safesave(plotsdir("equilibriumplot.pdf"), equilibriumplot)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Fourier transforms
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

fourierplot = with_theme(theme_latexfonts()) do 
    fig = Figure(; size = ( 500, 300 ))
    plotfourier!(
        fig, unforced6mfreqs, unforced400dfreqs, unforced25freqs, 
        unforced6mdensities, unforced400ddensities, unforced25densities, forced6mfreqs, 
        forced400dfreqs, forced25freqs, forced6mdensities, forced400ddensities, 
        forced25densities, simparms
    )
    fig
end

safesave(plotsdir("fourierplot.pdf"), fourierplot)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Non-pharmaceutical interventions plot 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

npisimulationplot = with_theme(theme_latexfonts()) do 
    fig = Figure(; size = ( 500, 390 ))
    plotnpi!(fig, npisim_phi0, npisim_phi5, npisim_phi13_2, npiparms)
    fig
end

safesave(plotsdir("npisimulationplot.pdf"), npisimulationplot)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Respiratory syncytial virus data 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rsvcases = with_theme(theme_latexfonts()) do 
    inds = findall(x -> x >= 50, crgtdata.StringencyIndex_Average)
    reduceday = crgtdata.Date[inds[1]]
    increaseday = crgtdata.Date[last(inds)]

    fig = Figure(; size = ( 500, 200 ))
    ax = Axis(fig[1, 1]; xticks=[ 2017, 2019, 2021, 2023 ])
    lines!(ax, saveat[2:end], data.Cases; color=:black, linewidth=1,)
    vspan!(ax, reduceday, increaseday, color=( :gray, 0.1 ))
    for x ∈ 2017:1:2023
        vlines!(
            ax, x; 
            color=RGBAf(0, 0, 0, 0.12), linestyle=( :dot, :dense ), linewidth=1,
        )
    end
    formataxis!(ax; trimspines=true, hidespines=( :t, :r ))
    Label(fig[1, 0], "Weekly incidence"; fontsize=11.84, rotation=π/2, tellheight=false)
    Label(fig[2, 1], "Year"; fontsize=11.84, tellwidth=false)
    colgap!(fig.layout, 1, 5)
    rowgap!(fig.layout, 1, 5)
    fig
end

safesave(plotsdir("rsvcases.pdf"), rsvcases)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Respiratory syncytial virus data with simulations
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

fittedparametersfig = let 
    fig = Figure(; size=( 800, 1000 ))
    plotvvector = [ 
        plotvals01, plotvals02, plotvals04, plotvals1, plotvals2, plotvals4, plotvals6 
    ]
    parametervector = [
        rsvparameters01, rsvparameters02, rsvparameters04, 
        rsvparameters1, rsvparameters2, rsvparameters4, rsvparameters6 
    ]
    plotfittedsimulations!(fig, plotvvector, parametervector, data, crgtdata, saveat)  
    
    fig
end

safesave(plotsdir("fittedparametersfig.pdf"), fittedparametersfig)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Stringency Index 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

stringencyplot = Figure(; size = ( 400, 400 ))
plotstringency!(stringencyplot, crgtdata, reduceday, increaseday) 

stringencyplot
safesave(plotsdir("stringencyplot.pdf"), stringencyplot)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Non-pharmaceutical interventions plot 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

npisimulationplot = Figure(; size = ( 800, 500 ))
plotnpi!(npisimulationplot, npisim_phi0, npisim_phi5, npisim_phi13_2, npiparms) 

npisimulationplot
safesave(plotsdir("npisimulationplot.pdf"), npisimulationplot)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Age-stratified respiratory syncytial virus data  
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rsvagefigure = with_theme(theme_latexfonts()) do 
    fig = Figure(; size = ( 500, 350 ))
    plotrsvage!(fig, agedata)
    fig
end

safesave(plotsdir("rsvagefigure.pdf"), rsvagefigure)
