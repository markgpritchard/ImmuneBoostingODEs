
using DrWatson
@quickactivate :ImmuneBoostingODEs

using CairoMakie, DataFrames, DifferentialEquations, Pigeons, Random, Turing

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load the data 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

include("rsvsetup.jl")

println("In the 12 months from 1 April each year")
for y ∈ 2016:2022 
    inds = findall(x -> y <= x < y + 1, data.AprilYear)
    println("    $(sum(data.Cases[inds])) cases in $y")
end 


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load results 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rsvparameters01 = loadrsvdata(0.1; ids=1:5) 
plotchains(rsvparameters01)
plotvals01 = fittedsimulationquantiles(rsvparameters01, 0.1, saveat, cbs)

rsvparameters02 = loadrsvdata(0.2; ids=1:5)
plotchains(rsvparameters02)
plotvals02 = fittedsimulationquantiles(rsvparameters02, 0.2, saveat, cbs)

rsvparameters04 = loadrsvdata(0.4; ids=1:5)
plotchains(rsvparameters04)
plotvals04 = fittedsimulationquantiles(rsvparameters04, 0.4, saveat, cbs)

rsvparameters1 = loadrsvdata(1.0; ids=1:5)
plotchains(rsvparameters1)
plotvals1 = fittedsimulationquantiles(rsvparameters1, 1.0, saveat, cbs)

rsvparameters2 = loadrsvdata(2.0; ids=1:5)
plotchains(rsvparameters2)
plotvals2 = fittedsimulationquantiles(rsvparameters2, 2.0, saveat, cbs)

rsvparameters4 = loadrsvdata(4.0; ids=1:5)
plotchains(rsvparameters4)
plotvals4 = fittedsimulationquantiles(rsvparameters4, 4.0, saveat, cbs)

rsvparameters10 = loadrsvdata(10.0; ids=1:5)
plotchains(rsvparameters10)
plotvals10 = fittedsimulationquantiles(rsvparameters10, 10.0, saveat, cbs)

fittedparametersfig = let
    plotvvector = [ 
        plotvals01, plotvals02, plotvals04, plotvals1, plotvals2, plotvals4, plotvals10 
    ]
    pv = [  # parametervector
        rsvparameters01, rsvparameters02, rsvparameters04, 
        rsvparameters1, rsvparameters2, rsvparameters4, rsvparameters10 
    ]
    γ = 48.7
    μ = 0.0087
    logomegavalues = log.([ 0.1, 0.2, 0.4, 1.0, 2.0, 4.0, 10.0 ])
    omegalabels = [ "0.1", "0.2", "0.4", "1.0", "2.0", "4.0", "10.0" ]
    
    fig = Figure(; size=( 800, 1200 ))
    ga = GridLayout(fig[1, 1])
    axs = [ Axis(ga[i, 1]) for i ∈ 1:7 ]
    for (i, v) ∈ enumerate(plotvvector)
        plotfittedsimulationquantiles!(axs[i], data, v, saveat)
        text!(axs[i], 2016.8, 550; text="ω=$(omegalabels[i])", fontsize=11.84)
    end
    gb = GridLayout(fig[1, 2])
    ax1 = Axis(gb[1, 1]; xticks=( logomegavalues, omegalabels ))
    scatter!(
        ax1, logomegavalues, [ quantile(v.log_density, 0.5) for v ∈ pv ]; 
        color=:blue
    )
    rangebars!(
        ax1, 
        logomegavalues, 
        [ quantile(v.log_density, 0.05) for v ∈ pv ], 
        [ quantile(v.log_density, 0.95) for v ∈ pv ];
        color=:blue,
    )
    ax2 = Axis(gb[2, 1]; xticks=( logomegavalues, omegalabels ))
    scatter!(
        ax2, 
        logomegavalues, 
        [ quantile(v.β0, 0.5) for v ∈ pv ] ./ (γ + μ); 
        color=:blue
    )
    rangebars!(
        ax2, 
        logomegavalues, 
        [ quantile(v.β0, 0.05) for v ∈ pv ] ./ (γ + μ), 
        [ quantile(v.β0, 0.95) for v ∈ pv ] ./ (γ + μ);
        color=:blue,
    )
    ax3 = Axis(gb[3, 1]; xticks=( logomegavalues, omegalabels ))
    scatter!(
        ax3, 
        logomegavalues, 
        [ quantile(v.β1, 0.5) for v ∈ pv ] .* [ quantile(v.β0, 0.5) for v ∈ pv ] ./ (γ + μ); 
        color=:blue
    )
    rangebars!(
        ax3, 
        logomegavalues, 
        [ quantile(v.β1, 0.05) for v ∈ pv ] .* [ quantile(v.β0, 0.05) for v ∈ pv ] ./ (γ + μ), 
        [ quantile(v.β1, 0.95) for v ∈ pv ] .* [ quantile(v.β0, 0.95) for v ∈ pv ] ./ (γ + μ);
        color=:blue,
    )
    ax4 = Axis(gb[4, 1]; xticks=( logomegavalues, omegalabels ))
    scatter!(
        ax4, logomegavalues, [ quantile(v.ψ, 0.5) for v ∈ pv ]; 
        color=:blue
    )
    rangebars!(
        ax4, 
        logomegavalues, 
        [ quantile(v.ψ, 0.05) for v ∈ pv ], 
        [ quantile(v.ψ, 0.95) for v ∈ pv ];
        color=:blue,
    )
    ax5 = Axis(gb[5, 1]; xticks=( logomegavalues, omegalabels ))
    scatter!(
        ax5, logomegavalues, 1 .- [ quantile(v.βreduction1, 0.5) for v ∈ pv ]; 
        color=:blue
    )
    rangebars!(
        ax5, 
        logomegavalues, 
        1 .- [ quantile(v.βreduction1, 0.05) for v ∈ pv ], 
        1 .- [ quantile(v.βreduction1, 0.95) for v ∈ pv ];
        color=:blue,
    )
    ax6 = Axis(gb[6, 1]; xticks=( logomegavalues, omegalabels ))
    scatter!(
        ax6, logomegavalues, [ quantile(v.detection, 0.5) for v ∈ pv ] .* 100; 
        color=:blue
    )
    rangebars!(
        ax6, 
        logomegavalues, 
        [ quantile(v.detection, 0.05) for v ∈ pv ] .* 100, 
        [ quantile(v.detection, 0.95) for v ∈ pv ] .* 100;
        color=:blue,
    )
    linkaxes!(axs...)
    for i ∈ 1:7 
        formataxis!(axs[i], hidex=(i != 7), hidexticks=(i != 7))
        if i != 7 hidespines!(axs[i], :b) end
    end
    Label(ga[1:7, 0], "Weekly incidence"; fontsize=11.84, rotation=π/2, tellheight=false)
    Label(ga[8, 1], "Year"; fontsize=11.84, tellwidth=false)
    colgap!(ga, 1, 5)
    rowgap!(ga, 7, 5)
    for (i, ax) ∈ enumerate([ ax1, ax2, ax3, ax4, ax5, ax6 ])
        formataxis!(ax, hidex=(i != 6), hidexticks=(i != 6))
        if i != 6 hidespines!(ax, :b) end
        if i >= 2
            setvalue!(ax, 1, 0)
        end
    end
    Label(gb[1, 0], "Log likelihood"; fontsize=11.84, rotation=π/2, tellheight=false)
    Label(gb[2, 0], "R₀"; fontsize=11.84, rotation=π/2, tellheight=false)
    Label(
        gb[3, 0], "Magnitude of seasonal forcing"; 
        fontsize=11.84, rotation=π/2, tellheight=false
    )
    Label(
        gb[4, 0], "Magnitude of natural\nimmune boosting, ψ"; 
        fontsize=11.84, rotation=π/2, tellheight=false
    )
    Label(
        gb[5, 0], "Transmission reduction from\nnon-pharmaceutical interventions"; 
        fontsize=11.84, rotation=π/2, tellheight=false
    )
    Label(gb[6, 0], "Proportion diagnosed, %"; fontsize=11.84, rotation=π/2, tellheight=false)
    Label(gb[7, 1], "Waning rate, ω"; fontsize=11.84, tellwidth=false)
    colgap!(gb, 1, 5)
    rowgap!(gb, 6, 5)

    fig
end
