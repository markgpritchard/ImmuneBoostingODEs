
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

rsvparameters01 = loadrsvdata(0.1)
plotchains(rsvparameters01)
plotvals01 = fittedsimulationquantiles(rsvparameters01, 0.1, saveat, cbs)

rsvparameters02 = loadrsvdata(0.2)
plotchains(rsvparameters02)
plotvals02 = fittedsimulationquantiles(rsvparameters02, 0.2, saveat, cbs)

rsvparameters04 = loadrsvdata(0.4)
plotchains(rsvparameters04)
plotvals04 = fittedsimulationquantiles(rsvparameters04, 0.4, saveat, cbs)

rsvparameters1 = loadrsvdata(1.0)
plotchains(rsvparameters1)
plotvals1 = fittedsimulationquantiles(rsvparameters1, 1.0, saveat, cbs)

rsvparameters2 = loadrsvdata(2.0)
plotchains(rsvparameters2)
plotvals2 = fittedsimulationquantiles(rsvparameters2, 2.0, saveat, cbs)

rsvparameters4 = loadrsvdata(4.0)
plotchains(rsvparameters4)
plotvals4 = fittedsimulationquantiles(rsvparameters4, 4.0, saveat, cbs)

rsvparameters10 = loadrsvdata(10.0)
plotchains(rsvparameters10)
plotvals10 = fittedsimulationquantiles(rsvparameters10, 10.0, saveat, cbs)

plotvvector = [ 
    plotvals01, plotvals02, plotvals04, plotvals1, plotvals2, plotvals4, plotvals10 
]
parametervector = [ 
    rsvparameters01, rsvparameters02, rsvparameters04, 
    rsvparameters1,rsvparameters2, rsvparameters4, rsvparameters10 
]
logomegavalues = log.([ 0.1, 0.2, 0.4, 1.0, 2.0, 4.0, 10.0 ])
omegalabels = [ "0.1", "0.2", "0.4", "1.0", "2.0", "4.0", "10.0" ]

fittedparametersfig = Figure(; size=( 800, 1200 ))
ga = GridLayout(fittedparametersfig[1, 1])
axs = [ Axis(ga[i, 1]) for i ∈ 1:7]
for (i, v) ∈ enumerate(plotvvector)
    plotfittedsimulationquantiles!(axs[i], data, v, saveat)
end
gb = GridLayout(fittedparametersfig[1, 2])
ax1 = Axis(gb[1, 1]; xticks=( logomegavalues, omegalabels ))
scatter!(ax1, logomegavalues, [ quantile(v.log_density, 0.5) for v ∈ parametervector ]; color=:blue)
rangebars!(
    ax1, 
    logomegavalues, 
    [ quantile(v.log_density, 0.05) for v ∈ parametervector ], 
    [ quantile(v.log_density, 0.95) for v ∈ parametervector ];
    color=:blue,
)
ax2 = Axis(gb[2, 1]; xticks=( logomegavalues, omegalabels ))
scatter!(
    ax2, logomegavalues, [ quantile(v.β0, 0.5) for v ∈ parametervector ] ./ (48.7 + 0.0087); 
    color=:blue
)
rangebars!(
    ax2, 
    logomegavalues, 
    [ quantile(v.β0, 0.05) for v ∈ parametervector ] ./ (48.7 + 0.0087), 
    [ quantile(v.β0, 0.95) for v ∈ parametervector ] ./ (48.7 + 0.0087);
    color=:blue,
)
ax3 = Axis(gb[3, 1]; xticks=( logomegavalues, omegalabels ))
scatter!(ax3, logomegavalues, [ quantile(v.β1, 0.5) for v ∈ parametervector ]; color=:blue)
rangebars!(
    ax3, 
    logomegavalues, 
    [ quantile(v.β1, 0.05) for v ∈ parametervector ], 
    [ quantile(v.β1, 0.95) for v ∈ parametervector ];
    color=:blue,
)
ax4 = Axis(gb[4, 1]; xticks=( logomegavalues, omegalabels ))
scatter!(ax4, logomegavalues, [ quantile(v.ψ, 0.5) for v ∈ parametervector ]; color=:blue)
rangebars!(
    ax4, 
    logomegavalues, 
    [ quantile(v.ψ, 0.05) for v ∈ parametervector ], 
    [ quantile(v.ψ, 0.95) for v ∈ parametervector ];
    color=:blue,
)
ax5 = Axis(gb[5, 1]; xticks=( logomegavalues, omegalabels ))
scatter!(ax5, logomegavalues, [ quantile(v.detection, 0.5) for v ∈ parametervector ]; color=:blue)
rangebars!(
    ax5, 
    logomegavalues, 
    [ quantile(v.detection, 0.05) for v ∈ parametervector ], 
    [ quantile(v.detection, 0.95) for v ∈ parametervector ];
    color=:blue,
)
linkaxes!(axs...)
for i ∈ 1:7 
    formataxis!(axs[i], hidex=(i != 7), hidexticks=(i != 7))
    if i != 7 hidespines!(axs[i], :b) end
end
Label(ga[1:7, 0], "Weekly incidence"; fontsize = 11.84, rotation = π/2, tellheight = false)
Label(ga[8, 1], "Year"; fontsize = 11.84, tellwidth = false)
for (i, ax) ∈ enumerate([ ax1, ax2, ax3, ax4, ax5 ])
    _sety = i == 1 ? -6000 : 0
    formataxis!(ax, hidex=(i != 5), hidexticks=(i != 5))
    if i != 5 hidespines!(ax, :b) end
    setvalue!(ax, 1, _sety)
end
Label(gb[1, 0], "Log probability"; fontsize = 11.84, rotation = π/2, tellheight = false)
Label(gb[2, 0], "R₀"; fontsize = 11.84, rotation = π/2, tellheight = false)
Label(gb[3, 0], "Magnitude of seasonal forcing"; fontsize = 11.84, rotation = π/2, tellheight = false)
Label(gb[4, 0], "Magnitude of natural immune boosting"; fontsize = 11.84, rotation = π/2, tellheight = false)
Label(gb[5, 0], "Proportion detected"; fontsize = 11.84, rotation = π/2, tellheight = false)
Label(gb[6, 1], "Waning rate, ω"; fontsize = 11.84, tellwidth = false)


fittedparametersfig

