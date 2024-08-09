
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

rsvparameters025 = loadrsvdata(0.25)
plotchains(rsvparameters025)
plotvals025 = fittedsimulationquantiles(rsvparameters025, 0.25, saveat, cbs)

rsvparameters1 = loadrsvdata(1.0)
plotchains(rsvparameters1)
plotvals1 = fittedsimulationquantiles(rsvparameters1, 1.0, saveat, cbs)

rsvparameters2 = loadrsvdata(2.0)
plotchains(rsvparameters2)
plotvals2 = fittedsimulationquantiles(rsvparameters2, 2.0, saveat, cbs)

rsvparameters4 = loadrsvdata(4.0)
plotchains(rsvparameters4)
plotvals4 = fittedsimulationquantiles(rsvparameters4, 4.0, saveat, cbs)

rsvparameters26 = loadrsvdata(26.0)
plotchains(rsvparameters26)
plotvals26 = fittedsimulationquantiles(rsvparameters26, 26.0, saveat, cbs)

plotvvector = [ plotvals01, plotvals025, plotvals1, plotvals2, plotvals4, plotvals26 ]
parametervector = [ 
    rsvparameters01, rsvparameters025, rsvparameters1, 
    rsvparameters2, rsvparameters4, rsvparameters26 
]
omegalabels = [ "0.1", "0.25", "1.0", "2.0", "4.0", "26.0" ]

fittedparametersfig = Figure(; size=( 800, 1200 ))
ga = GridLayout(fittedparametersfig[1, 1])
axs = [ Axis(ga[i, 1]) for i ∈ 1:6]
for (i, v) ∈ enumerate(plotvvector)
    plotfittedsimulationquantiles!(axs[i], data, v, saveat)
end
gb = GridLayout(fittedparametersfig[1, 2])
ax1 = Axis(gb[1, 1]; xticks=( 1:6, omegalabels ))
scatter!(ax1, 1:6, [ quantile(v.log_density, 0.5) for v ∈ parametervector ]; color=:blue)
rangebars!(
    ax1, 
    1:6, 
    [ quantile(v.log_density, 0.05) for v ∈ parametervector ], 
    [ quantile(v.log_density, 0.95) for v ∈ parametervector ];
    color=:blue,
)
ax2 = Axis(gb[2, 1]; xticks=( 1:6, omegalabels ))
scatter!(
    ax2, 1:6, [ quantile(v.β0, 0.5) for v ∈ parametervector ] ./ (48.7 + 0.0087); 
    color=:blue
)
rangebars!(
    ax2, 
    1:6, 
    [ quantile(v.β0, 0.05) for v ∈ parametervector ] ./ (48.7 + 0.0087), 
    [ quantile(v.β0, 0.95) for v ∈ parametervector ] ./ (48.7 + 0.0087);
    color=:blue,
)
ax3 = Axis(gb[3, 1]; xticks=( 1:6, omegalabels ))
scatter!(ax3, 1:6, [ quantile(v.β1, 0.5) for v ∈ parametervector ]; color=:blue)
rangebars!(
    ax3, 
    1:6, 
    [ quantile(v.β1, 0.05) for v ∈ parametervector ], 
    [ quantile(v.β1, 0.95) for v ∈ parametervector ];
    color=:blue,
)
ax4 = Axis(gb[4, 1]; xticks=( 1:6, omegalabels ))
scatter!(ax4, 1:6, [ quantile(v.ψ, 0.5) for v ∈ parametervector ]; color=:blue)
rangebars!(
    ax4, 
    1:6, 
    [ quantile(v.ψ, 0.05) for v ∈ parametervector ], 
    [ quantile(v.ψ, 0.95) for v ∈ parametervector ];
    color=:blue,
)
ax5 = Axis(gb[5, 1]; xticks=( 1:6, omegalabels ))
scatter!(ax5, 1:6, [ quantile(v.detection, 0.5) for v ∈ parametervector ]; color=:blue)
rangebars!(
    ax5, 
    1:6, 
    [ quantile(v.detection, 0.05) for v ∈ parametervector ], 
    [ quantile(v.detection, 0.95) for v ∈ parametervector ];
    color=:blue,
)
linkaxes!(axs...)
for i ∈ 1:6 
    formataxis!(axs[i], hidex=(i != 6), hidexticks=(i != 6))
    if i != 6 hidespines!(axs[i], :b) end
end
Label(ga[1:6, 0], "Weekly incidence"; fontsize = 11.84, rotation = π/2, tellheight = false)
Label(ga[7, 1], "Year"; fontsize = 11.84, tellwidth = false)


fittedparametersfig

