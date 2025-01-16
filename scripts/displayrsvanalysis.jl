
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
filter!(:chain => x -> x ∈ [ 2, 4, 5 ], rsvparameters01)
plotchains(rsvparameters01)
plotvals01 = fittedsimulationquantiles(rsvparameters01, 0.1, saveat, cbs)
R0_01 = let 
    beta0s = quantile(rsvparameters01.β0, [ 0.05, 0.5, 0.95 ])
    beta0s ./ (48.7 + 0.0087)
end
#3-element Vector{Float64}:
# 27.9398281195748
# 29.58401062939877
# 30.730225320740598

rsvparameters02 = loadrsvdata(0.2)
plotchains(rsvparameters02)
filter!(:chain => x -> x == 1, rsvparameters02)
plotchains(rsvparameters02)
plotvals02 = fittedsimulationquantiles(rsvparameters02, 0.2, saveat, cbs)

rsvparameters04 = loadrsvdata(0.4)
plotchains(rsvparameters04)
filter!(:chain => x -> x == 1, rsvparameters04)
plotchains(rsvparameters04)
plotvals04 = fittedsimulationquantiles(rsvparameters04, 0.4, saveat, cbs)
βreduction1_04 = 1 .- quantile(rsvparameters04.βreduction1, [ 0.95, 0.5, 0.05 ])
#3-element Vector{Float64}:
# 0.42537389746949894
# 0.4274937278262533
# 0.4304722062246449

rsvparameters1 = loadrsvdata(1.0)
plotchains(rsvparameters1)
filter!(:chain => x -> x != 1, rsvparameters1)
plotchains(rsvparameters1)
plotvals1 = fittedsimulationquantiles(rsvparameters1, 1.0, saveat, cbs)
beta1_1 = quantile(rsvparameters1.β1, [ 0.05, 0.5, 0.95 ])
#3-element Vector{Float64}:
# 0.1754460107922378
# 0.17788693495654073
# 0.1818164970772631
psi_1 = quantile(rsvparameters1.ψ, [ 0.05, 0.5, 0.95 ])
#3-element Vector{Float64}:
# 0.00010909782961219074
# 0.0013969680682125209
# 0.005838432275809601

rsvparameters2 = loadrsvdata(2.0)
plotchains(rsvparameters2)
filter!(:chain => x -> x ∈ [ 2, 3, 5 ], rsvparameters2)
plotchains(rsvparameters2)
plotvals2 = fittedsimulationquantiles(rsvparameters2, 2.0, saveat, cbs)
βreduction1_2 = 1 .- quantile(rsvparameters2.βreduction1, [ 0.95, 0.5, 0.05 ])
#3-element Vector{Float64}:
# 0.27889082920822683
# 0.2816906091386049
# 0.28507902779791494

rsvparameters4 = loadrsvdata(4.0)
plotchains(rsvparameters4)
filter!(:chain => x -> x ∈ [ 2, 3, 5 ], rsvparameters4)
plotchains(rsvparameters4)
plotvals4 = fittedsimulationquantiles(rsvparameters4, 4.0, saveat, cbs)

rsvparameters6 = loadrsvdata(6.0)
plotchains(rsvparameters6)
filter!(:chain => x -> x != 1, rsvparameters6)
plotchains(rsvparameters6)
plotvals6 = fittedsimulationquantiles(rsvparameters6, 6.0, saveat, cbs)
R0_6 = let 
    beta0s = quantile(rsvparameters6.β0, [ 0.05, 0.5, 0.95 ])
    beta0s ./ (48.7 + 0.0087)
end
#3-element Vector{Float64}:
# 1.5808061809771425
# 1.592694730925088
# 1.6042283249278522
beta1_6 = quantile(rsvparameters6.β1, [ 0.05, 0.5, 0.95 ])
#3-element Vector{Float64}:
# 0.06894139362433213
# 0.07223981399834295
# 0.07529313242053985
psi_6 = quantile(rsvparameters6.ψ, [ 0.05, 0.5, 0.95 ])
#3-element Vector{Float64}:
# 288.26257743152075
# 301.2348000922061
# 314.895983917204
