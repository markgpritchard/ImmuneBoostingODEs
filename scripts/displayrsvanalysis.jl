
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
R0_01 = let 
    beta0s = quantile(rsvparameters01.β0, [0.05, 0.5, 0.95])
    beta0s ./ (48.7 + 0.0087)
end
#3-element Vector{Float64}:
# 24.094402995167428
# 27.09900529252659
# 32.83900102448048

rsvparameters02 = loadrsvdata(0.2)
plotchains(rsvparameters02)
plotvals02 = fittedsimulationquantiles(rsvparameters02, 0.2, saveat, cbs)

rsvparameters04 = loadrsvdata(0.4)
plotchains(rsvparameters04)
filter!(:chain => x -> x ∈ [ 1, 2, 5 ], rsvparameters04)
plotchains(rsvparameters04)
plotvals04 = fittedsimulationquantiles(rsvparameters04, 0.4, saveat, cbs)
βreduction1_04 = 1 .- quantile(rsvparameters04.βreduction1, [ 0.95, 0.5, 0.05 ])
#3-element Vector{Float64}:
# 0.42364255789478045      
# 0.4288752645900158       
# 0.44024896448040374  

rsvparameters1 = loadrsvdata(1.0)
plotchains(rsvparameters1)
filter!(:chain => x -> x != 1, rsvparameters1)
plotchains(rsvparameters1)
plotvals1 = fittedsimulationquantiles(rsvparameters1, 1.0, saveat, cbs)
beta1_1 = quantile(rsvparameters1.β1, [ 0.05, 0.5, 0.95 ])
#3-element Vector{Float64}:
# 0.17539991622191103
# 0.1778485156330903
# 0.18209665569439593
psi_1 = quantile(rsvparameters1.ψ, [ 0.05, 0.5, 0.95 ])
#3-element Vector{Float64}:
# 0.00010836786619787666
# 0.0013407470879928376
# 0.005730379235476766

rsvparameters2 = loadrsvdata(2.0)
plotchains(rsvparameters2)
filter!(:chain => x -> x != 4, rsvparameters2)
plotchains(rsvparameters2)
plotvals2 = fittedsimulationquantiles(rsvparameters2, 2.0, saveat, cbs)
βreduction1_2 = 1 .- quantile(rsvparameters2.βreduction1, [ 0.95, 0.5, 0.05 ])
#3-element Vector{Float64}:
# 0.25571515469505746
# 0.26337343069493213
# 0.28768665501096025

rsvparameters4 = loadrsvdata(4.0)
plotchains(rsvparameters4)
filter!(:chain => x -> x ∈ [ 2, 3, 5 ], rsvparameters4)
plotchains(rsvparameters4)
plotvals4 = fittedsimulationquantiles(rsvparameters4, 4.0, saveat, cbs)

rsvparameters10 = loadrsvdata(10.0)
plotchains(rsvparameters10)
filter!(:chain => x -> x ∈ [3, 4 ], rsvparameters10)
plotchains(rsvparameters10)
plotvals10 = fittedsimulationquantiles(rsvparameters10, 10.0, saveat, cbs)
R0_10 = let 
    beta0s = quantile(rsvparameters10.β0, [0.05, 0.5, 0.95])
    beta0s ./ (48.7 + 0.0087)
end
#3-element Vector{Float64}:
# 1.5601342033520065
# 1.570051382100839
# 1.5796336148006065
beta1_10 = quantile(rsvparameters10.β1, [ 0.05, 0.5, 0.95 ])
#3-element Vector{Float64}:
# 0.07060322642648463
# 0.07321548687540379
# 0.07582937390755169
psi_10 = quantile(rsvparameters10.ψ, [0.05, 0.5, 0.95])
#3-element Vector{Float64}:
# 792.804216547885
# 820.3010861906888
# 848.5342792386459
