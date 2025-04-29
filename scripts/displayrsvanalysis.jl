
using DrWatson
@quickactivate :ImmuneBoostingODEs

using CairoMakie, DataFrames, DifferentialEquations, Random, Turing

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

chain_omega_01 = load(datadir("sims", "chain_omega_0.1_nrounds_25.jld2"))
chain_omega_01df = DataFrame(chain_omega_01["chain"])
plotchains(chain_omega_01df)
plotvals01 = fittedsimulationquantiles(
    chain_omega_01df, 0.1, saveat, cbs, [ 0.025, 0.5, 0.975 ]
)

chain_omega_02 = load(datadir("sims", "chain_omega_0.2_nrounds_25.jld2"))
chain_omega_02df = DataFrame(chain_omega_02["chain"])
plotchains(chain_omega_02df)
plotvals02 = fittedsimulationquantiles(
    chain_omega_02df, 0.2, saveat, cbs, [ 0.025, 0.5, 0.975 ]
)

chain_omega_04 = load(datadir("sims", "chain_omega_0.4_nrounds_25.jld2"))
chain_omega_04df = DataFrame(chain_omega_04["chain"])
plotchains(chain_omega_04df)
plotvals04 = fittedsimulationquantiles(
    chain_omega_04df, 0.4, saveat, cbs, [ 0.025, 0.5, 0.975 ]
)

chain_omega_1 = load(datadir("sims", "chain_omega_1.0_nrounds_25.jld2"))
chain_omega_1df = DataFrame(chain_omega_1["chain"])
plotchains(chain_omega_1df)
plotvals1 = fittedsimulationquantiles(
    chain_omega_1df, 1.0, saveat, cbs, [ 0.025, 0.5, 0.975 ]
)

chain_omega_2 = load(datadir("sims", "chain_omega_2.0_nrounds_25.jld2"))
chain_omega_2df = DataFrame(chain_omega_2["chain"])
plotchains(chain_omega_2df)
plotvals2 = fittedsimulationquantiles(
    chain_omegadf, 2.0, saveat, cbs, [ 0.025, 0.5, 0.975 ]
)

chain_omega_4 = load(datadir("sims", "chain_omega_4.0_nrounds_25.jld2"))
chain_omega_4df = DataFrame(chain_omega_4["chain"])
plotchains(chain_omega_4df)
plotvals4 = fittedsimulationquantiles(
    chain_omega_4df, 4.0, saveat, cbs, [ 0.025, 0.5, 0.975 ]
)

chain_omega_6 = load(datadir("sims", "chain_omega_6.0_nrounds_25.jld2"))
chain_omega_6df = DataFrame(chain_omega_6["chain"])
plotchains(chain_omega_6df)
plotvals6 = fittedsimulationquantiles(
    chain_omega_6df, 6.0, saveat, cbs, [ 0.025, 0.5, 0.975 ]
)
