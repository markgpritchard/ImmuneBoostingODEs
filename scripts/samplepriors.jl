
using DataFrames, DifferentialEquations, Random, Turing

include("rsvsetup.jl")
include("rsvfitmodel.jl")

Random.seed!(1729)

if isfile(datadir("sims", "priorsdict.jld2"))
    @info "Loading prior values"
    priorsdict = load(datadir("sims", "priorsdict.jld2"))
else
    priorsdict = let 
        pd = Dict{String, Chains}()
        names = [ "0.1", "0.2", "0.4", "1.0", "2.0", "4.0", "6.0" ]
        omegas = [ 0.1, 0.2, 0.4, 1.0, 2.0, 4.0, 6.0 ]
        for i ∈ 1:7 
            @info "Sampling for ω=$(omegas[i])"
            prob = fittedsimulationsetup(saveat)
            m = fitmodel(data.Cases, prob, cbs, saveat; omega=omegas[i])
            p = sample(m, Prior(), MCMCThreads(), 250, 4)
            push!(pd, names[i] => p)
        end
        pd
    end
    safesave(datadir("sims", "priorsdict.jld2"), priorsdict)
end
