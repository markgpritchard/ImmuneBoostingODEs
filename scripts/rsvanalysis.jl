
using DrWatson

@quickactivate "ImmuneBoostingODEs"
using Pkg 
Pkg.instantiate()

@quickactivate :ImmuneBoostingODEs
using DataFrames, DifferentialEquations, Pigeons, Random, Turing

testrun = true 

if length(ARGS) == 3 
    omega = parse(Float64, ARGS[1])
    id = parse(Int, ARGS[2])
    n_rounds = min(10, parse(Int, ARGS[3]))
else
    omega = 2.0
    id = 1 
    if testrun 
        n_rounds = 4 
    else
        n_rounds = 10
    end
end

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load the data 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

include("rsvsetup.jl")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Fitting parameters 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

prob = fittedsimulationsetup(saveat)

include("rsvfitmodel.jl")

function fitmodel_target(incidence=data.Cases, prob=prob, cbs=cbs, saveat=saveat; kwargs...)
    return Pigeons.TuringLogPotential(fitmodel(incidence, prob, cbs, saveat; kwargs...))
end

const FitmodelType = typeof(fitmodel_target())

function Pigeons.initialization(target::FitmodelType, rng::AbstractRNG, ::Int64)
    result = DynamicPPL.VarInfo(
        rng, target.model, DynamicPPL.SampleFromPrior(), DynamicPPL.PriorContext()
    )
    DynamicPPL.link!!(result, DynamicPPL.SampleFromPrior(), target.model)

    Pigeons.update_state!(result, :β0, 1, 2.0)
    Pigeons.update_state!(result, :β1, 1, 0.1)
    Pigeons.update_state!(result, :ϕ, 1, 0.0)
    Pigeons.update_state!(result, :ψ, 1, 0.0)
    Pigeons.update_state!(result, :βreduction1, 1, 0.9)
    Pigeons.update_state!(result, :βreduction1, 1, 1.0)
    Pigeons.update_state!(result, :detection, 1, 0.02)

    return result
end

const seed = (round(Int, omega * 100) + id)

fitted_pt = pigeons( ;
    target=fitmodel_target(; omega), 
    n_rounds=0,
    n_chains=10,
    multithreaded=true,
    record=[ traces; record_default() ],
    seed,
    variational=GaussianReference(),
)

new_pt = fitted_pt

for i ∈ 1:n_rounds
    filename = "rsvparameters_omega_$(omega)_seed_$(seed)_id_$(id)_nrounds_$(i).jld2"
    nextfilename = "rsvparameters_omega_$(omega)_seed_$(seed)_id_$(id)_nrounds_$(i + 1).jld2"
    isfile(datadir("sims", nextfilename)) && continue
    if isfile(datadir("sims", filename))
        global new_pt = load(datadir("sims", filename))["pt"]
    else
        pt = increment_n_rounds!(new_pt, 1)
        global new_pt = pigeons(pt)
        new_chains = Chains(new_pt)
        resultdict = Dict(
            "chain" => new_chains, 
            "pt" => new_pt, 
            "n_rounds" => i, 
            "n_chains" => 6,
        )
        safesave(datadir("sims", filename), resultdict)
    end
end
