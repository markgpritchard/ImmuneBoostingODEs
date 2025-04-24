
using DrWatson

@quickactivate :ImmuneBoostingODEs
using DataFrames, DifferentialEquations, DynamicPPL, Optim, Pigeons, Random, Turing
import Pigeons: initialization

testrun = true 

if length(ARGS) == 3 
    const omega = parse(Float64, ARGS[1])
    const id = parse(Int, ARGS[2])
    n_rounds = parse(Int, ARGS[3])
    @assert id >= 1
else
    const omega = 0.1#2.0
    const id = 1 
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

const constcases = data.Cases

function optimmodel(paramvector)
    # take any real values in paramvector and transform into "acceptable" values 
    β0, ψ, _detection = paramvector
    #β0 = exp(_β0)
    #ψ = 10_000 * exp(_ψ) / (1 + exp(_ψ))
    detection = 0.001 + 0.1 * exp(_detection) / (1 + exp(_detection))
    if β0 < 48.7 || ψ < 0 || detection < 0.001
        println("fail′ with x=$paramvector, β0=$β0, ψ=$ψ, detection=$detection")
        return Inf
    end
    #if minimum(paramvector) < 0 
    #    println("fail′ with x=$paramvector")
    #    return 2.0^63 
    #end
    p = SirnsParameters(
        β0, 
        0.1,  # β1 
        0.0,  # ϕ 
        48.7,  # γ 
        0.0087,  # μ 
        ψ, 
        omega,  # declared globally above
        β0, 
        0.8 * β0, 
        0.9 * β0
    )
    u0 = sirns_u0(0.01, 2e-5; p, equalrs=true, t0=1996.737)  # 10 years before data collection
    sol = memosolver(
        prob, Vern9(; lazy=false); 
        p, u0, callback=cbs, saveat, save_idxs=[ 8 ], 
        abstol=1e-15, maxiters=1e8, verbose=false,
    )
    if sol.retcode != :Success
        println("fail with x=$paramvector, β0=$β0, ψ=$ψ, detection=$detection")
        return 2.0^63 + abs2(β0 + ψ)
    end
    cumulativecases = modelcompartments(sol, 1)
    incidentcases = casespertimeblock(cumulativecases) .* 5_450_000 .* detection
    return sum(
        [ 
            constcases[i] < 1 ? 
                abs2(constcases[i] - incidentcases[i]) :
                abs2(constcases[i] - incidentcases[i]) / constcases[i]
            for i ∈ eachindex(constcases) 
        ]
    )
end

result = optimize(optimmodel, [ 974.0, 0.1, -2.4941 ], NelderMead())

const initvalues = Optim.minimizer(result)

# Success with x=[326082.84130647726, 26841.64263241941, -3.933717402561181e7], β0=326082.84130647726, ψ=26841.64263241941, detection=0.001
# Result = 25893.340013631245

function fitmodel_target(incidence=data.Cases, prob=prob, cbs=cbs, saveat=saveat; kwargs...)
    return Pigeons.TuringLogPotential(fitmodel(incidence, prob, cbs, saveat; kwargs...))
end

const FitmodelType = typeof(fitmodel_target())

function Pigeons.initialization(target::FitmodelType, rng::AbstractRNG, ::Int64)
    result = DynamicPPL.VarInfo(
        rng, target.model, DynamicPPL.SampleFromPrior(), DynamicPPL.PriorContext()
    )
    β0, ψ, _detection = paramvector
    detection = 0.001 + 0.1 * exp(_detection) / (1 + exp(_detection))

    DynamicPPL.link!!(result, DynamicPPL.SampleFromPrior(), target.model)

    Pigeons.update_state!(result, :β0, 1, β0)
    Pigeons.update_state!(result, :β1, 1, 0.1)
    Pigeons.update_state!(result, :ϕ, 1, 0.0)
    Pigeons.update_state!(result, :ψ, 1, ψ)
    Pigeons.update_state!(result, :βreduction1, 1, 0.8 * β0)
    Pigeons.update_state!(result, :βreduction1, 1, 0.9 * β0)
    Pigeons.update_state!(result, :detection, 1, detection)

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
            "n_chains" => 10,
        )
        safesave(datadir("sims", filename), resultdict)
    end
end
