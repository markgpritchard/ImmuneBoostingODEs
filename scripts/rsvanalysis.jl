
using DrWatson
@quickactivate :ImmuneBoostingODEs
using DataFrames, DifferentialEquations, Distributions, Pigeons, Random, Turing

testrun = true 

if length(ARGS) == 2 
    id = parse(Int, ARGS[1])
    n_rounds = min(10, parse(Int, ARGS[2]))
else
    id = 1 
    if testrun 
        n_rounds = 4 
    else
        n_rounds = 10
    end
end

include("rsvanalysisconsts.jl")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Model set-up 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

p = SirrrsParameters(0.2, contacts_alllocations, 48.7, births[1], mortality, 1.0, 0.913)

u0 = let 
    props = [ [ 0.3 ]; ones(20) .* 0.01; ones(4) .* 0.005; ones(2) .* 0.007; [ 0.0125 ] ]
    u = zeros(6, 28)
    for a ∈ 1:28
        u[1, a] = props[a] * pops[a] 
        u[2, a] = 0.01 * pops[a] 
        for i ∈ 3:5 
            u[i, a] = (pops[a] - props[a] * pops[a] - 0.01 * pops[a]) / 3
        end
    end
    u
end

tspan = ( 2006.581, 2023.5233 )

prob = ODEProblem(sirrrs!, u0, tspan, p)
sol = solve(
    prob, Vern9(lazy=false); 
    p, u0, callback=cbs, saveat, abstol=1e-15, maxiters=5e4, 
)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Parameter fitting model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

@model function fitmodel(data, prob, u0, cbs, mortality)
    τ ~ truncated(Exponential(0.1), 0.0, 10.0)
    ψ ~ truncated(Exponential(1), 0.0, 100.0) 
    ϕ ~ Beta(1, 1)

    p = SirrrsParameters(τ, contacts_alllocations, 48.7, births[1], mortality, ψ, 0.913)
    sol = memosolver(
        prob, Vern9(lazy=false); 
        p, u0, callback=cbs, saveat, abstol=1e-15, maxiters=5e4,
    )
    if sol.retcode != :Success
        @info "Adding logprob -Inf when τ=$τ, ψ=$ψ, ϕ=$ϕ, sol.retcode=$(sol.retcode)"
        Turing.@addlogprob! -Inf
        return nothing
    end
    
    incidentcases = makeoutputmatrix_incidentcases(sol)

    for i ∈ axes(incidentcases, 1), j ∈ axes(incidentcases, 2)
        i == 1 && continue  # first row of incident cases always 0
        data[i, j] ~ Poisson(incidentcases[i, j] * ϕ + 1e-10)
    end
end

function fitmodel_target(
    data=intagedatamatrix, prob=prob, u0=u0, cbs=cbs, mortality=mortality
)
    return Pigeons.TuringLogPotential(fitmodel(data, prob, u0, cbs, mortality))
end

const FitmodelType = typeof(fitmodel_target())

function Pigeons.initialization(target::FitmodelType, rng::AbstractRNG, ::Int64)
    result = DynamicPPL.VarInfo(
        rng, target.model, DynamicPPL.SampleFromPrior(), DynamicPPL.PriorContext()
    )
    DynamicPPL.link!!(result, DynamicPPL.SampleFromPrior(), target.model)

    Pigeons.update_state!(result, :τ, 1, 0.2)
    Pigeons.update_state!(result, :ψ, 1, 1.0)
    Pigeons.update_state!(result, :ϕ, 1, 0.005)

    return result
end


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Fitting
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

fitted_pt = pigeons( ;
    target=fitmodel_target(), 
    n_rounds=0,
    n_chains=16,
    multithreaded=true,
    record=[ traces; record_default() ],
    seed=id,
    variational=GaussianReference(),
)

@assert Pigeons.variable(fitted_pt.replicas[1].state, :τ) == [ 0.2 ]
@assert Pigeons.variable(fitted_pt.replicas[1].state, :ψ) == [ 1 ]

new_pt = fitted_pt

for i ∈ 1:n_rounds
    filename = "rsvparameters_id_$(id)_nrounds_$(i).jld2"
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
            "n_chains" => 16,
        )
        safesave(datadir("sims", filename), resultdict)
    end
end
