
using DrWatson
@quickactivate :ImmuneBoostingODEs
using DataFrames, DifferentialEquations, Pigeons, Random, Turing

#using CairoMakie

testrun = true 

if length(ARGS) == 3 
    omega = parse(Float64, ARGS[1])
    id = parse(Int, ARGS[2])
    n_rounds = min(12, parse(Int, ARGS[3]))
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

# RSV data from Scotland
data = processrsvdata("respiratory_scot.csv", "rsv.csv")

# Age-specific data
agedata = processagedata("respiratory_age.csv", "rsv_age.csv")

# Data from Oxford Covid-19 Government Response Tracker
crgtdata = processcrgtvdata("OxCGRT_compact_subnational_v1.csv", "crgt.csv")

# To avoid splitting outbreaks, count cases from April each year 
let 
    april1value = MONTHDAYS[4] / 365
    offsetdate = data.Date .- april1value
    aprilyear = round.(Int, offsetdate, RoundDown)
    aprilfractiondate = offsetdate - aprilyear
    insertcols!(data, :AprilYear => aprilyear)
    insertcols!(data, :AprilFractionDate => aprilfractiondate)
    # insert cumulative cases since last April 
    cumulativecases = Vector{Float64}(undef, size(data, 1))
    cumulativecases[1] = data.Cases[1]
    for i ∈ axes(data, 1)
        i == 1 && continue
        if data.AprilYear[i] == data.AprilYear[i-1]
            cumulativecases[i] = data.Cases[i] + cumulativecases[i-1]
        else 
            cumulativecases[i] = data.Cases[i]
        end 
    end 
    insertcols!(data, :AprilCumulativeCases => cumulativecases)
end 
println("In the 12 months from 1 April each year")
for y ∈ 2016:2022 
    inds = findall(x -> y <= x < y + 1, data.AprilYear)
    println("    $(sum(data.Cases[inds])) cases in $y")
end 


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Fitting parameters 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

saveat = let 
    savefirst = data.Date[1] - 7 / 365  # to allow calculation of new cases in first week
    [ [ savefirst ]; data.Date ]
end

## Callbacks  
cbs = let 
    # When is the infection parameter expected to change?
    # Find dates (as fractions of year) when Strigency Index goes above then below 50
    inds = findall(x -> x >= 50, crgtdata.StringencyIndex_Average)
    reduceday = crgtdata.Date[inds[1]]
    increaseday = crgtdata.Date[last(inds)]
    save_positions = ( false, false )
    resetcb = PresetTimeCallback(reduceday + 1e-9, reducetransmission!; save_positions)
    rescb = PresetTimeCallback(increaseday, restoretransmission!; save_positions)
    CallbackSet(resetcb, rescb)
end

prob = let 
    p = SirnsParameters(50.0, 0.1, 0.0, 48.7, 0.0087, 0.0, 2.0, 50.0, 50.0, 50.0) 
    u0 = sirns_u0(0.01, 2e-5; p, equalrs=true, t0=1996.737)
    tspan = ( 1996.737, last(saveat) )
    ODEProblem(sirns!, u0, tspan, p)
end

@model function fitmodel(
    incidence, prob, cbs, saveat;
    betazeroprior=truncated(Exponential(100), 0, 4870),  # truncated at R0 = 100
    betaoneprior=Uniform(0, 0.9),
    phiprior=Uniform(0, 2π),
    psiprior=truncated(Exponential(1), 0, 1000),
    betareduction1prior=Beta(4, 1),
    betareduction2prior=Beta(9, 1),
    omega=2.0,
    detectionprior=Beta(1, 98)
)
    β0 ~ betazeroprior
    β1 ~ betaoneprior
    ϕ ~ phiprior
    γ = 48.7  # generation time 7.5 days
    μ = 0.0087  # Scotland's birth rate = 48000 / 5.5e6
    ψ ~ psiprior
    ω = omega
    βreduction1 ~ betareduction1prior
    βreduction2 ~ betareduction2prior
    detection ~ detectionprior

    p = SirnsParameters(β0, β1, ϕ, γ, μ, ψ, ω, β0, βreduction1 * β0, βreduction2 * β0)
    u0 = sirns_u0(0.01, 2e-5; p, equalrs=true, t0=1996.737)  # 10 years before data collection

    sol = memosolver(
        prob, Vern9(; lazy=false); 
        p, u0, callback=cbs, saveat, abstol=1e-15, maxiters=1e7,
    )
    if sol.retcode != :Success
        @info "Adding logprob -Inf when p=$p, detection=$detection"
        Turing.@addlogprob! -Inf
        return nothing
    end

    cumulativecases = modelcompartments(sol, :cc)
    incidentcases = casespertimeblock(cumulativecases) .* 5_450_000 .* detection

    for i ∈ eachindex(incidentcases)
        incidence[i] ~ Poisson(incidentcases[i] + 1e-10)
    end
end

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

fitted_pt = pigeons( ;
    target=fitmodel_target(; omega), 
    n_rounds=0,
    n_chains=6,
    multithreaded=true,
    record=[ traces; record_default() ],
    seed=id,
    variational=GaussianReference(),
)

new_pt = fitted_pt

for i ∈ 1:n_rounds
    filename = "rsvparameters_omega_$(omega)_id_$(id)_nrounds_$(i).jld2"
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
