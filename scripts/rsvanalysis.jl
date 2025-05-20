
using DrWatson

@quickactivate :ImmuneBoostingODEs

using DataFrames
using DifferentialEquations
using Memoization
using Optimization
#using OptimizationBBO
using OptimizationOptimJL
#using OptimizationPolyalgorithms
using OrdinaryDiffEq
using Random
#using ReverseDiff
#using SciMLSensitivity
using Turing
using Zygote

#include("samplepriors.jl")

testrun = true 

if length(ARGS) == 2 
    n_rounds = parse(Int, ARGS[1])
    optimizationsolvermaxiters = parse(Int, ARGS[2])
else
    n_rounds = testrun ? 25 : 10_000
    optimizationsolvermaxiters = testrun ? 25_000 : 1e6
end

println("Starting with n_rounds=$n_rounds, optimizationsolvermaxiters=$optimizationsolvermaxiters")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load the data 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

include("rsvsetup.jl")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Fitting parameters 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function loss(
    parms;  # a vector containing the following in order:
        # logr0, 
        # logitβ1, 
        # ϕ, 
        # logψ, 
        # logω, 
        # logitbetaprimemultiplier 
        # logitfinalbetaprime, 
        # logitproportiondetected, 
        # <any non-negative number>
        # logitS0max
        # logitI0, 
    data, 
    prob, 
    callback, 
    saveat,
    equalrs=true,
    t0=1996.737,
    alg=Vern9(; lazy=false),
    abstol=1e-15,
    maxiters=1e8,
    gamma=48.7,
    mu=0.0087,
) 
    if parms[3] < -π || parms[3] > π || parms[9] < 0
        return Inf
    end

    I0 = ImmuneBoostingODEs._logistic(parms[11])
    S0 = min(ImmuneBoostingODEs._logistic(parms[10]), 1 - I0)
    u0 = sirns_u0_transformedp(S0, I0; p=parms, equalrs, t0)
    sol = solve(prob, alg; p=parms, u0, callback, saveat, abstol, maxiters,)

    if !SciMLBase.successful_retcode(sol)
        @warn "$(sol.retcode) with parms=$parms, u0=$u0"
        return Inf
    end

    cumulativecases = modelcompartments(sol, 8)
    incidentcases = casespertimeblock(cumulativecases) .* 5_450_000
    loss = sum(
        [
            -log(
                pdf(
                    NegativeBinomial(
                        0.25,
                        0.25 / (0.25 + incidentcases[t] * ImmuneBoostingODEs._logistic(parms[8]))
                    ),
                    data[t]
                )
            )
            for t ∈ eachindex(data)
        ]
    ) - abs2(ImmuneBoostingODEs._logistic(parms[7]) - parms[9])  # so that parms[9] has some influence on `loss` 
    return loss
end

function optimizesirns(
    data,
    parms;
    callback,
    saveat,
    equalrs=true,
    t0=1996.737,
    tspan=( 1996.737, last(saveat) ),
    gamma=48.7,
    mu=0.0087,
    alg=Vern9(; lazy=false),
    odesolverabstol=1e-15,
    odesolvermaxiters=5e7,
    optimizationsolvermaxiters=1e5,
    adtype=Optimization.AutoZygote(),
    lb=[ -1, -5.3, -2, -4, -0.7, -2, 0.24, -6.6, 0, -4.3, -10.3 ],
    ub=[ 3, 0, 2, 3, 1.1, 0.2, 4.2, -2.6, 1, 4.3, -1.7 ],
    nt=10,
    rt=0.975,
    r_expand=2.0,
    verbosity=3,
)
    I0 = ImmuneBoostingODEs._logistic(parms[11])
    S0 = min(ImmuneBoostingODEs._logistic(parms[10]), 1 - I0)
    u0 = sirns_u0_transformedp(S0, I0; p=parms, equalrs, t0)
    prob = ODEProblem(transformedsirns!, u0, tspan, parms)
    sol = solve(
        prob, alg; 
        p=parms, u0, callback, saveat, abstol=odesolverabstol, maxiters=odesolvermaxiters,
    )
    
    if !SciMLBase.successful_retcode(sol)
        @error "$(sol.retcode) with parms=$parms, u0=$u0"
        return nothing
    end
    
    optf = Optimization.OptimizationFunction(
        (x, p) -> loss(
            x; 
            callback, 
            data, 
            prob, 
            saveat,  
            equalrs, 
            t0, 
            alg, 
            abstol=odesolverabstol, 
            maxiters=odesolvermaxiters,
        ), 
        adtype
    )
    optprob = Optimization.OptimizationProblem(optf, parms; lb, ub)
    result_ode = Optimization.solve(
        optprob, Optim.SAMIN(; nt, rt, r_expand, verbosity); 
        maxiters=optimizationsolvermaxiters
    )

    return result_ode
end

adjustedparams(p) = [ p[1:8]; 4.0; p[10]; p[11] .+ 6 ]

initial_params = Vector{Vector{Float64}}(undef, 8)

Threads.@threads for i ∈ 1:8 
    Random.seed!(n_rounds + optimizationsolvermaxiters + i)

    lb = [ -1.0, -Inf, -0.8, -4.0, -Inf, -2.0, -Inf, -3.892, 0.0, -4.3, -10.3 ]
    ub = [ 3.0, Inf, 0.8, 3.0, Inf, 0.2, Inf, -3.892, 1.0, 4.3, -1.7 ]
    
    if isodd(i)
        lb[2] = ub[2] = log(0.01 / 0.99)
    else
        lb[2] = ub[2] = log(0.25 / 0.75)
    end

    if i ∈ [ 1, 2, 5, 6 ]
        lb[5] = ub[5] = log(2)
    else
        lb[5] = ub[5] = log(0.5)
    end

    if i <= 4
        lb[7] = ub[7] = log(0.75 / 0.25)
    else
        lb[7] = ub[7] = log(0.95 / 0.05)
    end

    _add = (ub .- lb) .* 0.5
    ip = optimizesirns(
        data.Cases, lb .+ _add;
        callback=optimcbs, saveat, optimizationsolvermaxiters, lb, ub, verbosity=0
    )
    initial_params[i] = adjustedparams(ip.minimizer)
    @info "initial_params[$i] $(ip.retcode)"
end

tspan = ( 1996.737, last(saveat) )
initialp = SirnsParameters(
    2.0 * (48.7 + 0.0087),  # β0::T
    0.1,  # β1::T
    0.0,  # ϕ::T
    48.7,  # γ::Float64
    0.0087,  # μ::Float64 
    1.0,  # ψ::T
    1.0,  # ω::T
    2.0 * (48.7 + 0.0087),  # originalβ0::T
    0.5,  # betaprimemultiplier::T
    0.5,  # finalbetaprime::T
    0.5,  # proportiondetected::T
)
u0 = sirns_u0(0.01, 2e-5; p=initialp, equalrs=true, t0=1996.737)
prob = ODEProblem(sirns!, u0, tspan, initialp)

include("rsvfitmodel.jl")
#=
Random.seed!(1729)
priorschain = sample(fitmodel(data.Cases, prob, cbs, saveat), Prior(), MCMCThreads(), 250, 4)
priorschaindf = DataFrame(priorschain)
plotchains(priorschaindf)

using CairoMakie

priormodeloutputs = Vector{Vector{<:Union{Float64, Missing}}}(undef, size(priorschaindf, 1))

for i ∈ axes(priorschaindf, 1)
    r0 = exp(priorschaindf.logr0[i])
    p = SirnsParameters(
        r0 * (48.7 + 0.0087),  # β0::T
        _logistic(priorschaindf.logitβ1[i]),  # β1::T
        priorschaindf.ϕ[i],  # ϕ::T
        48.7,  # γ::Float64
        0.0087,  # μ::Float64 
        exp(priorschaindf.logψ[i]),  # ψ::T
        exp(priorschaindf.logω[i]),  # ω::T
        r0 * (48.7 + 0.0087),  # originalβ0::T
        _logistic(priorschaindf.logitbetaprimemultiplier[i]),  # betaprimemultiplier::T
        _logistic(priorschaindf.logitfinalbetaprime[i]),  # finalbetaprime::T
        _logistic(priorschaindf.logitproportiondetected[i]),  # proportiondetected::T
    )
    I0 = _logistic(priorschaindf.transformedlogitI0[i] - 6)
    S0 = min(ImmuneBoostingODEs._logistic(priorschaindf.logitS0max[i]), 1 - I0)
    u0 = sirns_u0(S0, I0; p, equalrs=true, t0=1996.737)
    tspan = ( 1996.737, last(saveat) )
    
    prob = ODEProblem(sirns!, u0, tspan, p) 
    sol = solve(
        prob, Vern9(; lazy=false); 
        p, u0, callback=cbs, saveat, abstol=1e-15, maxiters=1e8,
    )

    if SciMLBase.successful_retcode(sol)
        cumulativecases = modelcompartments(sol, 8)
        incidentcases = casespertimeblock(cumulativecases .* 5_450_000 * p.proportiondetected)
        priormodeloutputs[i] = incidentcases
    else
        @warn "$(sol.retcode) with p=$p, u0=$u0"
        priormodeloutputs[i] = missings(353)
    end
end

medianoutput = Vector{Float64}(undef, 353)
lcioutput = Vector{Float64}(undef, 353)
ucioutput = Vector{Float64}(undef, 353)
for t ∈ 1:353 
    #lc, me, uc = quantile(skipmissing([ priormodeloutputs[i][t] for i ∈ axes(priorschaindf, 1) ]), [ 0.025, 0.5, 0.975 ])
    lc, me, uc = quantile(
        skipmissing([ priormodeloutputs[i][t] for i ∈ axes(priorschaindf, 1) ]), 
        [ 0.05, 0.5, 0.95 ]
    )
    medianoutput[t] = me
    lcioutput[t] = lc 
    ucioutput[t] = uc 
end

fig = Figure(; size=( 500, 500 ))
ga = GridLayout(fig[1, 1])
ax = Axis(ga[1, 1]; xticks=2017:2:2023,)# yticks=0:200:600)
band!(ax, data.Date, lcioutput, ucioutput; color=( COLOURVECTOR[1], 0.5 ))
lines!(ax, data.Date, medianoutput; color=COLOURVECTOR[1], linewidth=1,)
scatter!(ax, data.Date, data.Cases; color=:black, markersize=3)


fig

=#

 
chain = sample(
    fitmodel(data.Cases, prob, cbs, saveat),
    Turing.NUTS(0.65),
    MCMCThreads(),
    n_rounds,
    8;
    initial_params
)
#=
chaindf = DataFrame(chain)
plotchains(chaindf)
=#
chaindict = Dict(
    "chain" => chain,
    "initial_params" => initial_params,
    "n_rounds" => n_rounds,
    "optimizationsolvermaxiters" => optimizationsolvermaxiters,
)

safesave(datadir("sims", "chaindict_nrounds_$(n_rounds).jld2"), chaindict)
#=
modeloutputs = Vector{Vector{Float64}}(undef, size(chaindf, 1))

for i ∈ axes(chaindf, 1)
    r0 = exp(chaindf.logr0[i])
    p = SirnsParameters(
        r0 * (48.7 + 0.0087),  # β0::T
        _logistic(chaindf.logitβ1[i]),  # β1::T
        chaindf.ϕ[i],  # ϕ::T
        48.7,  # γ::Float64
        0.0087,  # μ::Float64 
        exp(chaindf.logψ[i]),  # ψ::T
        exp(chaindf.logω[i]),  # ω::T
        r0 * (48.7 + 0.0087),  # originalβ0::T
        _logistic(chaindf.logitbetaprimemultiplier[i]),  # betaprimemultiplier::T
        _logistic(chaindf.logitfinalbetaprime[i]),  # finalbetaprime::T
        _logistic(chaindf.logitproportiondetected[i]),  # proportiondetected::T
    )
    I0 = _logistic(chaindf.transformedlogitI0[i] - 6)
    S0 = min(ImmuneBoostingODEs._logistic(chaindf.logitS0max[i]), 1 - I0)
    u0 = sirns_u0(S0, I0; p, equalrs=true, t0=1996.737)
    tspan = ( 1996.737, last(saveat) )
    
    prob = ODEProblem(sirns!, u0, tspan, p) 
    sol = solve(
        prob, Vern9(; lazy=false); 
        p, u0, callback=cbs, saveat, abstol=1e-15, maxiters=1e8,
    )

    if SciMLBase.successful_retcode(sol)
        cumulativecases = modelcompartments(sol, 8)
        incidentcases = casespertimeblock(cumulativecases .* 5_450_000 .* p.proportiondetected)
        modeloutputs[i] = incidentcases
    else
        @warn "$(sol.retcode) with p=$p, u0=$u0"
        modeloutputs[i] = missings(353)
    end
end

medianoutput = Vector{Float64}(undef, 353)
lcioutput = Vector{Float64}(undef, 353)
ucioutput = Vector{Float64}(undef, 353)
for t ∈ 1:353 
    lc, me, uc = quantile(
        skipmissing([ modeloutputs[i][t] for i ∈ axes(chaindf, 1) ]), 
        [ 0.025, 0.5, 0.975 ]
    )
    medianoutput[t] = me
    lcioutput[t] = lc 
    ucioutput[t] = uc 
end

fig = Figure(; size=( 500, 500 ))
ga = GridLayout(fig[1, 1])
ax = Axis(ga[1, 1]; xticks=2017:2:2023,)# yticks=0:200:600)
lines!(ax, data.Date, medianoutput; color=COLOURVECTOR[1], linewidth=1,)
band!(ax, data.Date, lcioutput, ucioutput; color=( COLOURVECTOR[1], 0.5 ))
scatter!(ax, data.Date, data.Cases; color=:black, markersize=3)


fig
=#

println("Completed with n_rounds=$n_rounds, optimizationsolvermaxiters=$optimizationsolvermaxiters")






chaindict2 = load(datadir("sims", "chaindict_nrounds_25.jld2"))
chaindf = DataFrame(chaindict2["chain"])
