
using DrWatson

@quickactivate :ImmuneBoostingODEs

using DataFrames
using DifferentialEquations
using Memoization
using Optimization
using OptimizationOptimJL
using OrdinaryDiffEq
using Random
using Turing
using Zygote

testrun = true 

if length(ARGS) == 2
    const r0 = parse(Float64, ARGS[1])
    const omega = parse(Float64, ARGS[2])
    n_rounds = 2000
    optimizationsolvermaxiters = 1_000_000
else
    const r0 = 2.0
    const omega = 0.5
    n_rounds = testrun ? 25 : 2000
    optimizationsolvermaxiters = testrun ? 10_000 : 1_000_000
end

println("Starting with r0=$r0, omega=$omega, n_rounds=$n_rounds, optimizationsolvermaxiters=$optimizationsolvermaxiters")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load the data 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

include("rsvsetup.jl")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Fitting parameters 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function transformedsirns_fixedr0omega!(du, u, p, t)
    newparms = transformparameters(p; r0, omega)
    sirns!(du, u, newparms, t)
end

function loss(
    parms;  # a vector containing the following in order:
        # logitβ1, 
        # ϕ, 
        # logγ,
        # logψ, 
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
    r0,
    omega,
) 
    if parms[2] < -π || parms[2] > π || parms[8] < 0
        return Inf
    end

    I0 = ImmuneBoostingODEs._logistic(parms[10])
    S0 = min(ImmuneBoostingODEs._logistic(parms[9]), 1 - I0)
    u0 = sirns_u0_transformedp(S0, I0; p=parms, r0, omega, equalrs, t0)
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
                        0.25 / (0.25 + incidentcases[t] * ImmuneBoostingODEs._logistic(parms[7]))
                    ),
                    data[t]
                )
            )
            for t ∈ eachindex(data)
        ]
    ) - abs2(ImmuneBoostingODEs._logistic(parms[6]) - parms[8])  # so that parms[8] has some influence on `loss` 
    return loss
end

function optimizesirns(
    data,
    parms;
    callback,
    saveat,
    equalrs=true,
    t0=1996.737,
    tspan=( saveat[1], last(saveat) ),
    mu=0.0087,
    alg=Vern9(; lazy=false),
    odesolverabstol=1e-15,
    odesolvermaxiters=5e7,
    optimizationsolvermaxiters=1e5,
    adtype=Optimization.AutoZygote(),
    lb=[ -5.3, -2, 2.8, -0.7, -2, 0.24, -6.6, 0, -4.3, -10.3 ],
    ub=[ 0, 2, 5.2, 1.1, 0.2, 4.2, -2.6, 1, 4.3, -1.7 ],
    nt=10,
    rt=0.975,
    r_expand=2.0,
    verbosity=3,
    r0,
    omega,
)
    I0 = ImmuneBoostingODEs._logistic(parms[9])
    S0 = min(ImmuneBoostingODEs._logistic(parms[8]), 1 - I0)
    u0 = sirns_u0_transformedp(S0, I0; p=parms, r0, omega, equalrs, t0)
    prob = ODEProblem(transformedsirns_fixedr0omega!, u0, tspan, parms)
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
            r0,
            omega,
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

adjustedparams(p) = [ p[1:7]; 4.0; p[9]; p[10] .+ 6 ]

initial_params = Vector{Vector{Float64}}(undef, 8)

Threads.@threads for i ∈ 1:8 
    Random.seed!(n_rounds + optimizationsolvermaxiters + round(Int, omega) + i)

    lb = [ -Inf, -0.8, 3.3, -4.0, -Inf, -Inf, -3.892, 0.0, -4.3, -10.3 ]
    ub = [ Inf, 0.8, 4.4, 3.0, Inf, Inf, -3.892, 1.0, 4.3, -1.7 ]
    
    if isodd(i)
        lb[1] = ub[1] = log(0.01 / 0.99)
    else
        lb[1] = ub[1] = log(0.25 / 0.75)
    end

    if i ∈ [ 1, 2, 5, 6 ]
        lb[5] = ub[5] = log(0.1 / 0.9)
    else
        lb[5] = ub[5] = log(0.9 / 0.1)
    end

    if i <= 4
        lb[6] = ub[6] = log(0.75 / 0.25)
    else
        lb[6] = ub[6] = log(0.95 / 0.05)
    end

    _add = (ub .- lb) .* 0.5
    ip = optimizesirns(
        data.Cases, lb .+ _add;
        callback=optimcbs,
        saveat, 
        optimizationsolvermaxiters, 
        lb, 
        ub, 
        verbosity=0, 
        r0, 
        omega,
    )
    initial_params[i] = adjustedparams(ip.minimizer)
    @info "initial_params[$i]=$(initial_params[i]), $(ip.retcode)"
end

chaindictinit = Dict(
    "initial_params" => initial_params,
    "n_rounds" => n_rounds,
    "optimizationsolvermaxiters" => optimizationsolvermaxiters,
    "r0" => r0,
    "omega" => omega,
)

safesave(
    datadir("sims", "chaindictinit_nrounds_$(n_rounds)_r0_$(r0)_omega_$omega.jld2"), 
    chaindictinit
)

tspan = ( saveat[1], last(saveat) )
initialp = SirnsParameters(
    r0 * (48.7 + 0.0087),  # β0::S
    0.1,  # β1::T
    0.0,  # ϕ::T
    48.7,  # γ::T
    0.0087,  # μ::Float64 
    1.0,  # ψ::T
    omega,  # ω::S
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

Random.seed!(n_rounds + optimizationsolvermaxiters + round(Int, omega))

chain = sample(
    fitmodel(data.Cases, prob, cbs, saveat; r0, omega),
    #Turing.NUTS(0.65; adtype=AutoReverseDiff(false)),
    Turing.NUTS(0.65),
    MCMCThreads(),
    n_rounds,
    4;
    #initial_params,
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
    "r0" => r0,
    "omega" => omega,
)

safesave(
    datadir("sims", "chaindict_nrounds_$(n_rounds)_r0_$(r0)_omega_$omega.jld2"), 
    chaindict
)
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

println("Completed with r0=$r0, omega=$omega, n_rounds=$n_rounds, optimizationsolvermaxiters=$optimizationsolvermaxiters")



chaindict15_02 = load(datadir("sims", "chaindict_nrounds_2000_r0_1.5_omega_0.2.jld2"))
chaindf15_02 = DataFrame(chaindict15_02["chain"])
plotchains(chaindf15_02)

chaindict15_05 = load(datadir("sims", "chaindict_nrounds_2000_r0_1.5_omega_0.5.jld2"))
chaindf15_05 = DataFrame(chaindict15_05["chain"])
plotchains(chaindf15_05)

chaindict15_1 = load(datadir("sims", "chaindict_nrounds_2000_r0_1.5_omega_1.0.jld2"))
chaindf15_1 = DataFrame(chaindict15_1["chain"])
plotchains(chaindf15_1)

chaindict15_2 = load(datadir("sims", "chaindict_nrounds_2000_r0_1.5_omega_2.0.jld2"))
chaindf15_2 = DataFrame(chaindict15_2["chain"])
plotchains(chaindf15_2)

chaindict2_05 = load(datadir("sims", "chaindict_nrounds_1000_r0_2.0_omega_0.5.jld2"))
chaindf2_05 = DataFrame(chaindict2_05["chain"])
plotchains(chaindf2_05)

chaindict2_1 = load(datadir("sims", "chaindict_nrounds_1000_r0_2.0_omega_1.0.jld2"))
chaindf2_1 = DataFrame(chaindict2_1["chain"])
plotchains(chaindf2_1)

chaindict2_2 = load(datadir("sims", "chaindict_nrounds_1000_r0_2.0_omega_2.0.jld2"))
chaindf2_2 = DataFrame(chaindict2_2["chain"])
plotchains(chaindf2_2)

chaindict3_05 = load(datadir("sims", "chaindict_nrounds_1000_r0_3.0_omega_0.5.jld2"))
chaindf3_05 = DataFrame(chaindict3_05["chain"])
plotchains(chaindf3_05)

chaindict3_1 = load(datadir("sims", "chaindict_nrounds_1000_r0_3.0_omega_1.0.jld2"))
chaindf3_1 = DataFrame(chaindict3_1["chain"])
plotchains(chaindf3_1)

chaindict3_2 = load(datadir("sims", "chaindict_nrounds_1000_r0_3.0_omega_2.0.jld2"))
chaindf3_2 = DataFrame(chaindict3_2["chain"])
plotchains(chaindf3_2)

chaindict5_05 = load(datadir("sims", "chaindict_nrounds_1000_r0_5.0_omega_0.5.jld2"))
chaindf5_05 = DataFrame(chaindict5_05["chain"])
plotchains(chaindf5_05)

chaindict5_1 = load(datadir("sims", "chaindict_nrounds_1000_r0_5.0_omega_1.0.jld2"))
chaindf5_1 = DataFrame(chaindict5_1["chain"])
plotchains(chaindf5_1)

chaindict5_2 = load(datadir("sims", "chaindict_nrounds_1000_r0_5.0_omega_2.0.jld2"))
chaindf5_2 = DataFrame(chaindict5_2["chain"])
plotchains(chaindf5_2)

chaindict9_05 = load(datadir("sims", "chaindict_nrounds_1000_r0_9.0_omega_0.5.jld2"))
chaindf9_05 = DataFrame(chaindict9_05["chain"])
plotchains(chaindf9_05)

chaindict9_1 = load(datadir("sims", "chaindict_nrounds_1000_r0_9.0_omega_1.0.jld2"))
chaindf9_1 = DataFrame(chaindict9_1["chain"])
plotchains(chaindf9_1)

chaindict9_2 = load(datadir("sims", "chaindict_nrounds_1000_r0_9.0_omega_2.0.jld2"))
chaindf9_2 = DataFrame(chaindict9_2["chain"])
plotchains(chaindf9_2)

